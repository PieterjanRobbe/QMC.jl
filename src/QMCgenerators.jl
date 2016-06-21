# QMC generator type
abstract QMCgenerator

# lattice sequence
type LatSeq{s} <: QMCgenerator
  	z::Vector{UInt32} # generating vector
	nmax::UInt32 # max number of points
	k::UInt32 # current point number
end

# randomly shifted lattice sequence
type RandLatSeq{s,q} <: QMCgenerator
	latSeq::LatSeq{s} # underlying lattice rule
	shifts::Matrix{Float64} # random shifts (s by q matrix)
	k::UInt32 # current point number (duplicate)
end

# constructors
function LatSeq(s::Integer)
	(s > 0 && s <= 3600) || error("Can only generate 2^20 points in 3600 dimensions. 
		Supply your own generating vector to increase.")	
	if s <= 250
		LatSeq{s}(CKN_250_20, 2^20, 0)
	else
		LatSeq{s}(K_3600_20, 2^20, 0)
	end
end

function LatSeq{N<:Integer}(z::Vector{N}, s::Integer, nmax::Integer)
 	(s > 0 && s <= length(z)) || error("Can only generate points $(nmax) in $(length(z)) dimensions. 
 		Supply your own generating vector to increase.")
 	LatSeq{s}(z, nmax, 0)
end

function RandLatSeq{N<:Integer}(q::N)
	l = LatSeq(250)
	shifts = rand(250,q)
	RandLatSeq{250,q}(l, shifts, 0)
end

function RandLatSeq{N<:Integer}(s::N, q::N)
	l = LatSeq(s)
	shifts = rand(s,q)
	RandLatSeq{s,q}(l, shifts, 0)
end

function RandLatSeq{s,N<:Integer}(l::LatSeq{s}, q::N)
	shifts = rand(s,q)
	RandLatSeq{s,q}(l, shifts, 0)
end

# functions
ndims{s}(l::LatSeq{s}) = s

nshifts{s}(l::LatSeq{s}) = 1

ndims{s,q}(r::RandLatSeq{s,q}) = s

nshifts{s,q}(r::RandLatSeq{s,q}) = q

function show(io::IO, l::LatSeq)
    print(io, "$(ndims(l))-dimensional lattice sequence")
end

function show(io::IO, r::RandLatSeq)
    print(io, "$(ndims(r))-dimensional randomized lattice sequence with $(nshifts(r)) shifts")
end

# get the (unmod'ed) k'th point of the lattice sequence
# a point has dimension (s,), where s is the dimenion of the lattice sequence
function getPoint{s}(l::LatSeq{s}, k::UInt32)
	@assert k < l.nmax
	phi_k = convert(Float64,bitreverse32(k $ k >>> 1)) 
	x = zeros(s)
	for i in 1:s
		@inbounds x[i] = phi_k*recipid*l.z[i]
	end
	return x
end

# get the (unmod'ed) k'th point of the randomized lattice sequence
# a point has dimension (s,q) with s the dimension of the lattice sequence and
# q is the number of shifts
function getPoint{s,q}(r::RandLatSeq{s,q}, k::UInt32)
	x = getPoint(r.latSeq, k)
	return x .+ r.shifts
end

# general functions that apply to all QMC generators

# get the next point of the QMC generator
function next!(g::QMCgenerator)
	x = mod(getPoint(g, g.k),1)
	g.k += 1
	return x
end

# get the next n points of the QMC generator
function next!{N<:Integer}(g::QMCgenerator, n::N)
	x = Array{Float64,nshifts(g)}[] # turns out to be faster and more memeff than preallocating
	for i in 1:n
		@inbounds push!(x, next!(g))
	end
	return x
end

# make the lattice sequence iterable
start(g::QMCgenerator) = g

next(g::QMCgenerator, g_::QMCgenerator) = (next!(g), g_)

done(g::QMCgenerator, g_::QMCgenerator) = false # infinite loop

# reset the counter of the QMC generator
function reset!(g::QMCgenerator)
	g.k = 0
end

# reverse all bits of unsigned 32-bit integer
# see http://graphics.stanford.edu/~seander/bithacks.html#ReverseParallel
function bitreverse32(v::UInt32)
	v = ((v >> 1) & 0x55555555) | ((v & 0x55555555) << 1)  # swap odd and even bits
	v = ((v >> 2) & 0x33333333) | ((v & 0x33333333) << 2)  # swap consecutive pairs
	v = ((v >> 4) & 0x0f0f0f0f) | ((v & 0x0f0f0f0f) << 4)  # swap nibbles 
	v = ((v >> 8) & 0x00ff00ff) | ((v & 0x00ff00ff) << 8)  # swap bytes
	v = ( v >> 16             ) | ( v               << 16) # swap 2-byte long pairs
	return v
end

# constant
const recipid = exp2(-32)
