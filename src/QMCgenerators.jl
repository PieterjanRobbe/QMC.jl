# QMC generator type in s dimensions
abstract QMCgenerator{s}

#
# RANDOMLY SHIFTED LATTICE RULES
#

# lattice sequence (32 bit)
type LatSeq{s,U<:Unsigned} <: QMCgenerator{s}
  	z::Vector{U} # generating vector
	nmax::U # max number of points
	k::U # current point number
	recipid::Float64
end

# constructor
function LatSeq{U<:Unsigned,N<:Integer}(z::Vector{U}, s::N; nmax::U = typemax(U))
 	(s > 0 && s <= length(z)) || error("Can only generate $(nmax) points in $(length(z)) dimensions. 
 		Supply your own generating vector to increase.")
 	t = length(bin(typemax(U)))
 	LatSeq{s,U}(z, nmax, zero(U), exp2(-t))
end

# default generating vectors
function LatSeq(s::Integer)
	if s <= 250
		LatSeq(CKN_250_20, s, nmax = convert(UInt32,2^20))
	elseif s <=3600
		LatSeq(K_3600_32, s)
	else
		error("Cannot generate lattice points beyond dimension 3600, supply your own generating vector.")
	end
end

# functions
ndims{s}(l::LatSeq{s}) = s

function show(io::IO, l::LatSeq)
    print(io, "$(ndims(l))-dimensional lattice sequence")
end

# get the k'th point of the lattice sequence
# a point has dimension (s,), where s is the dimenion of the lattice sequence
function getPoint{s,N<:Integer}(l::LatSeq{s}, k::N)
	k = convert(UInt32,k)
	@assert k < l.nmax
	phi_k = convert(Float64,reversebits(k $ k >>> 1)) 
	x = zeros(s)
	for i in 1:s
		@inbounds x[i] = phi_k*l.recipid*l.z[i]
	end
	return mod(x,1)
end

#
# INTERLACED POLYNOMIAL LATTICE RULES
#

# interlaced polynomial lattice sequence
type DigSeq{s,U<:Unsigned} <: QMCgenerator{s}
  	C::Matrix{U} # generating matrix
	nmax::U # max number of points
	k::U # current point number
	cur::Vector{U} # previous point as integer
	recipid::Float64
end

# constructor
function DigSeq{U<:Unsigned,N<:Integer}(C::Matrix{U}, s::N; nmax::U = typemax(U))
 	(s > 0 && s <= size(C,1)) || error("Can only generate $(nmax) points in $(size(C,1)) dimensions. 
 		Supply your own generating matrix to increase.")
 	C = reversebits(C)#[reversebits(C[i,j]) for i in 1:size(C,1), j in 1:size(C,2)]
 	t = maximum([ length(bin(C[1,i])) for i in 1:size(C,2) ])
 	DigSeq{s,U}(C, nmax, zero(U), zeros(U,s), exp2(-t))
end

# default generating matrices
function DigSeq(s::Integer)
	if s <= 32
		DigSeq(N_32_32, s)
	elseif s <=1111
		DigSeq(JK_1111_32, s)
	else
		error("Cannot generate lattice points beyond dimension 1111, supply your own generating matrix.")
	end
end

# functions
ndims{s}(l::DigSeq{s}) = s

function show(io::IO, l::DigSeq)
    print(io, "$(ndims(l))-dimensional digital sequence")
end

# get the k'th point of the digital sequence
# a point has dimension (s,), where s is the dimenion of the digital sequence
function getPoint{s,U<:Unsigned}(d::DigSeq{s}, k::U)
	k = convert(U,k)
	k > 0 || return zeros(s) # return zero if k == 0
	@assert k < d.nmax
	c = length(bin(((k $ (k-1)) + 1) >> 1))
	x = zeros(s)
	for i in 1:s
		d.cur[i] = d.cur[i] $ d.C[i,c]
		x[i] = d.recipid * d.cur[i]
	end
	return x
end

#
# RANDOMIZATION
#

# random shift wrapper
type RandWrapper{s,q,T}
	generator::QMCgenerator{s} # underlying QMC generator
	shifts::Matrix{T} # random shifts (s by q matrix)
	shift_fcn::Function # function to shift this lattice sequence
end

# constructors
function RandWrapper{s,U,N}(generator::LatSeq{s,U}, q::N)
	shifts = rand(s,q)
	RandWrapper{s,q,Float64}(generator, shifts, lat_shift)
end

function RandWrapper{s,U,N}(generator::DigSeq{s,U}, q::N)
	shifts = convert(Array{U},floor(rand(s,q)*M))
	RandWrapper{s,q,U}(generator, shifts, dig_shift)
end

# functions
ndims{s,q}(r::RandWrapper{s,q}) = s

nshifts{s,q}(r::RandWrapper{s,q}) = q

function show(io::IO, r::RandWrapper)
    print(io, "$(ndims(r))-dimensional randomized sequence with $(nshifts(r)) shifts")
end

# get the k'th point of the randomized QMC generator
# a point has dimension (s,q) with s the dimension of the QMC generator and q the number of shifts
function getPoint{s,q,N<:Integer}(r::RandWrapper{s,q}, k::N)
	x = getPoint(r.generator, k)
	return r.shift_fcn(x,r.shifts)
end

# function to shift the lattice sequence
function lat_shift{T}(x::Vector{T},shifts::Matrix{T})
	return mod(x .+ shifts,1)
end

# function to shift the digital sequence
function dig_shift{T,U}(x::Vector{T},shifts::Matrix{U})
	Ps = (shifts $ repmat(convert(Array{UInt32},floor(M*x)),1,size(shifts,2))) * recipid;
end

#
# GENERAL FUNCTIONS
#

# get the next point of the QMC generator
# a point has dimension (s,), where s is the dimenion of the QMC generator
function next(generator::QMCgenerator)
	x = getPoint(generator, generator.k)
	generator.k += 1
	return x
end

# get the next n points of the QMC generator
# a point has dimension (s,), where s is the dimenion of the lattice sequence
# all points are returned in a vector of size (n,)
function next{N<:Integer}(generator::QMCgenerator, n::N)
	x = Vector{Float64}[] # turns out to be faster and more memeff than preallocating
	for i in 1:n
		@inbounds push!(x, next(generator))
	end
	return x
end

# get the next point of the randomized wrapper
# a point has dimension (s,q), where s is the dimenion of the QMC generator and q the number of shifts
function next(r::RandWrapper)
	x = getPoint(r, r.generator.k)
	r.generator.k += 1
	return x
end

# get the next n points of the QMC generator
# a point has dimension (s,q), where s is the dimenion of the lattice sequence and q the number of shifts
# all points are returned in a vector of size (n,)
function next{N<:Integer}(r::RandWrapper, n::N)
	x = Array{Float64,2}[]
	for i in 1:n
		@inbounds push!(x, next(r))
	end
	return x
end

# reset the counter of the QMC generator
function reset(g::QMCgenerator)
	g.k = 0
end

# reset the counter of the QMC generator
function reset(r::RandWrapper)
	reset(r.generator)
end

# make the generator iterable
start(g::QMCgenerator) = g

next(g::QMCgenerator, g_::QMCgenerator) = (next(g), g_)

done(g::QMCgenerator, g_::QMCgenerator) = false # infinite loop

start(r::RandWrapper) = r

next(r::RandWrapper, r_::RandWrapper) = (next(r), r_)

done(r::RandWrapper, r_::RandWrapper) = false # infinite loop

# reverse bits of UInt type
reversebits{U<:Unsigned}(n::U)=parse(U,reverse(bits(n)),2)
@vectorize_1arg Unsigned reversebits

# constant
const recipid = exp2(-32)
const M = 2^32