# QMC generator type in s dimensions
abstract QMCgenerator{s}

#
# RANDOMLY SHIFTED LATTICE RULES
#

# lattice sequence (32 bit)
type LatSeq{s,U<:Unsigned,V<:AbstractVector} <: QMCgenerator{s}
  	z::V # generating vector (vector of U's, see http://docs.julialang.org/en/release-0.4/manual/performance-tips/#avoid-fields-with-abstract-containers)
	nmax::U # max number of points
	k::U # current point number
	recipid::Float64

	LatSeq(z::AbstractVector{U},nmax::U,k::U,recipid::Float64) = new(z,nmax,k,recipid)
end

# constructor
function LatSeq{U<:Unsigned,N<:Integer}(z::AbstractVector, s::N, nmax::U = typemax(U))
 	(s > 0 && s <= length(z)) || error("Can only generate $(nmax) points in $(length(z)) dimensions. 
 		Supply your own generating vector to increase.")
 	t = length(bin(typemax(U)))
 	LatSeq{s,eltype(z),typeof(z)}(z, nmax, zero(U), exp2(-t)) :: LatSeq{s,U}
end

# default generating vectors
function LatSeq(s::Integer)
	(s > 0 && s <= 3600) || error("Cannot generate lattice points beyond dimension 3600, supply your own generating vector.")
	if s <= 250
		LatSeq(CKN_250_20, s, convert(UInt32,2^20))
	else
		LatSeq(K_3600_32, s, typemax(UInt32))
	end
end

# get the k'th point of the lattice sequence
# a point has dimension (s,), where s is the dimenion of the lattice sequence
function getPoint{s,N<:Integer}(l::LatSeq{s}, k::N)
	k = convert(UInt32,k)
	#@assert k < l.nmax
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
type DigSeq{s,U<:Unsigned,M<:AbstractMatrix} <: QMCgenerator{s}
  	C::M # generating matrix
	nmax::U # max number of points
	k::U # current point number
	recipid::Float64

	DigSeq(C::AbstractMatrix{U},nmax::U,k::U,recipid::Float64) = new(C,nmax,k,recipid)
end

# constructor
function DigSeq{U<:Unsigned,N<:Integer}(C::AbstractMatrix, s::N, nmax::U = typemax(U))
 	(s > 0 && s <= size(C,1)) || error("Can only generate $(nmax) points in $(size(C,1)) dimensions. 
 		Supply your own generating matrix to increase.")
 	C = reversebits(C)
 	t = maximum([ length(bin(C[1,i])) for i in 1:size(C,2) ])
 	DigSeq{s,eltype(C),typeof(C)}(C, nmax, zero(U), exp2(-t))
end

# default generating matrices
function DigSeq(s::Integer)
	if s <= 32
		DigSeq(N_32_32, s, typemax(UInt32))
	elseif s <=1111
		DigSeq(JK_1111_32, s, typemax(UInt32))
	else
		error("Cannot generate lattice points beyond dimension 1111, supply your own generating matrix.")
	end
end

# get the k'th point of the digital sequence
# a point has dimension (s,), where s is the dimenion of the digital sequence
function getPoint{s,N<:Integer,U<:Unsigned}(d::DigSeq{s,U}, k::N)
	k = convert(U,k)
	k > 0 || return zeros(s) # return zero if k == 0
	@assert k < d.nmax
	cur = zeros(U,s)
	for i in 1:length(bits(k))
		if ( k & (1 << (i - 1) ) ) != 0
			for j in 1:s
				cur[j] $= d.C[j,i]
			end
		end
	end
	x = zeros(s)
	for i in 1:s
		x[i] = d.recipid * cur[i]
	end
	return x
end

#
# RANDOMIZATION
#

# random shift wrapper
type RandWrapper{s,q,T,M<:AbstractMatrix}
	generator::QMCgenerator{s} # underlying QMC generator
	shifts::M # random shifts (s by q matrix)

	RandWrapper(generator::QMCgenerator{s},shifts::AbstractMatrix{T}) = new(generator,shifts)
end

# constructors
function RandWrapper{s,U,V,N}(generator::LatSeq{s,U,V}, q::N)
	shifts = rand(s,q)
	RandWrapper{s,q,eltype(shifts),typeof(shifts)}(generator, shifts)
end

function RandWrapper{s,U,M,N}(generator::DigSeq{s,U,M}, q::N)
	shifts = convert(Array{U},floor(rand(s,q)*typemax(U)))
	RandWrapper{s,q,eltype(shifts),typeof(shifts)}(generator, shifts)
end

# get the k'th point of the randomized lattice sequence
# a point has dimension (s,q) with s the dimension of the QMC generator and q the number of shifts
function getPoint{s,q,M,N<:Integer}(r::RandWrapper{s,q,Float64,M}, k::N)
	x = getPoint(r.generator, k)
	return mod(x .+ r.shifts,1)
end

# get the k'th point of the randomized digital sequence
# a point has dimension (s,q) with s the dimension of the QMC generator and q the number of shifts
function getPoint{s,q,U<:Unsigned,M,N<:Integer}(r::RandWrapper{s,q,U,M}, k::N)
	x = getPoint(r.generator, k)
	return (r.shifts $ repmat(convert(Array{UInt32},floor(typemax(U)*x)),1,size(r.shifts,2))) * r.generator.recipid;
end

#
# GENERAL FUNCTIONS
#

# get number of dimensions
ndims{s}(g::QMCgenerator{s}) = s

ndims{s,q}(r::RandWrapper{s,q}) = s

nshifts{s,q}(r::RandWrapper{s,q}) = q# functions

function show(io::IO, l::DigSeq)
    print(io, "$(ndims(l))-dimensional digital sequence")
end

function show(io::IO, r::RandWrapper)
    print(io, "$(ndims(r))-dimensional randomized sequence with $(nshifts(r)) shifts")
end

function show(io::IO, l::LatSeq)
    print(io, "$(ndims(l))-dimensional lattice sequence")
end

# get the next point of the QMC generator
# a point has dimension (s,), where s is the dimenion of the QMC generator
function next(g::QMCgenerator)
	x = getPoint(g, g.k)
	g.k += 1
	return x
end

# get the next n points of the QMC generator
# a point has dimension (s,), where s is the dimenion of the lattice sequence
# all points are returned in a vector of size (n,)
function next{N<:Integer}(g::QMCgenerator, n::N)
	x = Vector{Float64}[] # turns out to be faster and more memeff than preallocating
	for i in 1:n
		@inbounds push!(x, next(g))
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

@vectorize_1arg Unsigned reversebits # awesome!