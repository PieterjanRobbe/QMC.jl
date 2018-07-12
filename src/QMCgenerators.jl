## QMCgenerators.jl: QMC point set generators

"""
QMCgenerator{s}

Supertype for a QMC generator that generates points in `s` dimensions.
"""
abstract type QMCgenerator{s} end

## Randomly shifted lattice rules ##
"""
LatSeq{s,U,V}

`s`-dimensional lattice sequence based on the unsigned integer type `U`. `V` is the type of the generating vector, i.e., `V=Vector{U}`.
"""
mutable struct LatSeq{s,U<:Unsigned,V<:AbstractVector} <: QMCgenerator{s}
    z::V # generating vector 
    nmax::U # max number of points
    k::U # current point number

    recipid::Float64 # internal
end

"""
LatSeq(z, s, nmax)

Create a lattice sequence in `s` dimensions from the generating vector `z` that can generate at most `nmax` points.

# Examples

```jldoctest
julia> LatSeq([0x00000001,0x00000022],2,0x00000037) # fibonacci lattice with 55 points
2-dimensional lattice sequence
```
"""
function LatSeq(z::AbstractVector, s::N where N<:Integer, nmax::U = typemax(U)) where U<:Unsigned
    s <= 0 && throw(BoundsError("s must be positive"))
    s > length(z) && throw(BoundsError("Can only generate $(nmax) points in $(length(z)) dimensions. 
                                       Supply a different generating vector to increase."))
    eltype(z) == U || throw(TypeError("Generating vector element type and type of `nmax` do not agree."))
    t = length(bin(typemax(U)))
    LatSeq{s,eltype(z),typeof(z)}(z, nmax, zero(U), exp2(-t))
end

# using default generating vectors
"""
LatSeq(s)

Create a lattice sequence in `s` dimensions using a default generating vector.

# Examples

```jldoctest
julia> LatSeq(112)
112-dimensional lattice sequence
```
"""
function LatSeq(s::Integer)
    if s <= 250
        LatSeq(CKN_250_20, s, convert(UInt32,2^20))
    elseif s <= 3600
        LatSeq(K_3600_32, s, typemax(UInt32))
    else
        throw(BoundsError("Cannot generate lattice points beyond dimension 3600, supply your own generating vector."))
    end
end

"""
getPoint(lat, k)

Get the `k`-th point of the lattice sequence `lat`. A point is a vector of length `s`, where `s` is the dimension of the lattice sequence.

# Examples

```jldoctest
julia> lat = LatSeq(8)
8-dimensional lattice sequence

julia> getPoint(lat,10)
8-element Array{Float64,1}:
0.9375
0.3125
0.8125
0.9375
0.4375
0.5625
0.4375
0.8125
```
"""
function getPoint(lat::LatSeq{s}, k::N where N<:Integer) where s
    k = convert(UInt32,k)
    @assert k < lat.nmax
    phi_k = convert(Float64,reversebits(xor(k, k >>> 1))) 
    x = zeros(s)
    for i in 1:s
        @inbounds x[i] = phi_k*lat.recipid*lat.z[i]
    end
    return mod.(x,1)
end

## Interlaced polynomial lattice rules ##
"""
DigSeq{s,U,M}

`s`-dimensional lattice sequence based on the unsigned integer type `U`. `M` is the type of the generating matrix, i.e., `V=Matrix{U}`.
"""
mutable struct DigSeq{s,U<:Unsigned,M<:AbstractMatrix} <: QMCgenerator{s}
    C::M # generating matrix
    nmax::U # max number of points
    k::U # current point number

    recipid::Float64 # internal
end

"""
DigSeq(C, s, nmax)

Create a digital sequence in `s` dimensions from the generating matrix `C` that can generate at most `nmax` points.

"""
function DigSeq(C::AbstractMatrix, s::N where N<:Integer, nmax::U = typemax(U)) where U<:Unsigned
    s <= 0 && throw(BoundsError("s must be positive"))
    s > size(C,1) && throw(BoundsError("Can only generate $(nmax) points in $(size(C,1)) dimensions. 
                                       Supply a different generating vector to increase."))
    eltype(C) == U || throw(TypeError("Generating vector element type and type of `nmax` do not agree."))
    C = reversebits.(C)
    t = maximum([ length(bin(C[1,i])) for i in 1:size(C,2) ])
    DigSeq{s,eltype(C),typeof(C)}(C, nmax, zero(U), exp2(-t))
end

# using default generating matrices
"""
DigSeq(s)

Create a digital sequence in `s` dimensions using a default generating matrix.

# Examples

```jldoctest
julia> dig = DigSeq(8)
8-dimensional digital sequence 
```
"""
function DigSeq(s::Integer)
    if s <= 32
        DigSeq(N_32_32, s, typemax(UInt32))
    elseif s <=1111
        DigSeq(JK_1111_32, s, typemax(UInt32))
    else
        throw(BoundsError("Cannot generate lattice points beyond dimension 1111, supply your own generating matrix."))
    end
end

"""
getPoint(dig, k)

Get the `k`-th point of the digital sequence `dig`. A point is a vector of length `s`, where `s` is the dimension of the digital sequence.

# Examples
```jldoctest
julia> dig = DigSeq(8)
8-dimensional digital sequence

julia> getPoint(dig,10)
8-element Array{Float64,1}:
0.994143
0.712787
0.971628
0.689584
0.532585
0.664304
0.837412
0.64264 
```
"""
function getPoint(d::DigSeq{s,U}, k::N where N<:Integer) where {s,U<:Unsigned}
    k = convert(U,k)
    k > 0 || return zeros(s) # return zero if k == 0
    @assert k < d.nmax
    cur = zeros(U,s)
    for i in 1:length(bits(k))
        if ( k & (1 << (i - 1) ) ) != 0
            for j in 1:s
                cur[j] ⊻= d.C[j,i]
            end
        end
    end
    x = zeros(s)
    for i in 1:s
        x[i] = d.recipid * cur[i]
    end
    return x
end

## Randomisation ##

"""
RandWrapper{s,q,T,M}

A randomized version of the underlying QMC generator in `s` dimensions using `q` shifts. The shifts are a matrix `M` with element type `T`.
"""
mutable struct RandWrapper{s,q,T,M<:AbstractMatrix}
    generator::QMCgenerator{s} # underlying QMC generator
    shifts::M # random shifts (s by q matrix)
end

"""
RandWrapper(lat, q)

Randomize the lattice sequence `lat` (or general QMCgenerator) using `q` shifts.

# Examples
```jldoctest
julia> lat = LatSeq(8)
8-dimensional lattice sequence

julia> ran = RandWrapper(lat,16)
8-dimensional randomized sequence with 16 shifts
```
"""
function RandWrapper(generator::QMCgenerator{s}, q::N where N<:Integer) where s
    shifts = rand(s,q)
    RandWrapper{s,q,eltype(shifts),typeof(shifts)}(generator, shifts)
end

"""
RandWrapper(dig, q)

Randomize the digital sequence `dig` using `q` shifts.

# Examples
```jldoctest
julia> dig = DigSeq(8)
8-dimensional digital sequence

julia> ran = RandWrapper(dig,16)
8-dimensional randomized sequence with 16 shifts
```
"""
function RandWrapper(generator::DigSeq{s,U,M}, q::N where N<:Integer) where {s,U,M}
    shifts = convert(Array{U},floor.(rand(s,q)*typemax(U)))
    RandWrapper{s,q,eltype(shifts),typeof(shifts)}(generator, shifts)
end

"""
getPoint(ran, k)

Get the `k`-th point of the randomized lattice sequence `rand`. A point is a matrix of size `s`+`q`, where `s` is the dimension of the sequence, and `q` is the number of shifts.

# Examples
```jldoctest
julia> lat = LatSeq(2)
2-dimensional lattice sequence

julia> ran = RandWrapper(lat,4)
2-dimensional randomized sequence with 4 shifts

julia> getPoint(ran, 10)
2×4 Array{Float64,2}:
0.516178  0.834656  0.935372  0.757624
0.711024  0.345201  0.485886  0.672963
```
"""
function getPoint(r::RandWrapper{s,q,Float64,M}, k::N where N<:Integer) where {s,q,M}
    x = getPoint(r.generator, k)
    return mod.(x .+ r.shifts,1)
end

"""
getPoint(ran, k)

Get the `k`-th point of the randomized digital sequence `rand`. A point is a matrix of size `s`-by-`q`, where `s` is the dimension of the sequence, and `q` is the number of shifts.

# Examples
```jldoctest
julia> dig = DigSeq(2)
2-dimensional digital sequence

julia> ran = RandWrapper(dig,4)
2-dimensional randomized sequence with 4 shifts

julia> getPoint(ran, 10)
2×4 Array{Float64,2}:
0.562022  0.734445  0.923817  0.757195
0.334069  0.212996  0.896181  0.214044
```
"""
function getPoint(r::RandWrapper{s,q,U,M}, k::N where N<:Integer) where {s,q,U<:Unsigned,M}
    x = getPoint(r.generator, k)
    return xor.(r.shifts, repmat(convert(Array{UInt32},floor.(typemax(U)*x)),1,size(r.shifts,2))) * r.generator.recipid
end

## General behavior ##

# get number of dimensions
ndims(g::QMCgenerator{s}) where s = s

ndims(r::RandWrapper{s,q}) where {s,q} = s

nshifts(r::RandWrapper{s,q}) where {s,q} = q

# show
function show(io::IO, l::DigSeq)
    print(io, "$(ndims(l))-dimensional digital sequence")
end

function show(io::IO, r::RandWrapper)
    print(io, "$(ndims(r))-dimensional randomized sequence with $(nshifts(r)) shifts")
end

function show(io::IO, l::LatSeq)
    print(io, "$(ndims(l))-dimensional lattice sequence")
end

"""
next(gen)

Get the next point of the sequence. A point is a vector of length `s`, where `s` is the dimension of the sequence.
"""
function next(g::QMCgenerator)
    x = getPoint(g, g.k)
    g.k += 1
    return x
end

"""
next(gen, n)

Get the next n points of the sequence. A point is a vector of length `s`, where `s` is the dimension of the sequence. All points are returned in a vector of size (n,).
"""
function next(g::QMCgenerator, n::N where N<:Integer)
    x = Vector{Float64}[] # turns out to be faster and more memeff than preallocating
    for i in 1:n
        @inbounds push!(x, next(g))
    end
    return x
end

"""
next(ran)

Get the next point of the randomized sequence. A point is a matrix of size `s`-by-`q`, where `s` is the dimension of the sequence and `q` is the number of shifts.
"""
function next(r::RandWrapper)
    x = getPoint(r, r.generator.k)
    r.generator.k += 1
    return x
end

"""
next(ran, n)

Get the next n points of the randomized sequence. A point is a matrix of size `s`-by-`q`, where `s` is the dimension of the sequence and `q` is the number of shifts. All points are returned in a vector of size (n,).
"""
function next(r::RandWrapper, n::N where N<:Integer)
    x = Array{Float64,2}[]
    for i in 1:n
        @inbounds push!(x, next(r))
    end
    return x
end

"""
reset(gen)

Reset the counter of the QMC generator.
"""
function reset(g::QMCgenerator)
    g.k = 0
end

"""
reset(ran)

Reset the counter of the randomized QMC generator.
"""
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
reversebits(n::U) where {U<:Unsigned}=parse(U,reverse(bits(n)),2)
