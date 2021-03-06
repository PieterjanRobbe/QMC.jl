# The QMC module for Julia
[![Build Status](https://travis-ci.org/PieterjanRobbe/QMC.jl.png)](https://travis-ci.org/PieterjanRobbe/QMC.jl)

This module provides an implementation of some standard Quasi-Monte Carlo (QMC) point set generators.
Currently, only rank-1 lattice rules and randomized rank-1 lattice rules in base 2 are supported.

The module generates quasi-random points in `s` dimensions that are distributed *better than random*.
These points are used in [Quasi-Monte Carlo methods](https://en.wikipedia.org/wiki/Quasi-Monte_Carlo_method) to speed up 
the convergence of traditional [Monte Carlo methods](https://en.wikipedia.org/wiki/Monte_Carlo_method).

## Rank-1 lattice rules

QMC methods are equal-weight cubature rules to approximate high-dimensional integrals over the unit cube `[0,1]^s`. For rank-1 lattice rules, the quadrature points are chosen as multiples of a deterministic *generating vector* `z`, mapped back to the unit cube by taking the `mod 1`. The typical best convergence that can be achieved using these lattice rules is `O(N^−1)`, where `N` is the number of points.

Note that our lattice rules are in fact (extensible) lattice sequences, since we apply a *radical inverse* reordening to every point in the sequence. By doing this, the total number of points doesn't need to be known in advance.

## Installation

```julia
] add https://github.com/PieterjanRobbe/QMC.jl
```

## Usage

Now you should be able to 

```julia
using QMC
```

A lattice sequence in `s = 2` dimensions can be constructed as

```julia
lat = LatSeq(2)
```

The first point of the sequence can be obtained as

```julia
next(lat)
```
It is a `Vector{Float64}` of length `s` containing

```
[0.0, 0.0]
```

We can also ask for multiple points at once, using the second argument, i.e.,

```julia
next(lat,10)
```

This gives the next 10 points of the sequence:

```
[0.5,0.5]      
[0.75,0.25]    
[0.25,0.75]    
[0.375,0.125]  
[0.875,0.625]  
[0.625,0.875]  
[0.125,0.375]  
[0.1875,0.0625]
[0.6875,0.5625]
[0.9375,0.3125]
```

We provide two different generating vectors `z`:

* the `250`-dimensional generating vector from Cools, R., Kuo, F. Y., and Nuyens, D., "Constructing embedded lattice rules for multivariate integration." *SIAM Journal on Scientific Computing* 28.6 (2006): 2162-2188
* a `3600`-dimensional generating vector with product weights `1/j^2` from [Kuo's website](http://web.maths.unsw.edu.au/~fkuo/lattice/)

Depending on the size of `s`, a suitable generating vector is used. If needed, the user can supply their own generating vector `z` (in `Unsigned` format) with

```julia
lat = LatSeq(z,s)
```
There is an optional argument `nmax` that gives the maximum number of points in the sequence. By default, this is `typemax(z[1])`. The module also offers randomized lattice rules (see [test/QMCtest.jl](test/QMCtest.jl) for a typical use).

```julia
lat = LatSeq(250)
randlat = RandWrapper(lat,16) # randomized lattice rule in 250 dimensions with 16 shifts
next(randlat)
```
This last command returns a `250x16` matrix containing the first `250`-dimensional point of the lattice rule shifted by 16 uniformly distributed random shifts. Similar as above, we can request `N` points at once using

```julia
points = next(randlat,N)
```

Any point set generator can be reset using

```julia
reset(lat) # next point of the sequnce will be [0., 0.]
```

The code also allows the use of (interlaced) polynomial lattice rules (IPLR). These lattice rules are based on generating matrices. The points can be generated as a digital sequence with a finite field `F` in base `b`. For practical reasons, `b = 2`. Given the generating matrix `C` of a polynomial lattice rule, or the generating matrix `B` of an interlaced polynomial lattice rule, we construct a sequence of quasi-random points. The point generators are able to construct higher-order point sets with a theoretical convergence rate of `O(N^(-alpha))`, where `alpha>0`.

We provide two default generating matrices:

* a simple `32`-dimensional generating matrix from Nuyens' [Magic point shop](https://people.cs.kuleuven.be/~dirk.nuyens/qmc-generators/) 
* a `1111`-dimensional `32`-bit generating matrix from the [Joe and Kuo paper in ACM Transactions om Mathematical Software](https://people.cs.kuleuven.be/~dirk.nuyens/qmc-generators/sobolmats/)


The use of these lattice rules is the same as for rank-1 lattice rules, i.e., we can initialize a generator and ask points with

```julia
dig = DigSeq(25)
next(dig)
```

The `RandWrapper` now does digital shifting under the hood, but it can be used in the same way as before.

```julia
dig = DigSeq(25)
randdig = RandWrapper(dig,8)
next(randdig)
```

Other generating vectors or generating matrices can be found or constructed using the following resources:

* [Magic point shop](https://people.cs.kuleuven.be/~dirk.nuyens/qmc-generators/)
* [QMC4PDE](https://people.cs.kuleuven.be/~dirk.nuyens/qmc4pde/)
* [Kuo's homepage](kuo generating vectors)

## A simple example

Let's set up a simple lattice sequence and investigate the projection of certain dimensions.
First, construct a lattice sequence and generate some points:

```julia
using PyPlot

lat = LatSeq(10) # 10-dimensional lattice sequence
points = next(lat,1024)
```
Next, we plot dimension `9` against dimension `4`:

```julia
points = hcat(points...)' # to "flat" array
figure(1), subplot(111,aspect="equal")
plot(points[:,9], points[:,4], "b.") # a good projection
```
![projection of dimension 9 versus dimension 4](figures/9_versus_4.png "projection of dimension 9 versus dimension 4")

This is a good projection, since the points distribute the space `[0,1]^s` evenly. Next, look at the projection of dimension `6` against dimension `2`:

```julia
figure(2), subplot(111,aspect="equal")
plot(points[:,6], points[:,2], "b.") # a bad projection
```
![projection of dimension 6 versus dimension 2](figures/6_versus_2.png "projection of dimension 9 versus dimension 4")

Now, we clearly see a pattern in the distribution of the lattice points. Other bad (but nice) projections are `1` versus `9`, `8` versus `5` or `10` versus `4`.
