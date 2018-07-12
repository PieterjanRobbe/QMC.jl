using QMC
using Test
using SpecialFunctions
using Random
using StatsBase

srand(1208)

function approx_pi_qmc(n::Integer)
	lat = LatSeq(2)
	count = 0
	for i in 1:n
		x = next(lat)
		if x[1]*x[1] + x[2]*x[2] < 1
			count += 1
		end
	end
	return 4*count/n
end

function approx_pi_qmc_2(n::Integer)
	dig = DigSeq(2)
	count = 0
	for i in 1:n
		x = next(dig)
		if x[1]*x[1] + x[2]*x[2] < 1
			count += 1
		end
	end
	return 4*count/n
end

function norminv(x::AbstractArray{T}) where T
	return sqrt(2)*erfinv.(2*x .-1)
end


println("============================")
println("=    TEST CONSTRUCTORS     =")
println("============================")
d = [2 250 2500] # valid
n = [100 5 10 50]
for i in 1:length(d)
	println("Constructing lattice sequence in $(d[i]) dimensions...")
	l = LatSeq(d[i])
	println("Taking $(n[i]) points of the lattice sequence...")
	x = next(l,n[i])
	@test(size(x)==(n[i],))
end
d = [2 250 1111] # valid
n = [100 5 10 50]
for i in 1:length(d)
	println("Constructing digital sequence in $(d[i]) dimensions...")
	s = DigSeq(d[i])
	println("Taking $(n[i]) points of the lattice sequence...")
	x = next(s,n[i])
	@test(size(x)==(n[i],))
end
d = [-1 0 3601] # invalid
for i in 1:length(d)
	println("Constructing lattice sequence in $(d[i]) dimensions...")
	@test_throws BoundsError LatSeq(d[i])
end
d = [-1 0 1112] # invalid
for i in 1:length(d)
	println("Constructing digital sequence in $(d[i]) dimensions...")
	@test_throws BoundsError DigSeq(d[i])
end
d = [2 150] # valid
q = [8 16]
for i in 1:length(d)
	println("Constructing randomized lattice sequence in $(d[i]) dimensions with $(q[i]) shifts...")
	l = LatSeq(d[i])
	r = RandWrapper(l,q[i])
	println("Taking a point of the randomized lattice sequence...")
	x = next(r)
	@test(size(x)==(d[i],q[i]))
end
d = [2 150] # valid
q = [8 16]
for i in 1:length(d)
	println("Constructing randomized lattice sequence in $(d[i]) dimensions with $(q[i]) shifts...")
	s = DigSeq(d[i])
	r = RandWrapper(s,q[i])
	println("Taking a point of the randomized digital sequence...")
	x = next(r)
	@test(size(x)==(d[i],q[i]))
end

println("============================")
println("=   APPROXIMATION OF PI    =")
println("============================")
n = 2 .^(2:16)
v = [ 4.000000000000000, 3.000000000000000, 3.250000000000000, 3.250000000000000, 3.187500000000000,
      3.125000000000000, 3.140625000000000, 3.140625000000000, 3.136718750000000, 3.140625000000000,
      3.142578125000000, 3.144042968750000, 3.140869140625000, 3.141479492187500, 3.141540527343750 ]
for i in 1:length(n)
	println("Computing integral with $(n[i]) points...")
	@test approx_pi_qmc(n[i]) ≈v[i] atol=1e-4
end

println("=============================")
println("= APPROXIMATION OF PI (DIG) =")
println("=============================")
n = 2 .^(2:16)
v = [ 2.000000000000000, 2.000000000000000, 2.500000000000000, 2.500000000000000, 2.562500000000000,
	  2.531250000000000, 2.968750000000000, 3.085937500000000, 3.121093750000000, 3.128906250000000,
	  3.142578125000000, 3.141601562500000, 3.142333984375000, 3.143310546875000, 3.14190673828125 ]
for i in 1:length(n)
	println("Computing integral with $(n[i]) points...")
	@test approx_pi_qmc_2(n[i]) ≈v[i] atol=1e-4
end

println("============================")
println("=     COMPUTE INTEGRAL     =")
println("============================")
# see Keister, Bradley D. "Multidimensional quadrature algorithms." Computers in Physics 10.2 (1996): 119-128.
d = [9, 25, 60, 80, 100]
exact = [-71.633234291 -1.356914e6 4.89052986e14 6.78878724e19 4.57024396e24]
for i in 1:length(d)
	println("Computing cos-integral in $(d[i]) dimensions...")
	N = 0
	I = 0.
	l = LatSeq(d[i])
	skip = next(l)
	for x in l
		I = I + cos(sqrt(sum(norminv(x).^2/2)))
		N = N + 1
		if abs(pi^(d[i]/2)/N*I-exact[i])/abs(exact[i]) < 1e-3
			break
		end
	end
	println("  Converged after $(N) points!")
end

println("============================")
println("=  COMPUTE INTEGRAL (DIG)  =")
println("============================")
d = [9, 25, 60, 80, 100]
exact = [-71.633234291 -1.356914e6 4.89052986e14 6.78878724e19 4.57024396e24]
for i in 1:length(d)
	println("Computing cos-integral in $(d[i]) dimensions...")
	N = 0
	I = 0.
	l = DigSeq(d[i])
	skip = next(l)
	for x in l
		I = I + cos(sqrt(sum(norminv(x).^2/2)))
		N = N + 1
		if abs(pi^(d[i]/2)/N*I-exact[i])/abs(exact[i]) < 1e-3
			break
		end
	end
	println("  Converged after $(N) points!")
end

println("============================")
println("=    COMPUTE INTEGRAL 2    =")
println("============================")
# same, but now with randomized sequence
d = [9, 25, 60, 80, 100]
q = 16
for i in 1:length(d)
	println("Computing cos-integral in $(d[i]) dimensions...")
	N = 0
	I = zeros(q)
	l = LatSeq(d[i])
	r = RandWrapper(l,q)
	skip = next(r)
	for x in r
		for j in 1:q
			@inbounds I[j] = I[j] + cos(sqrt(sum(norminv(x[:,j]).^2/2)))
		end
		N = N + 1
		if std(pi^(d[i]/2)/N*I)/abs(mean(pi^(d[i]/2)/N*I)) < 5e-4
			break
		end
	end
	println("  Converged after $(N) points!")
end

println("============================")
println("= COMPUTE INTEGRAL 2 (DIG) =")
println("============================")
# same, but now with randomized sequence
d = [9, 25, 60, 80, 100]
q = 16
for i in 1:length(d)
	println("Computing cos-integral in $(d[i]) dimensions...")
	N = 0
	I = zeros(q)
	l = DigSeq(d[i])
	r = RandWrapper(l,q)
	skip = next(r)
	for x in r
		for j in 1:q
			@inbounds I[j] = I[j] + cos(sqrt(sum(norminv(x[:,j]).^2/2)))
		end
		N = N + 1
		if std(pi^(d[i]/2)/N*I)/abs(mean(pi^(d[i]/2)/N*I)) < 5e-4
			break
		end
	end
	println("  Converged after $(N) points!")
end
