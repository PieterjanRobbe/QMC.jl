println(LOAD_PATH)
println(ENV["TRAVIS_JULIA_VERSION"])

# obtain absolute path to module
PATH = LOAD_PATH[find([length(search(LOAD_PATH[i],"QMC")) for i in 1:length(LOAD_PATH)].!=0)][1]

# the 250-dimensional generating vector of the Cools, Kuo and Nuyens paper in SIAM SISC for 2^20 points
const CKN_250_20 = collect(readdlm(PATH*"/data/CKN_250_20.txt",UInt32))

# 3600-dimensional generating vector with product weights 1/j^2 from 
# http://web.maths.unsw.edu.au/~fkuo/lattice/
const K_3600_20 = collect(readdlm(PATH*"/data/K_3600_20.txt",UInt32))
