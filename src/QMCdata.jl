# obtain absolute path to module
PATH = LOAD_PATH[find([length(search(LOAD_PATH[i],"QMC")) for i in 1:length(LOAD_PATH)].!=0)][1]

# the 250-dimensional generating vector of the Cools, Kuo and Nuyens paper in SIAM SISC for 2^20 points
<<<<<<< HEAD
const CKN_250_20 = collect(readdlm(PATH*"/data/CKN_250_20.txt",UInt32))

# 3600-dimensional generating vector with product weights 1/j^2 from 
# http://web.maths.unsw.edu.au/~fkuo/lattice/
const K_3600_32 = collect(readdlm(PATH*"/data/K_3600_32.txt",UInt32))

# 32 dimensional generating matrix from the Niederreiter Xing directory on
# https://people.cs.kuleuven.be/~dirk.nuyens/qmc-generators/nxmats/
const N_32_32 = reshape(collect(readdlm(PATH*"/data/N_32_32.txt",UInt32)),(32,32))

# Sobol matrix of Joe and Kuo paper in ACM Transactions on Mathematical Software
const JK_1111_32 = reshape(collect(readdlm(PATH*"/data/JK_1111_32.txt",UInt32)),(1111,32))
=======
const CKN_250_20 = vec(readdlm(PATH*"/data/CKN_250_20.txt",UInt32))

# 3600-dimensional generating vector with product weights 1/j^2 from 
# http://web.maths.unsw.edu.au/~fkuo/lattice/
const K_3600_20 = vec(readdlm(PATH*"/data/K_3600_20.txt",UInt32))
>>>>>>> 9b9d4a66b77977997fbe65eb2c0f79c2807135b0
