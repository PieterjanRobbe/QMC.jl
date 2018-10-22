module QMC

# dependencies
using DelimitedFiles

# import statements
import Base: show, next, reset, iterate

# export statements
export LatSeq, DigSeq, RandWrapper, getPoint, next, reset, iterate # from QMCgenerators.jl

# include statements
include("QMCdata.jl")

include("QMCgenerators.jl")

end # module
