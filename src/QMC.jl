module QMC

# import statements
import Base: show, next, reset, start, done

# export statements
export LatSeq, DigSeq, RandWrapper, next, reset, start, done # from QMCgenerators.jl

# include statements
include("QMCdata.jl")

include("QMCgenerators.jl")

end # module
