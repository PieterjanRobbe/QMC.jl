module QMC

# import statements
import Base: show, next, start, done

# export statements
export LatSeq, RandLatSeq, next!, reset!, getPoint, ndims, nshifts, start, next, done # from QMCgenerators.jl

# include statements
include("QMCdata.jl")

include("QMCgenerators.jl")

end # module
