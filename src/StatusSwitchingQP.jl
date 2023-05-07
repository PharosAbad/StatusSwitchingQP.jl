module StatusSwitchingQP
using LinearAlgebra

export Status, IN, DN, UP, OE, EO, Event
export LP, QP


include("./types.jl")

include("./Simplex.jl")
using .Simplex
export SimplexLP

#include("./SSLP.jl")
#using .SSLP


#include("./SSQP.jl")
#using .SSQP
#export solveQP, QP


end
