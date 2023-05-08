"Status Switching Method for Quadratic Programming and Linear Programming"
module StatusSwitchingQP
using LinearAlgebra


include("./types.jl")
export Status, IN, DN, UP, OE, EO, Event
export LP, QP, Settings

include("./utils.jl")
export getRowsGJ, getRowsGJr



include("./Simplex.jl")
using .Simplex
export SimplexLP

include("./SSLP.jl")
using .SSLP
export solveLP


include("./SSQP.jl")
using .SSQP
export solveQP


end
