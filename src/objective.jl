"""
objective.jl
"""

############################################################################################
#                              OBJECTIVES                                                  #
############################################################################################

abstract type AbstractObjective end
state_dim(obj::AbstractObjective) = throw(ErrorException("state_dim not implemented"))
control_dim(obj::AbstractObjective) = throw(ErrorException("control_dim not implemented"))

"""$(TYPEDEF)
Objective: stores stage cost(s) and terminal cost functions

Constructors:
```julia
Objective(cost, cost_term, N)
```
"""
struct Objective{Tsc,Ttc} <: AbstractObjective
    stage_cost::Tsc
    terminal_cost::Ttc
    N::Int
end

# constructors
function Objective(stage_cost::Tsc, terminal_cost::Ttc, N::Int, checks=true) where {
    Tsc <: CostFunction, Ttc <: CostFunction}
    if checks
        @assert stage_cost.terminal != true
        @assert terminal_cost.terminal == true
    end
    return Objective{Tsc,Ttc}(stage_cost, terminal_cost, N)
end

# methods
Base.copy(obj::Objective) = Objective(copy(obj.stage_cost), copy(obj.terminal_cost), obj.N)
Base.show(io::IO, obj::Objective) = print(io,"Objective")

@inline control_dim(obj::Objective) = control_dim(obj.stage_cost[1])
@inline state_dim(obj::Objective) = state_dim(obj.stage_cost[1])


# ############################################################################################
# #                            Quadratic Objectives (Expansions)
# ############################################################################################
# const QuadraticObjective{n,m,T} = Objective{QuadraticCost{n,m,T,Matrix{T},Matrix{T}}}
# const QuadraticExpansion{n,m,T} = Objective{<:QuadraticCostFunction{n,m,T}}
# const DiagonalCostFunction{n,m,T} = Union{DiagonalCost{n,m,T},QuadraticCost{n,m,T,<:Diagonal,<:Diagonal}}

# function QuadraticObjective(n::Int, m::Int, N::Int, isequal::Bool=false)
#     Objective([QuadraticCost{Float64}(n,m, terminal=(k==N) && !isequal) for k = 1:N])
# end

# function QuadraticObjective(obj::QuadraticObjective, model::AbstractModel)
#     # Create QuadraticObjective linked to error cost expansion
#     @assert RobotDynamics.state_diff_size(model) == size(model)[1]
#     return obj
# end

# function QuadraticObjective(obj::QuadraticObjective, model::LieGroupModel)
#     # Create QuadraticObjective partially linked to error cost expansion
#     @assert length(obj[1].q) == RobotDynamics.state_diff_size(model)
#     n,m = size(model)
#     costfuns = map(obj.cost) do costfun
#         Q = zeros(n,n)
#         R = costfun.R
#         H = zeros(m,n)
#         q = zeros(n)
#         r = costfun.r
#         c = costfun.c
#         QuadraticCost(Q,R,H,q,r,c, checks=false, terminal=costfun.terminal)
#     end
#     Objective(costfuns)
# end


# # Convenience constructors
# @doc raw"""
#     LQRObjective(Q, R, Qf, xf, N)

# Create an objective of the form
# `` (x_N - x_f)^T Q_f (x_N - x_f) + \sum_{k=0}^{N-1} (x_k-x_f)^T Q (x_k-x_f) + u_k^T R u_k``

# Where `eltype(obj) <: DiagonalCost` if `Q`, `R`, and `Qf` are
#     `Union{Diagonal{<:Any,<:StaticVector}}, <:StaticVector}`
# """
# function LQRObjective(Q::AbstractArray, R::AbstractArray, Qf::AbstractArray,
#         xf::AbstractVector, N::Int; checks=true, uf=@SVector zeros(size(R,1)))
#     @assert size(Q,1) == length(xf)
#     @assert size(Qf,1) == length(xf)
#     @assert size(R,1) == length(uf)
#     n = size(Q,1)
#     m = size(R,1)
#     H = SizedMatrix{m,n}(zeros(m,n))
#     q = -Q*xf
#     r = -R*uf
#     c = 0.5*xf'*Q*xf + 0.5*uf'R*uf
#     qf = -Qf*xf
#     cf = 0.5*xf'*Qf*xf

#     ℓ = QuadraticCost(Q, R, H, q, r, c, checks=checks)
#     ℓN = QuadraticCost(Qf, R, H, q, r, cf, checks=false, terminal=true)

#     Objective(ℓ, ℓN, N)
# end

# function LQRObjective(
#         Q::Union{<:Diagonal, <:AbstractVector},
#         R::Union{<:Diagonal, <:AbstractVector},
#         Qf::Union{<:Diagonal, <:AbstractVector},
#         xf::AbstractVector, N::Int;
#         uf=(@SVector zeros(size(R,1))),
#         checks=true)
#     n,m = size(Q,1), size(R,1)
#     @assert size(Q,1) == length(xf)
#     @assert size(Qf,1) == length(xf)
#     @assert size(R,1) == length(uf)
#     Q,R,Qf = Diagonal(Q), Diagonal(R), Diagonal(Qf)
#     q = -Q*xf
#     r = -R*uf
#     c = 0.5*xf'*Q*xf + 0.5*uf'R*uf
#     qf = -Qf*xf
#     cf = 0.5*xf'*Qf*xf

#     ℓ = DiagonalCost(Q, R, q, r, c, checks=checks, terminal=false)

#     ℓN = DiagonalCost(Qf, R, qf, r, cf, checks=false, terminal=true)

#     Objective(ℓ, ℓN, N)
# end
