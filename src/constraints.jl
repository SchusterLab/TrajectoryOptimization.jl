"""
constraints.jl

Notes:
The jacobian! method is meant for Markovian constraints
while the jacobian_copy! method is intended for a direct solver
which accomodates more general constraints
"""

import RobotDynamics: state_dim, control_dim

const NULL_MAT = Array{Float64,2}(undef, 0, 0)
const NULL_VEC = Array{Float64,1}(undef, 0)

# general

@enum ConstraintSense begin
    EQUALITY = 1
    INEQUALITY = 2
end

abstract type AbstractConstraint end

"Only a function of states and controls at a single knotpoint"
abstract type StageConstraint <: AbstractConstraint end
"Only a function of states at a single knotpoint"
abstract type StateConstraint <: StageConstraint end
"Only a function of controls at a single knotpoint"
abstract type ControlConstraint <: StageConstraint end
"Only a function of states and controls at two adjacent knotpoints"
abstract type CoupledConstraint <: AbstractConstraint end
"Only a function of states at adjacent knotpoints"
abstract type CoupledStateConstraint <: CoupledConstraint end
"Only a function of controls at adjacent knotpoints"
abstract type CoupledControlConstraint <: CoupledConstraint end

"""
	GoalConstraint{P,T}

Constraint of the form ``x_g = a``, where ``x_g`` can be only part of the state
vector.

# Constructors:
```julia
GoalConstraint(xf::AbstractVector)
GoalConstraint(xf::AbstractVector, inds)
```
where `xf` is an n-dimensional goal state. If `inds` is provided,
only `xf[inds]` will be used.
"""
struct GoalConstraint{Tx,Ti,Txp,Tup,Tp,Tpx,Tpu} <: StateConstraint
    n::Int
    m::Int
    # constraint length
    p::Int
    # goal state vector
    xf::Tx
    # indices into state vector
    inds::Ti
    # tmp for cost_derivatives!
    XP_tmp::Txp
    UP_tmp::Tup
    p_tmp::Vector{Tp}
    # stores for jacobian!
    Cx::Tpx
    Cu::Tpu
    # misc
    const_jac::Bool
    state_expansion::Bool
    control_expansion::Bool
    coupled_expansion::Bool
    direct::Bool
    sense::ConstraintSense
end

# constructors
function GoalConstraint(n::Int, m::Int, xf::Tx, inds::Ti, M, V; direct::Bool=false) where {Tx,Ti}
    xf_ = xf[inds]
    p = length(inds)
    XP_tmp = M(zeros(n, p))
    UP_tmp = M(zeros(n, p))
    p_tmp = [V(zeros(p)) for i = 1:2]
    Cx = M(zeros(p, n))
    Cu = M(zeros(p, m))
    Txp = typeof(XP_tmp)
    Tup = typeof(UP_tmp)
    Tp = typeof(p_tmp[1])
    Tpx = typeof(Cx)
    Tpu = typeof(Cu)
    const_jac = true
    state_expansion = true
    control_expansion = false
    coupled_expansion = false
    sense = EQUALITY
    con = GoalConstraint{Tx,Ti,Txp,Tup,Tp,Tpx,Tpu}(
        n, m, p, xf_, inds, XP_tmp, UP_tmp, p_tmp, Cx, Cu, const_jac,
        state_expansion, control_expansion, coupled_expansion, direct, sense
    )
    jacobian!(con.Cx, con.Cu, con, NULL_VEC, NULL_VEC, 0)
    return con
end

# evaluation
function evaluate!(c::AbstractVector, con::GoalConstraint,
                   X::AbstractVector, U::AbstractVector, k::Int)
    for (i, j) in enumerate(con.inds)
        c[i] = X[k][j] - con.xf[j] 
    end
    return nothing
end

function jacobian!(Cx::AbstractMatrix, Cu::AbstractMatrix, con::GoalConstraint, X::AbstractVector,
                   U::AbstractVector, k::Int)
    # ASSUMPTION: Cu .= 0
    T = eltype(Cx)
    for (i, j) in enumerate(con.inds)
	    Cx[i, j] = one(T)
    end
    return nothing
end

function jacobian_copy!(D::AbstractMatrix, con::GoalConstraint,
                        X::AbstractVector, U::AbstractVector, k::Int,
                        c_ginds::AbstractVector, x_ginds::AbstractVector,
                        u_ginds::AbstractVector)
    # ASSUMPTION: D[c_ginds, i] .= 0 for i ∉ x_ginds[k]
    for (i, j) in enumerate(con.inds)
        D[c_ginds[i], x_ginds[k][j]] = 1.
    end
    return nothing
end

# methods
@inline Base.length(con::GoalConstraint) = con.p

function max_violation_info(con::GoalConstraint, c::AbstractVector, k::Int)
    max_viol = -Inf
    info_str = ""
    for i in con.inds
        viol = abs(c[i])
        if viol > max_viol
            info_str = "GoalConstraint x[$i] k=$k"
            max_viol = viol
        end
    end
    return max_viol, info_str
end



"""
DynamicsConstraint - constraint for explicit dynamics
"""
struct DynamicsConstraint{T,Tir,Tm,Tix,Tiu,Tx,Txx,Txu,Txz} <: CoupledConstraint
    n::Int
    m::Int
    ir::Tir
    model::Tm
    ts::Vector{T}
    ix::Tix
    iu::Tiu
    # store for evaluate!
    x_tmp::Tx
    # store for jacobian_copy!
    A::Txx
    B::Txu
    AB::Txz
    # misc
    const_jac::Bool
    direct::Bool
    sense::ConstraintSense
end

# constructors
function DynamicsConstraint(n::Int, m::Int, ir::Tir, model::Tm, ts::Vector{T},
                            ix::Tix, iu::Tiu, M, V; direct::Bool=true) where {Tir,Tm,T,Tix,Tiu}
    x_tmp = V(zeros(n))
    A = M(zeros(n, n))
    B = M(zeros(n, m))
    AB = M(zeros(n, n+m))
    Tx = typeof(x_tmp)
    Txx = typeof(A)
    Txu = typeof(B)
    Txz = typeof(AB)
    const_jac = false
    sense = EQUALITY
    con = DynamicsConstraint{T,Tir,Tm,Tix,Tiu,Tx,Txx,Txu,Txz}(
        n, m, ir, model, ts, ix, iu, x_tmp, A, B, AB, const_jac, direct, sense
    )
    return con
end

# evaluation
function evaluate!(c::AbstractVector, con::DynamicsConstraint, X::AbstractVector,
                   U::AbstractVector, k::Int)
    RD.discrete_dynamics!(con.x_tmp, con.ir, con.model, X[k - 1], U[k - 1], con.ts[k - 1],
                          con.ts[k] - con.ts[k - 1])
    for i = 1:con.n
        c[i] = con.x_tmp[i] - X[k][i]
    end
    return nothing
end

function jacobian!(Cx::AbstractMatrix, Cu::AbstractMatrix, con::DynamicsConstraint,
                   X::AbstractVector, U::AbstractVector, k::Int)
    throw("not implemented")
    return nothing
end

function jacobian_copy!(D::AbstractMatrix, con::DynamicsConstraint,
                        X::AbstractVector, U::AbstractVector, k::Int,
                        c_ginds::AbstractVector, x_ginds::AbstractVector,
                        u_ginds::AbstractVector)
    # ASSUMPTION: D[c_ginds, i] .= 0 for i ∉ x_ginds[k - 1] U u_ginds[k - 1] U x_ginds[k]
    RD.discrete_jacobian!(con.AB, con.A, con.B, con.ir, con.model,
                          X[k - 1], U[k - 1], con.ts[k - 1],
                          con.ts[k] - con.ts[k - 1], con.ix, con.iu)
    D[c_ginds, x_ginds[k - 1]] .= con.A
    D[c_ginds, u_ginds[k - 1]] .= con.B
    for (i, j) in enumerate(c_ginds)
        D[j, x_ginds[k][i]] = 1.
    end
end

# methods
@inline Base.length(con::DynamicsConstraint) = con.n

function max_violation_info(con::DynamicsConstraint, c::AbstractVector, k::Int)
    max_viol = -Inf
    info_str = ""
    for i = 1:con.n
        viol = abs(c[i])
        if viol > max_viol
            info_str = "DynamicsConstraint x[$i] k=$k"
            max_viol = viol
        end
    end
    return max_viol, info_str
end


"""
	BoundConstraint{Tz,Tiu,Til,Tinds}

Linear bound constraint on states and controls
# Constructors
```julia
BoundConstraint(n, m; x_min, x_max, u_min, u_max)
```
Any of the bounds can be ±∞. The bound can also be specifed as a single scalar, which applies the bound to all state/controls.
"""
struct BoundConstraint{Tx,Tu,Tixu,Tixl,Tiuu,Tiul,Ti,Txp,Tup,Tp,Tpx,Tpu} <: StageConstraint
    n::Int
    m::Int
    # constraint length
    p::Int
    x_max::Tx
    x_min::Tx
    u_max::Tu
    u_min::Tu
    x_max_inds::Tixu
    x_min_inds::Tixl
    u_max_inds::Tiuu
    u_min_inds::Tiul
    inds::Ti
    XP_tmp::Txp
    UP_tmp::Tup
    p_tmp::Vector{Tp}
    Cx::Tpx
    Cu::Tpu
    const_jac::Bool
    state_expansion::Bool
    control_expansion::Bool
    coupled_expansion::Bool
    direct::Bool
    sense::ConstraintSense
end

# constructor
function BoundConstraint(n::Int, m::Int, x_max::Tx, x_min::Tx,
                         u_max::Tu, u_min::Tu, M, V; direct::Bool=false, checks=true) where {
                             Tx<:AbstractVector,Tu<:AbstractVector}
    if checks
        @assert all(x_max .>= x_min)
        @assert all(u_max .>= u_min)
    end
    # TODO: the construction for these indices can be done in a less
    # disgusting way by finding all finite among x_max, u_max, etc. individually
    # get constraint indices
    b = V([-x_max; -u_max; x_min; u_min])
    inds = V(findall(isfinite, b))
    # indices for evaluation (constraint ind, vector ind)
    state_expansion = control_expansion = false
    x_max_inds = Tuple{Int,Int}[]
    x_min_inds = Tuple{Int,Int}[]
    u_max_inds = Tuple{Int,Int}[]
    u_min_inds = Tuple{Int,Int}[]
    for (c_ind, z_ind) in enumerate(inds)
        if z_ind <= n
            state_expansion = true
            x_ind = z_ind
            insert!(x_max_inds, length(x_max_inds) + 1, (c_ind, x_ind))
        elseif n < z_ind <= n + m
            control_expansion = true
            u_ind = z_ind - n
            insert!(u_max_inds, length(u_max_inds) + 1, (c_ind, u_ind))
        elseif n + m < z_ind <= n + m + n
            state_expansion = true
            x_ind = z_ind - n - m
            insert!(x_min_inds, length(x_min_inds) + 1, (c_ind, x_ind))
        else # nm + n < i
            control_expansion = true
            u_ind = z_ind - n - m - n
            insert!(u_min_inds, length(u_min_inds) + 1, (c_ind, u_ind))
        end
    end
    x_max_inds = V(x_max_inds)
    x_min_inds = V(x_min_inds)
    u_max_inds = V(u_max_inds)
    u_min_inds = V(u_min_inds)
    coupled_expansion = state_expansion && control_expansion
    # tmps for jacobian!
    p = length(inds)
    XP_tmp = M(zeros(n, p))
    UP_tmp = M(zeros(m, p))
    p_tmp = [V(zeros(p)), V(zeros(p))]
    Cx = M(zeros(p, n))
    Cu = M(zeros(p, m))
    const_jac = true
    sense = INEQUALITY
    # types
    Tixu = typeof(x_max_inds)
    Tixl = typeof(x_min_inds)
    Tiuu = typeof(u_max_inds)
    Tiul = typeof(u_min_inds)
    Ti = typeof(inds)
    Txp = typeof(XP_tmp)
    Tup = typeof(UP_tmp)
    Tp = typeof(p_tmp[1])
    Tpx = typeof(Cx)
    Tpu = typeof(Cu)
    # construct
    con = BoundConstraint{Tx,Tu,Tixu,Tixl,Tiuu,Tiul,Ti,Txp,Tup,Tp,Tpx,Tpu}(
        n, m, p, x_max, x_min, u_max, u_min, x_max_inds, x_min_inds, u_max_inds,
        u_min_inds, inds, XP_tmp, UP_tmp, p_tmp, Cx, Cu, const_jac, state_expansion,
        control_expansion, coupled_expansion, direct, sense
    )
    # initialize
    jacobian!(Cx, Cu, con, x_max, x_max, 0)
    return con
end

# evaluation
function evaluate!(c::AbstractVector, con::BoundConstraint, X::AbstractVector,
                   U::AbstractVector, k::Int; log=false)
    for (i, j) in con.x_max_inds
        c[i] = X[k][j] - con.x_max[j]
        # if k == 79 && log
        #     println("c[$i]: $(c[i]), X[$k][$j]: $(X[k][j]), x_max[j]: $(con.x_max[j])")
        # end
    end
    for (i, j) in con.u_max_inds
        c[i] = U[k][j] - con.u_max[j]
    end
    for (i, j) in con.x_min_inds
        c[i] = con.x_min[j] - X[k][j]
        # if k == 79 && log
        #     println("c[$i]: $(c[i]), X[$k][$j]: $(X[k][j]), x_min[j]: $(con.x_min[j])")
        # end
    end
    for (i, j) in con.u_min_inds
        c[i] = con.u_min[j] - U[k][j]
    end
    return nothing
end

function jacobian!(Cx::AbstractMatrix, Cu::AbstractMatrix, con::BoundConstraint,
                   X::AbstractVector, U::AbstractVector, k::Int)
    for (i, j) in con.x_max_inds
        Cx[i, j] = 1.
    end
    for (i, j) in con.u_max_inds
        Cu[i, j] = 1.
    end
    for (i, j) in con.x_min_inds
        Cx[i, j] = -1.
    end
    for (i, j) in con.u_min_inds
        Cu[i, j] = -1.
    end
    return nothing
end

function jacobian_copy!(D::AbstractMatrix, con::BoundConstraint,
                        X::AbstractVector, U::AbstractVector, k::Int,
                        c_ginds::AbstractVector, x_ginds::AbstractVector,
                        u_ginds::AbstractVector)
    for (i, j) in con.x_max_inds
        D[c_ginds[i], x_ginds[k][j]] = 1.
    end
    for (i, j) in con.u_max_inds
        D[c_ginds[i], u_ginds[k][j]] = 1.
    end
    for (i, j) in con.x_min_inds
        D[c_ginds[i], x_ginds[k][j]] = -1.
    end
    for (i, j) in con.u_min_inds
        D[c_ginds[i], u_ginds[k][j]] = -1.
    end
    return nothing
end

# methods
@inline Base.length(con::BoundConstraint) = con.p

function max_violation_info(con::BoundConstraint, c::AbstractVector, k::Int)
    max_viol = -Inf
    info_str = ""
    for (i, j) in con.x_max_inds
        viol = c[i]
        if viol > max_viol
            max_viol = viol
            info_str = "BoundConstraint x_max[$j] k=$k"
        end
    end
    for (i, j) in con.u_max_inds
        viol = c[i]
        if viol > max_viol
            max_viol = viol
            info_str = "BoundConstraint u_max[$j] k=$k"
        end
    end
    for (i, j) in con.x_min_inds
        viol = c[i]
        if viol > max_viol
            max_viol = viol
            info_str = "BoundConstraint x_min[$j] k=$k"
        end
    end
    for (i, j) in con.u_min_inds
        viol = c[i]
        if viol > max_viol
            max_viol = viol
            info_str = "BoundConstraint u_min[$j] k=$k"
        end
    end
    return max_viol, info_str
end
