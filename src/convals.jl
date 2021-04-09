"""
convals.jl
"""

mutable struct ConstraintParams{T}
	ϕ::T  	    # penalty scaling parameter
	μ0::T 	    # initial penalty parameter
	μ_max::T    # max penalty parameter
	λ_max::T    # max Lagrange multiplier
end

function ConstraintParams(ϕ::T1 = 10, μ0::T2 = 1.0, μ_max::T3 = 1e8, λ_max::T4 = 1e8) where {T1,T2,T3,T4}
	T = promote_type(T1,T2,T3,T4)
	ConstraintParams(T(ϕ), T(μ0), T(μ_max), T(λ_max))
end

"""
	ConVal{C,V,M,W}

Holds information about a constraint
"""
struct ConVal{C,Tc,Tcx,Tcu,T}
    con::C
    c::Tc # constraint function value
    Cx::Tcx # cf derivative w.r.t. x
    Cu::Tcu # cf derivative w.r.t. u
    λ::Tc # dual
    μ::Tc # penalty multiplier
    a::Tc # active constraints
    params::ConstraintParams{T}
end

# constructors
function ConVal(con::C, n::Int, m::Int, V, M, penalty_scaling::T,
                penalty_initial::T, penalty_max::T, dual_max::T) where {C, T}
    p = length(con)
    c = V(zeros(p))
    sense_ = sense(con)
    if sense_ isa StateConstraint
        Cu = nothing
        Cx = M(zeros(p, n))
    elseif sense_ isa ControlConstraint
        Cu = M(zeros(p, m))
        Cx = nothing
    else
        Cu = M(zeros(p, m))
        Cx = M(zeros(p, n))
    end
    λ = V(zeros(p))
    μ = V(fill(penalty_initial, p))
    a = V(ones(Bool, p))
    params = ConstraintParams(penalty_scaling, penalty_initial,
                              penalty_max, dual_max)
    Tc = typeof(c)
    Tcx = typeof(Cx)
    Tcu = typeof(Cu)
    return ConVal{C,Tc,Tcx,Tcu,T}(con, c, Cx, Cu, λ, μ, a, params)
end

# methods
@inline violation(::Equality, v) = norm(v, Inf)
@inline violation(::Inequality, v) = maximum(v)
@inline violation(::Equality, v::Real) = abs(v)
@inline violation(::Inequality, v::Real) = v > 0 ? v : 0.0

function violation(conval::ConVal)
    s = sense(cval.con)
    return violation(s, conval.c)
end

function dual_update!(conval::ConVal)
    # TODO in-place ops
    λ_max = conval.params.λ_max
    λ_min = sense(conval.con) == Equality() ? -λ_max : zero(λ_max)
    conval.λ .= clamp.(conval.λ + conval.μ .* conval.c, λ_min, λ_max)
    return nothing
end

function penalty_update!(conval::ConVal)
    # TODO in place opts
    conval.μ .= clamp.(conval.params.ϕ .* conval.μ, 0., conval.params.μ_max)
end

# @inline get_data(A::AbstractArray) = A
# @inline get_data(A::SizedArray) = A.data
# 
# # ***
# function ConVal(n::Int, m::Int, con::C, inds::Ti, 
# 		jac::Tj, vals::Tv, M, V, iserr::Bool=false) where {C,Ti,Tj,Tv}
#     if !iserr && size(gen_jacobian(con)) != size(jac[1])
# 	throw(DimensionMismatch("size of jac[i] $(size(jac[1])) does not match "
#                                 "the expected size of $(size(gen_jacobian(con)))"))
#     end
#     vals2 = deepcopy(vals)
#     p = length(con)
#     P = length(vals)
#     views = [gen_views(∇c, con, n, m) for ∇c in jac]
#     ∇x = [v[1] for v in views]
#     ∇u = [v[2] for v in views]
#     c_max = V(zeros(P))
#     is_const = [V(zeros(Bool,P)), V(zeros(Bool,P))]
#     Tx = typeof(∇x)
#     Tu = typeof(∇u)
#     Tc = typeof(c_max)
#     Tb = typeof(is_const[1])
#     return ConVal{C,Ti,Tv,Tj,Tx,Tu,Tc,Tb}(con, inds, vals, vals2, jac, ∇x, ∇u, c_max,
#                                           is_const, iserr)
# end

# function ConVal(n::Int, m::Int, cval::ConVal, M, V)
#     # create a ConVal for the "raw" Jacobians, if needed
#     # 	otherwise return the same ConVal
#     if cval.iserr
# 	p = length(cval.con)
# 	ws = widths(cval.con, n, m)
# 	jac = [M(zeros(p,w)) for k in cval.inds, w in ws]
# 	cval_ = ConVal(n, m, cval.con, cval.inds, jac, cval.vals, M, V, false)
#     else
#         cval_ = cval
#     end
#     return cval_
# end

# function ConVal(n::Int, m::Int, con::AbstractConstraint, inds::AbstractVector, iserr::Bool=false)
#     C, c = gen_convals(n,m,con,inds)
#     ConVal(n, m, con, inds, C, c)
# end

# function _index(cval::ConVal, k::Int)
# 	if k ∈ cval.inds
# 		return k - cval.inds[1] + 1
# 	else
# 		return 0
# 	end
# end

# function evaluate!(cval::ConVal, Z::AbstractTrajectory)
# 	evaluate!(cval.vals, cval.con, Z, cval.inds)
# end

# # ***
# function jacobian!(cval::ConVal, Z::AbstractTrajectory, init::Bool=false)
#     if cval.iserr
# 	throw(ErrorException("Can't evaluate Jacobians directly on the error state Jacobians"))
#     else
# 	jacobian!(cval.jac, cval.con, Z, cval.inds)
#     end
# end

# function ∇jacobian!(G, cval::ConVal, Z::AbstractTrajectory, λ, init::Bool=false)
# 	∇jacobian!(G, cval.con, Z, λ, cval.inds, cval.is_const[2], init)
# end

@inline norm_violation(::Equality, v, p=2) = norm(v,p)

@inline function norm_violation(::Inequality, v, p=2)
	# TODO: try this with LazyArrays?
	if p == 1
		a = zero(eltype(v))
		for x in v
			a += max(x,0)
		end
		return a
	elseif p == 2
		a = zero(eltype(v))
		for x in v
			a += max(x, 0)^2
		end
		return sqrt(a)
	elseif p == Inf
		return maximum(v)
	else
		throw(ArgumentError("$p is not a valid norm value. Must be 1,2 or Inf"))
	end
end

function norm_violation(cval::ConVal, p=2)
	norm_violation!(cval, p)
	return norm(cval.c_max, p)
end

function norm_violation!(cval::ConVal, p=2)
	s = sense(cval.con)
	for i in eachindex(cval.inds)
		cval.c_max[i] = norm_violation(s, cval.vals[i], p)
	end
end

function norm_dgrad!(cval::ConVal, Z::AbstractTrajectory, p=1)
	for (i,k) in enumerate(cval.inds)
		zs = RobotDynamics.get_z(cval.con, Z, k)
		mul!(cval.vals2[i], cval.jac[i,1], zs[1])
		if length(zs) > 1
			mul!(cval.vals2[i], cval.jac[i,2], zs[2], 1.0, 1.0)
		end
		cval.c_max[i] = norm_dgrad(cval.vals[i], cval.vals2[i], p)
	end
	return nothing
end
"""
	    dgrad(x, dx, p=1)
Directional derivative of `norm(x, p)` in the direction `dx`
"""
function norm_dgrad(x, dx, p=1)
    g = zero(eltype(x))
    if p == 1
        @assert length(x) == length(dx)
        g = zero(eltype(x))
        for i in eachindex(x)
            if x[i] < 0
                g += -dx[i]
            elseif x[i] > 0
                g += dx[i]
            else
                g += abs(dx[i])
            end
        end
    else
        throw("Directional derivative of $p-norm isn't implemented yet")
    end
    return g
end

function norm_residual!(res, cval::ConVal, λ::Vector{<:AbstractVector}, p=2)
    for (i,k) in enumerate(cval.inds)
	mul!(res[i], cval.jac[i,1], λ[i])
	if size(cval.jac,2) > 1
	    mul!(res[i], cval.jac[i,2], λ[i], 1.0, 1.0)
	end
	cval.c_max[i] = norm(res[i], p)
    end
    return nothing
end

function error_expansion!(errval::ConVal, conval::ConVal, model::AbstractModel, G)
    if errval.jac !== conval.jac
	for (i,k) in enumerate(conval.inds)
	    mul!(errval.∇x[i], conval.∇x[i], get_data(G[k]))
	    errval.∇u[i] .= conval.∇u[i]
	end
    end
end

function error_expansion!(con::AbstractConstraint, err, jac, G)
    mul!(err, jac, G)
end

# ***
function gen_convals(n̄::Int, m::Int, con::AbstractConstraint, inds::AbstractVector, V, M)
    # n̄ is the state diff size
    p = length(con)
    ws = widths(con, n̄,m)
    C = [M(zeros(p,w)) for k in inds, w in ws]
    c = [V(zeros(p)) for k in inds]
    return C, c
end

function gen_convals(D::AbstractMatrix, d::AbstractVector, cinds, zinds,
                     con::AbstractConstraint, inds)
    P = length(inds)
    p = length(con)
	n,m = get_dims(con, length(zinds[1]))
    ws = widths(con, n, m)

    C = [begin
		view(D, cinds[i], zinds[k+(j-1)][1:ws[j]])
	end for (i,k) in enumerate(inds), j = 1:length(ws)]
    c = [view(d, cinds[i]) for i in 1:P]
    return C,c
end

function gen_convals(blocks::Vector, cinds, con::AbstractConstraint, inds)
    # assumes all cinds are contiguous indices (i.e. can be represented as a UnitRange)
    C1 = map(enumerate(inds)) do (i,k)
        nm = size(blocks[k].Y,2)
		if con isa StateConstraint
			iz = 1:width(con)
		elseif con isa ControlConstraint
			m = control_dim(con)
			n = nm - m
			iz = n .+ (1:m)
		else
			iz = 1:nm
		end
		ic = cinds[i][1]:cinds[i][end]
		n1 = size(blocks[k].D2, 1)
        view(blocks[k].Y, n1 .+ (ic), iz)
    end
	C2 = map(enumerate(inds)) do (i,k)
		if con isa StageConstraint
			w = size(blocks[k].Y,2)
			view(blocks[k].Y,1:0,1:w)
		else
			w = size(blocks[k+1].Y,2)
			n = state_dim(con)
			view(blocks[k+1].Y,1:n,1:w)
		end
	end
	C = [C1 C2]
    c = map(enumerate(inds)) do (i,k)
		ic = cinds[i][1]:cinds[i][end]
        view(blocks[k].y, ic)
    end
    return C,c
end
