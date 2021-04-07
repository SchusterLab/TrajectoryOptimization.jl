# export
#     cost

############################################################################################
#                              COST METHODS                                                #
############################################################################################

"""
	stage_cost(cost::CostFunction, z::AbstractKnotPoint)

Evaluate the cost at a knot point, and automatically handle terminal
knot point, multiplying by dt as necessary."""
function stage_cost(cost::CostFunction, z::AbstractKnotPoint)
    if is_terminal(z)
        stage_cost(cost, state(z))
    else
        stage_cost(cost, state(z), control(z))*z.dt
    end
end

"""
	cost(obj::Objective, Z::Traj)
	cost(obj::Objective, dyn_con::DynamicsConstraint{Q}, Z::Traj)

Evaluate the cost for a trajectory. If a dynamics constraint is given,
    use the appropriate integration rule, if defined.
"""
function cost(obj::AbstractObjective, Z::AbstractTrajectory)
    cost!(obj, Z)
    J = get_J(obj)
    return sum(J)
end

# ForwardDiff-able method
function cost(obj::AbstractObjective, Z::AbstractTrajectory)
    J = 0.
    for k = 1:length(obj)
        J += stage_cost(obj[k], Z[k])
    end
    return J
end

# Default to no integration
cost(obj, dyn_con::DynamicsConstraint{Q}, Z) where Q<:QuadratureRule = cost(obj, Z)

"Evaluate the cost for a trajectory (non-allocating)"
@inline function cost!(obj::Objective, Z::AbstractTrajectory)
    map!(stage_cost, obj.J, obj.cost, Z)
end


############################################################################################
#                              COST EXPANSIONS                                             #
############################################################################################

"""
	cost_gradient!(E::Objective, obj::Objective, Z, init)

Evaluate the cost gradient along the entire tracjectory `Z`, storing the result in `E`.

If `init == true`, all gradients will be evaluated, even if they are constant.
"""
function cost_gradient!(E::Objective, obj::Objective, Z::AbstractTrajectory, init::Bool=false)
    is_const = E.const_grad
    N = length(Z)
    for k in eachindex(Z)
        if init || !is_const[k]
			if is_terminal(Z[k])
            	is_const[k] = gradient!(E.cost[k], obj.cost[k], state(Z[k]))
			else
            	is_const[k] = gradient!(E.cost[k], obj.cost[k], state(Z[k]), control(Z[k]))
			end
            dt_x = k < N ? Z[k].dt :  one(Z[k].dt)
            dt_u = k < N ? Z[k].dt : zero(Z[k].dt)
            E[k].q .*= dt_x
            E[k].r .*= dt_u
        end
    end
end

"""
	cost_hessian!(E::Objective, obj::Objective, Z, init)

Evaluate the cost hessian along the entire tracjectory `Z`, storing the result in `E`.

If `init == true`, all hessian will be evaluated, even if they are constant. If false,
they will only be evaluated if they are not constant.
"""
function cost_hessian!(E::Objective, obj::Objective, Z::AbstractTrajectory, init::Bool=false, rezero::Bool=false)
    is_const = E.const_hess
    N = length(Z)
    for k in eachindex(Z)
        if init || !is_const[k]
			if rezero
				E[k].Q .= 0
				E[k].R .= 0
				!is_blockdiag(E[k]) && (E[k].H .= 0)
			end
			if is_terminal(Z[k])
            	is_const[k] = hessian!(E.cost[k], obj.cost[k], state(Z[k]))
			else
            	is_const[k] = hessian!(E.cost[k], obj.cost[k], state(Z[k]), control(Z[k]))
			end
            dt_x = k < N ? Z[k].dt :  one(Z[k].dt)
            dt_u = k < N ? Z[k].dt : zero(Z[k].dt)
            E[k].Q .*= dt_x
            E[k].R .*= dt_u
            E[k].H .*= dt_u
        end
    end
end

"""
	cost_expansion!(E::Objective, obj::Objective, Z, [init, rezero])

Evaluate the 2nd order Taylor expansion of the objective `obj` along the trajectory `Z`,
storing the result in `E`.

If `init == false`, the expansions will only be evaluated if they are not constant.

If `rezero == true`, all expansions will be multiplied by zero before taking the expansion.
"""
function cost_expansion!(E::Objective, obj::Objective, Z::Traj, init::Bool=false, rezero::Bool=false)
    cost_gradient!(E, obj, Z, init)
    cost_hessian!(E, obj, Z, init, rezero)
    return nothing
end

function error_expansion!(E::Objective, Jexp::Objective, model::AbstractModel, Z::Traj, G, tmp=G[end])
    @assert E === Jexp "E and Jexp should be the same object for AbstractModel"
    return nothing
end

function error_expansion!(E::Objective, Jexp::Objective, model::LieGroupModel, Z::Traj, G, tmp=G[end])
    for k in eachindex(E.cost)
        error_expansion!(E.cost[k], Jexp.cost[k], model, Z[k], G[k], tmp)
    end
end

function error_expansion!(E::QuadraticCost, cost::QuadraticCost, model, z::AbstractKnotPoint,
        G, tmp)
	E.Q .= 0
    RobotDynamics.∇²differential!(E.Q, model, state(z), cost.q)
    if size(model)[1] < 15
        G = SMatrix(G)
        E.H .= SMatrix(cost.H) * G
        E.q .= G'cost.q
        E.Q .+= G'cost.Q*G
    else
        mul!(E.H, cost.H, G)
        mul!(E.q, Transpose(G), cost.q)
        mul!(tmp, cost.Q, G)
        mul!(E.Q, Transpose(G), tmp, 1.0, 1.0)
    end
    return nothing
end

struct StaticExpansion{T,N,M,NN,MM,NM} <: AbstractExpansion{T}
	x::SVector{N,T}
	xx::SMatrix{N,N,T,NN}
	u::SVector{M,T}
	uu::SMatrix{M,M,T,MM}
	ux::SMatrix{M,N,T,NM}
end

function StaticExpansion(E::AbstractExpansion)
	StaticExpansion(SVector(E.x), SMatrix(E.xx),
		SVector(E.u), SMatrix(E.uu), SMatrix(E.ux))
end

function StaticExpansion(x,xx,u,uu,ux)
	StaticExpansion(SVector(x), SMatrix(xx), SVector(u), SMatrix(uu), SMatrix(ux))
end

function static_expansion(cost::QuadraticCost)
	StaticExpansion(cost.q, cost.Q, cost.r, cost.R, cost.H)
end

"""
	dgrad(E::QuadraticExpansion, dZ::Traj)

Calculate the derivative of the cost in the direction of `dZ`, where `E` is the current
quadratic expansion of the cost.
"""
function dgrad(E::QuadraticExpansion{n,m,T}, dZ::Traj)::T where {n,m,T}
	g = zero(T)
	N = length(E)
	for k = 1:N-1
		g += dot(E[k].q, state(dZ[k])) + dot(E[k].r, control(dZ[k]))
	end
	g += dot(E[N].q, state(dZ[N]))
	return g
end

"""
	dhess(E::QuadraticCost, dZ::Traj)

Calculate the scalar 0.5*dZ'G*dZ where G is the hessian of cost
"""
function dhess(E::QuadraticExpansion{n,m,T}, dZ::Traj)::T where {n,m,T}
	h = zero(T)
	N = length(E)
	for k = 1:N-1
		x = state(dZ[k])
		u = control(dZ[k])
		h += dot(x, E[k].Q, x) + dot(u, E[k].R, u)
	end
	x = state(dZ[N])
	h += dot(x, E[N].Q, x)
	return 0.5*h
end

"""
	norm_grad(E::QuadraticExpansion, p=2)

Norm of the cost gradient
"""
function norm_grad(E::QuadraticExpansion{n,m,T}, p=2)::T where {n,m,T}
	J = get_J(E)
	for (k,cost) in enumerate(E)
		J[k] = norm(cost.q, p)
		if !cost.terminal
			J[k] = norm(SA[norm(cost.r, p) J[k]], p)
		end
	end
	return norm(J, p)
end
