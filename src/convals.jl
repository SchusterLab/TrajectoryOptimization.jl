"""
convals.jl
"""

mutable struct ConstraintParams{T}
	ϕ::T  	    # penalty scaling parameter
	μ0::T 	    # initial penalty parameter
	μ_max::T    # max penalty parameter
	λ_max::T    # max Lagrange multiplier
end

function ConstraintParams(ϕ::T1=10., μ0::T2=1., μ_max::T3=1e8, λ_max::T4=1e8) where {T1,T2,T3,T4}
	T = promote_type(T1, T2, T3, T4)
	return ConstraintParams{T}(T(ϕ), T(μ0), T(μ_max), T(λ_max))
end

"""
	ConVal{C,V,M,W}

Holds information about a constraint
"""
struct ConVal{C,Tc,Tic,T}
    con::C
    # constraint function value
    c::Tc
    # dual
    λ::Tc
    # penalty multiplier
    μ::Tc
    # active constraints
    a::Tc
    # index for this conval's `c` in
    # the global linearized constraint
    c_ginds::Tic
    params::ConstraintParams{T}
end

# constructors
function ConVal(con::C, n::Int, m::Int, c_ginds::Tic,
                M, V, penalty_scaling::T, penalty_initial::T, penalty_max::T, dual_max::T
                ) where {C,T,Tic,Tix,Tiu}
    p = length(con)
    c = V(zeros(p))
    λ = V(zeros(p))
    μ = V(fill(penalty_initial, p))
    a = V(ones(Bool, p))
    params = ConstraintParams(penalty_scaling, penalty_initial,
                              penalty_max, dual_max)
    Tc = typeof(c)
    return ConVal{C,Tc,Tic,T}(con, c, λ, μ, a, c_ginds, params)
end

# methods
function update_active!(a::AbstractVector, con::AbstractConstraint,
                        c::AbstractVector, λ::AbstractVector, tol::Real)
    if con.sense == equality
        a.= true
    elseif con.sense == inequality
        for i = 1:length(a)
            a[i] = ((c[i] >= tol) | (abs(λ[i]) > tol))
        end
    end
    return nothing
end

function violation(con::AbstractConstraint, c::AbstractVector)
    viol = 0.
    if con.sense == equality
        viol = norm(c, Inf)
    elseif con.sense == inequality
        viol = max(0., maximum(c))
    end
    return viol
end

function update_dual_penalty!(convals::Vector{Vector{ConVal}})
    for convals_ in convals
        for conval in convals_
            # update dual
            λ_max = conval.params.λ_max
            λ_min = conval.con.sense == equality ? -λ_max : zero(λ_max)
            for i in eachindex(conval.λ)
                conval.λ[i] = clamp(conval.λ[i] + conval.μ[i] * conval.c[i], λ_min, λ_max)
            end
            # update penalty
            for i in eachindex(conval.μ)
                conval.μ[i] = clamp(conval.params.ϕ * conval.μ[i], 0, conval.params.μ_max)
            end
        end
    end
end

function max_violation_penalty(convals::Vector{Vector{ConVal}})
    max_violation = 0.
    max_penalty = 0.
    for (k, convals_) in enumerate(convals)
        for conval in convals_
            viol = violation(conval.con, conval.c)
            max_violation = max(max_violation, viol)
            max_penalty = max(max_penalty, maximum(conval.μ))
        end
    end
    return max_violation, max_penalty
end


# build list of convals from constraint list
function convals_from_constraint_list(cons::ConstraintList)
    # indices for tracking global concatenated constraint
    c_gind = 0
    convals = Vector{ConVal}[]
    # build a list of convals at each knot point
    for k = 1:cons.N
        convals_ = ConVal[]
        # iterate over each constraint
        for i in 1:length(cons)
            con = cons.constraints[i]
            knot_points = cons.inds[i]
            # if this constraint is active at the current knot point,
            # add a conval for it to the current list
            if k in knot_points
                p = length(con)
                c_ginds = cons.V((1:p) .+ c_gind)
                c_gind += p
                conval = ConVal(con, cons.n, cons.m, c_ginds,
                                cons.M, cons.V, 0., 0., 0., 0.)
                push!(convals_, conval)
            end
        end
        push!(convals, convals_)
    end
    return convals
end
