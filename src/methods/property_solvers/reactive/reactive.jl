function equilibrate(model::ReactiveEoSModel,p,T,n0;z0=nothing)
    ν = stoichiometric_coefficient(model) # Exists
    Keq = ideal_Keq(model,T,n0,ν)
    if isnothing(z0)
        z0 = x0_equilibrium_conditions(model,p,T,n0)
    end
    μ_ref_type = reference_chemical_potential_type(model)
    μ_ref = reference_chemical_potential(model.eosmodel,p,T,μ_ref_type)
    f!(F,z) = equilibrium_conditions(model,F,p,T,n0,exp10.(z),ν,Keq,μ_ref)
    sol = Solvers.nlsolve(f!,log10.(z0),TrustRegion(Newton(), NWI()),NEqOptions(),ForwardDiff.Chunk{1}()) # Exists
    ξ = exp10.(Solvers.x_sol(sol))
    n = n0+ν*ξ
    return n
end

function equilibrium_conditions(model::ReactiveModel,F,p,T,n0,ξ,ν,Keq,μ_ref)
    n = n0+ν*ξ
    a = activity_impl(model.eosmodel,p,T,n,μ_ref,:unknown,:unknown,true,nothing)
    F[1:end] = sum(log.(a).*ν,dims=1) - log.(Keq)
    return F
end

function x0_equilibrium_conditions(model::ReactiveEoSModel,p,T,n0)
    z0 = zeros(length(model.components))
    z0[1] = log10(p)
    return z0
end

export equilibrate
