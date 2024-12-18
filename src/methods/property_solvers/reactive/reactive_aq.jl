function equilibrate(model::ReactiveAqModel,p,T,n0;z0=nothing)
    ν = stoichiometric_coefficient(model) # Exists
    Keq = ideal_Keq(model,T,n0,ν)
    println(Keq)
    if isnothing(z0)
        z0 = x0_equilibrium_conditions(model,p,T,n0)
    end
    f!(F,z) = equilibrium_conditions(model,F,p,T,n0,exp10.(z),ν,Keq)
    sol = Solvers.nlsolve(f!,log10.(z0),TrustRegion(Newton(), NWI()),NEqOptions(),ForwardDiff.Chunk{length(model.reactions)}()) # Exists
    println(sol)
    ξ = exp10.(Solvers.x_sol(sol))
    println(ξ)
    n = n0+ν*ξ
    return n
end

function pH(model::ReactiveAqModel,p,T,n0;z0=nothing)
    n = equilibrate(model, p, T, n0; z0 = z0)
    a = aqueous_activity(model.eosmodel,p,T,n)
    hydronium_id = find_hydronium_index(model)
    hydroxide_id = find_hydroxide_index(model)

    idx_w = find_water_indx(model)

    zref = ones(length(model.components)).*1e-30
    zref[idx_w] = 1.
    zref[hydronium_id] = 0.01801528*1.
    zref[hydroxide_id] = 0.01801528*1.

    zref ./= sum(zref)

    μref = chemical_potential(model.eosmodel,p,T,zref)
    μmix = chemical_potential(model.eosmodel,p,T,n)

    a = exp.((μmix .- μref)./R̄./T)

    return -log10(a[hydronium_id])
end

function equilibrium_conditions(model::ReactiveAqModel,F,p,T,n0,ξ,ν,Keq)
    n = n0+ν*ξ

    nreact = length(model.reactions)
    nspecies = length(model.components)

    μmix = chemical_potential(model.eosmodel,p,T,n)
    idx_w = find_water_indx(model)
    for i in 1:nreact
        μref = zeros(nspecies)

        # 1. Charged Reactants
        ireact_chrg = findall(x->ν[x,i]<0 && model.eosmodel.charge[x] != 0,1:nspecies)

        zref = ones(nspecies).*1e-30
        zref[idx_w] = 1.
        zref[ireact_chrg] .= 0.01801528*abs.(ν[ireact_chrg,i])
        zref ./= sum(zref)

        μref[ireact_chrg] = chemical_potential(model.eosmodel,p,T,zref)[ireact_chrg]

        # 2. Charged Products
        iprod_chrg = findall(x->ν[x,i]>0 && model.eosmodel.charge[x] != 0,1:nspecies)
        zref = ones(nspecies).*1e-30
        zref[idx_w] = 1.
        zref[iprod_chrg] .= 0.01801528*abs.(ν[iprod_chrg,i])
        zref ./= sum(zref)

        μref[iprod_chrg] = chemical_potential(model.eosmodel,p,T,zref)[iprod_chrg]

        # 3. Neutral species
        ineutral = findall(x->model.eosmodel.charge[x] == 0,1:nspecies)
        for j in ineutral
            zref = ones(nspecies).*1e-30
            zref[j] = 0.01801528*abs(ν[j,i])
            zref[idx_w] = 1.
            zref ./= sum(zref)

            μref[j] = chemical_potential(model.eosmodel,p,T,zref)[j]
        end

        a = exp.((μmix .- μref)./R̄./T)

        F[i] = sum(log.(a).*ν[:,i]) - log(Keq[i])
    end
    return F
end