using JuMP, HiGHS

#=
Original code by Thomas Moore
@denbigh
included in https://github.com/ypaul21/Clapeyron.jl/pull/56
=#
"""
    HELDTPFlash(;numphases = 2;max_steps = 1e4*(numphases-1),population_size =20,time_limit = Inf,verbose = false, logspace = false)

Method to solve non-reactive multicomponent flash problem by finding global minimum of Gibbs Free Energy via Differential Evolution.

User must assume a number of phases, `numphases`. If true number of phases is smaller than numphases, model should predict either (a) identical composition in two or more phases, or (b) one phase with negligible total number of moles. If true number of phases is larger than numphases, a thermodynamically unstable solution will be predicted.

The optimizer will stop at `max_steps` evaluations or at `time_limit` seconds

"""
Base.@kwdef struct HELDTPFlash <: TPFlashMethod
    max_steps::Int = 1e4
    time_limit::Float64 = Inf
    eps_λ::Float64 = 0.5
    eps_b::Float64 = 1e-2
    eps_x::Float64 = 1e-3
    eps_η::Float64 = 1e-3
    eps_μ::Float64 = 1e-6
    eps_g::Float64 = 1e-6
    p_tun::Float64 = 1e-3
    verbose::Bool = false
    logspace::Bool = false
end

index_reduction(flash::HELDTPFlash,z) = flash

function tp_flash_impl(model::EoSModel, p, T, n, method::HELDTPFlash)
    nc = length(model)
    p̃ = method.p_tun
    TYPE = typeof(p+T+first(n))
    # Stage I: Stability test and initialisation
    if method.verbose == true
        println("==========================================")
        println("Stage I: Stability test and initialisation")
        println("==========================================")
        println("----------------------------")
        println("Step 1: Stability test at n₀")
        println("----------------------------")
    end
    vₛ = log10(volume(model,p,T,n))
    μ₀ = VT_chemical_potential(model,10 .^vₛ,T,n)
    xₛ = prepend!(deepcopy(n),vₛ)

    d(x) = (eos(model,10 .^x[1],T,x[2:end])+p*10 .^x[1])/sum(x[2:end])-∑(x[2:end].*μ₀)
    if method.verbose == true
        println("Initial point found. Beginning tunneling")
    end
    i = 0
    while i<10nc
        f(x) = d(x)*exp(p̃/√(∑((xₛ[2:end]-x[2:end]).^2)))
        x0 = rand(length(n))
        x0 = x0./sum(x0)
        v0 = log10(x0_volume(model,p,T,x0))
        x0 = prepend!(x0,v0)
        r = Solvers.optimize(f,x0)
        xₙ = Solvers.x_sol(r)
        dₙ = d(xₙ)
        if dₙ<0
            if method.verbose == true
                println("Negative tangent found. Moving on to step 2.")
            end
            break
        end
        i+=1
    end
    if i==10nc
        if method.verbose == true
            println("No negative tangent found. Initial point is stable.")
            println("Terminating HELD")
        end
        return (xˢ,gibbs_free_energy(model,p,T,xˢ))
    end
    
    if method.verbose == true
        println("--------------------------------------")
        println("Step 2: Initialisation of dual problem")
        println("--------------------------------------")
    end
    k = 0
    UBDⱽ = gibbs_free_energy(model,p,T,n)/R̄/T
    (ℳ,G,λᴸ,λᵁ) = initial_candidate_phases(model,p,T,n)
    
    if method.verbose == true
        println("Iteration counter set to k="*string(k))
        println("Upper bound set to UBDⱽ="*string(UBDⱽ))
        println("ℳ initialised")
        println("===================================================")
        println("Stage II: Identification of candidate stable phases")
        println("===================================================")
        println("-------------------------------------------------------")
        println("Step 3: Solve the outer problem (OPₓᵥ) at iteration k="*string(k))
        println("-------------------------------------------------------")
    end

    OPₓᵥ = Model(HiGHS.Optimizer)
    @variable(OPₓᵥ, v)
    @variable(OPₓᵥ, λ[1:nc-1])
    @constraint(OPₓᵥ,v<=UBDⱽ)
    @constraint(OPₓᵥ,[i ∈ 1:length(G)],v<=G[i]+∑(λ.*(n[1:nc-1] .-ℳ[i,1:nc-1])))
    @constraint(OPₓᵥ,[i ∈ 1:nc-1],λᴸ[i]<=λ[i]<=λᵁ[i])
    @objective(OPₓᵥ, Max, v)
    println(OPₓᵥ)
    optimize!(OPₓᵥ)
    λˢ = JuMP.value.(λ)
    if method.verbose == true
        println("-------------------------------------------------------")
        println("Step 4: Solve the inner problem (IPₓᵥ) at iteration k="*string(k))
        println("-------------------------------------------------------")
    end
    Lⱽᵏ = deepcopy(UBDⱽ)+1
    i = 0
    while i < method.max_steps
        println(i)
        Lⱽ(x) = (eos(model,10 .^x[1],T,Fractions.FractionVector(x[2:end]))+p*10 .^x[1])/R̄/T+sum(λˢ[j]*(n[j]-x[j+1]) for j ∈ 1:nc-1)
        x0 = rand(length(n))
        x0 = x0./sum(x0)
        v0 = log10(x0_volume(model,p,T,x0))
        x0 = prepend!(x0,v0)
        r = Solvers.optimize(Lⱽ,x0[1:end-1])
        xᵏ = Solvers.x_sol(r)
        Lⱽᵏ = Lⱽ(xᵏ)
        if Lⱽᵏ<UBDⱽ
            ℳ = [ℳ;vcat(xᵏ[2:end],1-sum(xᵏ[2:end]),xᵏ[1])']
            println(ℳ)
            G = append!(G,Lⱽᵏ)
            break
        end
        i+=1
    end
    println(G)
end

function initial_candidate_phases(model,p,T,n)
    nc = length(n)
    x̂ = zeros(nc-1,nc+1)
    x̄ = zeros(nc-1,nc+1)
    Ĝ = zeros(nc-1)
    Ḡ = zeros(nc-1)
    λᵁ = zeros(nc)
    λᴸ = zeros(nc)
    μ̂ = zeros(nc-1,nc)
    μ̄ = zeros(nc-1,nc)
    for i ∈ 1:nc-1
        x̂[i,i] = n[i]/2
        x̄[i,i] = (1+n[i])/2
        for k ∈ 1:nc-1
            if k != i
                x̂[i,k] = (1-x̂[i,i])/(nc-1)
                x̄[i,k] = (1-x̄[i,i])/(nc-1)
            end
        end
        x̂[i,nc]=1-sum(x̂[i,:])
        x̂[i,end] = log10(volume(model,p,T,x̂[i,1:nc]))
        Ĝ[i] = VT_gibbs_free_energy(model,10 .^x̂[i,end],T,x̂[i,1:nc])/R̄/T
        μ̂[i,:] = VT_chemical_potential(model,10 .^x̂[i,end],T,x̂[i,1:nc])/R̄/T

        x̄[i,nc]=1-sum(x̄[i,:])
        x̄[i,end]=log10(volume(model,p,T,x̄[i,1:nc]))
        Ḡ[i] = VT_gibbs_free_energy(model,10 .^x̄[i,end],T,x̄[i,1:nc])/R̄/T
        μ̄[i,:] = VT_chemical_potential(model,10 .^x̄[i,end],T,x̄[i,1:nc])/R̄/T
    end
    μ = [μ̂;μ̄]
    μ = μ[:,1:nc-1]-μ[:,nc]
    λᵁ = maximum(μ[:,1:nc-1];dims=1)
    λᴸ = minimum(μ[:,1:nc-1];dims=1)
    return ([x̂;x̄],[Ĝ;Ḡ],λᴸ,λᵁ)
end


"""
    Obj_de_tp_flash(model,p,T,z,dividers,numphases,logspace = false)

Function to calculate Gibbs Free Energy for given partition of moles between phases.
This is a little tricky. 

We must find a way of uniquely converting a vector of numbers,
each in (0, 1), to a partition. We must be careful that 
the mapping is 1-to-1, not many-to-1, as if many inputs
map to the same physical state in a redundant way, there
will be multiple global optima, and the global optimization 
will perform poorly.

Our approach is to specify (numphases-1) numbers in (0,1) for
each species. We then scale these numbers systematically in order to partition
the species between the phases. Each set of (numphases - 1) numbers
will result in a unique partition of the species into the numphases
phases.
"""

export HELDTPFlash