"""
    slle_solubility(model::CompositeModel, p, T)

Calculates the phase boundary for solid-liquid-liquid equilibriumm of a **ternary** mixture, at a given temperature and pressure.
Returns a matrix containing the composition of the two liquids phases.

Can only function when solid and liquid models are specified within a CompositeModel and when the third component is the solute.
"""
function slle_solubility(model::CompositeModel,p,T)
    p = p*one(eltype(model))
    T = T*one(eltype(model))
    if length(model.components) != 3
        error("SLLE can only be obtained for ternary systems")
    end
    solute = model.components[end]
    idx_sol = solute.==model.components

    solid_r,idx_sol_r = index_reduction(model.solid,idx_sol)
    μsol = chemical_potential(solid_r,p,T,[1.])
    μ_ref = reference_chemical_potential(model.fluid,p,T,reference_chemical_potential_type(model.fluid))
    x0 = x0_slle_solubility(model,p,T,μsol)
    f!(F,x) = obj_slle_solubility(F,model.fluid,p,T,[exp10(x[1]),1-exp10(x[1])-exp10(x[2]),exp10(x[2])],[exp10(x[3]),1-exp10(x[3])-exp10(x[4]),exp10(x[4])],μsol,μ_ref)
    results = Solvers.nlsolve(f!,x0,LineSearch(Newton()))
    sol = exp10.(Solvers.x_sol(results))
    x1 = [sol[1],1-sol[1]-sol[2],sol[2]]
    x2 = [sol[3],1-sol[3]-sol[4],sol[4]]
    return (x1,x2)
end

function obj_slle_solubility(F,liquid,p,T,x1,x2,μsol,μ_ref)
    nc = length(liquid)
    γliq1 = activity_coefficient(liquid,p,T,x1,μ_ref = μ_ref)
    μliq1 = @. Rgas()*T*log(γliq1*x1)
    γliq2 = activity_coefficient(liquid,p,T,x2,μ_ref = μ_ref)
    μliq2 = @. Rgas()*T*log(γliq2*x2)

    F[1] = (μliq1[end] - μsol[1])/Rgas()/T
    F[2:nc+1] = @. (μliq2 - μliq1)/Rgas()/T
    return F
end

function x0_slle_solubility(model,p,T,μsol)
    x0 = zeros(length(model),2)
    for i in 1:length(model)-1
        z∞ = ones(length(model))*1e-3
        z∞[i] = 1.0
        z∞ ./= sum(z∞)
        γ∞ = activity_coefficient(model.fluid,p,T,z∞)
        z∞[end] = exp(μsol[1]/(Rgas()*T)-log(γ∞[end]))
        z∞ ./= sum(z∞)
        γ∞ = activity_coefficient(model.fluid,p,T,z∞)
        x = zeros(length(model))
        x[i] = 1.0
        x[end] = exp(μsol[1]/(Rgas()*T)-log(γ∞[end]))
        x[findfirst(x.==0)] = 1/γ∞[findfirst(x.==0)]
        x ./= sum(x)

        # x = x0_lle_init(model.fluid,p,T,z0,x)

        x0[:,i] = x
    end

    x0 = [x0[1,1],x0[3,1],x0[1,2],x0[3,2]]

    return log10.(x0)
end
