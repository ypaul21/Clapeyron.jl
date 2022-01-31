abstract type HVRuleModel <: MixingRule end

struct HVRule{γ} <: HVRuleModel
    components::Array{String,1}
    activity::γ
    references::Array{String,1}
end
@registermodel HVRule

"""
    HVRule{γ} <: HVRuleModel
    
    HVRule(components::Vector{String};
    activity = Wilson,
    userlocations::Vector{String}=String[],
    activity_userlocations::Vector{String}=String[],
    verbose::Bool=false)

## Input Parameters

None

## Model Parameters

None

## Input models 

- `activity`: Activity Model

## Description

Huron-Vidal Mixing Rule 
```
aᵢⱼ = √(aᵢaⱼ)(1-kᵢⱼ)
bᵢⱼ = (bᵢ +bⱼ)/2
b̄ = ∑bᵢⱼxᵢxⱼ
c̄ = ∑cᵢxᵢ
ā = b̄(∑[xᵢaᵢᵢαᵢ/(bᵢᵢ)] - gᴱ/λ)
if the model is Peng-Robinson:
    λ = 0.6232252401402305
if the model is Redlich-Kwong:
    λ = 0.6931471805599453
```
## References
1. Huron, M.-J., & Vidal, J. (1979). New mixing rules in simple equations of state for representing vapour-liquid equilibria of strongly non-ideal mixtures. Fluid Phase Equilibria, 3(4), 255–271. doi:10.1016/0378-3812(79)80001-1

"""
HVRule

export HVRule
function HVRule(components::Vector{String}; activity = Wilson, userlocations::Vector{String}=String[],activity_userlocations::Vector{String}=String[], verbose::Bool=false)
    init_activity = activity(components;userlocations = activity_userlocations,verbose)
    
    references = ["10.1016/0378-3812(79)80001-1"]
    model = HVRule(components, init_activity,references)
    return model
end

HVλ(::PRModel) =  1/(2*√(2))*log((2+√(2))/(2-√(2)))
HVλ(::RKModel) = log(2)

function mixing_rule(model::CubicModel,V,T,z,mixing_model::HVRuleModel,α,a,b,c)
    n = sum(z)
    invn2 = (one(n)/n)^2
    b̄ = dot(z,Symmetric(b),z) * invn2
    c̄ = dot(z,c)/n
    gᴱ = excess_gibbs_free_energy(mixing_model.activity,1e5,T,z)*invn
    ∑ab = sum(z[i]*a[i,i]*α[i]/b[i,i] for i ∈ @comps)*invn
    _λ = HVλ(mixing_model,model)
    ā = b̄*(∑ab-gᴱ/_λ)
    return ā,b̄,c̄
end

