struct UNIQUACParam <: EoSParam
    a::PairParam{Float64}
    r::SingleParam{Float64}
    q::SingleParam{Float64}
    q_p::SingleParam{Float64}
    Mw::SingleParam{Float64}
end

abstract type UNIQUACModel <: ActivityModel end

struct UNIQUAC{c<:EoSModel} <: UNIQUACModel
    components::Array{String,1}
    icomponents::UnitRange{Int}
    params::UNIQUACParam
    puremodel::Vector{c}
    absolutetolerance::Float64
    references::Array{String,1}
end
@registermodel UNIQUAC
export UNIQUAC

function UNIQUAC(components::Vector{String}; puremodel=PR,
    userlocations=String[], 
     verbose=false)
    params = getparams(components, ["Activity/UNIQUAC/UNIQUAC_like.csv", "properties/molarmass.csv","Activity/UNIQUAC/UNIQUAC_unlike.csv"]; userlocations=userlocations, asymmetricparams=["a"], ignore_missing_singleparams=["a"], verbose=verbose)
    a  = params["a"]
    r  = params["r"]
    q  = params["q"]
    q_p = params["q_p"]
    Mw  = params["Mw"]
    icomponents = 1:length(components)
    
    init_puremodel = [puremodel([components[i]]) for i in icomponents]
    packagedparams = UNIQUACParam(a,r,q,q_p,Mw)
    references = String[]
    model = UNIQUAC(components,icomponents,packagedparams,init_puremodel,1e-12,references)
    return model
end
#=
function lnγ_comb(model::UNIQUACModel,p,T,z)
    #Φ =  x.*r/sum(x[i]*r[i] for i ∈ @comps)
    #θ = x.*q/sum(x[i]*q[i] for i ∈ @comps)
    q = model.params.q.values
    r = model.params.r.values
    ∑z = sum(z)
    Φ_mean = dot(z,q)/∑z
    θ_mean = dot(z,r)/∑z
    _0 = zero(p+T+first(z))
    lnγ = fill(_0,length(model))
    for i in @comps
        zi = z[i]
        Φᵢ = zi*q[i]/Φ_mean
        θᵢ = zi*r[i]/θ_mean
        lnγ[i] =  log(Φᵢ/zi)+(1-Φᵢ/zi)-5*q[i]*(log(Φᵢ/θᵢ)+(1-Φᵢ/θᵢ))
    end
    return lnγ
end
=#
function activity_coefficient(model::UNIQUACModel,p,T,z)
    a = model.params.a.values
    q = model.params.q.values
    q_p = model.params.q_p.values
    r = model.params.r.values

    x = z ./ sum(z)

    Φ =  x.*r/sum(x[i]*r[i] for i ∈ @comps)
    θ = x.*q/sum(x[i]*q[i] for i ∈ @comps)
    θ_p = x.*q_p/sum(x[i]*q_p[i] for i ∈ @comps)
    τ = @f(Ψ)
    lnγ_comb = @. log(Φ/x)+(1-Φ/x)-5*q*(log(Φ/θ)+(1-Φ/θ))
    lnγ_res  = q_p.*(1 .-log.(sum(θ_p[i]*τ[i,:] for i ∈ @comps)) .-sum(θ_p[i]*τ[:,i]/sum(θ_p[j]*τ[j,i] for j ∈ @comps) for i ∈ @comps))
    return exp.(lnγ_comb+lnγ_res)
end

function Ψ(model::UNIFACModel,V,T,z)
    Tinv = 1/T
    a = model.params.a.values
    return @. exp(-a*Tinv)
end