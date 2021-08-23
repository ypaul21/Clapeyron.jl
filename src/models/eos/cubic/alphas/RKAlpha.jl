abstract type RKAlphaModel <: AlphaModel end

struct RKAlphaParam <: EoSParam
end


@newmodelsimple RKAlpha RKAlphaModel RKAlphaParam
is_splittable(::RKAlpha) = false

export RKAlpha
function RKAlpha(components::Vector{String}; userlocations::Vector{String}=String[], verbose::Bool=false)
    params = getparams(components, ["properties/critical.csv"]; userlocations=userlocations, verbose=verbose)
    acentricfactor = SingleParam(params["w"],"acentric factor")
    packagedparams = RKAlphaParam(acentricfactor)
    model = RKAlpha(packagedparams, verbose=verbose)
    return model
end

function α_function(model::CubicModel,V,T,z,alpha_model::RKAlphaModel)
    Tc = model.params.Tc.values
    α = @. 1 /√(T/Tc)
    return α
end