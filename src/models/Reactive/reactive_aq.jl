struct ReactiveAqParams{T} <: ParametricEoSParam{T}
    ΔGf::SingleParam{T}
    ΔHf::SingleParam{T}
end

struct ReactiveAqModel{T} <: ReactiveEoSModel where T <: EoSModel
    components::Vector{String}
    reactions::Vector
    params::ReactiveAqParams{Float64}
    eosmodel::T
end

reference_chemical_potential_type(model::ReactiveAqModel) = :aqueous

function ReactiveAqModel(_components, reactions::Vector;
                       model::EoSModel,
                       userlocations=String[],
                       verbose=false)
    components = format_components(_components)
    params = getparams(components, ["properties/formation_aqueous.csv"]; userlocations = userlocations, verbose = verbose)
    ΔGf = params["Gf"]
    ΔHf = params["Hf"]
    packagedparams = ReactiveAqParams{Float64}(ΔGf,ΔHf)
    return ReactiveAqModel(components,reactions,packagedparams,model)
end

function ideal_Keq(model::ReactiveAqModel,T,z,ν)
    T0 = 298.15
    ΔGf0 = model.params.ΔGf.values
    ΔHf = model.params.ΔHf.values
    ΔGf = ΔGf0/T0-ΔHf*(1/T-1/T0)
    ΔrG = sum(ΔGf.*ν,dims=1)
    Keq = exp.(-ΔrG/Rgas(model.eosmodel))
end

export ReactiveAqModel
