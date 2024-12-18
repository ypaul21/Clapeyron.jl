abstract type ReactiveEoSModel <: EoSModel end

struct ReactiveParams{T} <: ParametricEoSParam{T}
    ΔHf::SingleParam{T}
    Sf::SingleParam{T}
end

struct ReactiveModel{T} <: ReactiveEoSModel
    components::Vector{String}
    reactions::Vector
    params::ReactiveParams{Float64}
    eosmodel::T
end

function ReactiveModel(_components, reactions::Vector;
                       model::EoSModel,
                       userlocations=String[],
                       verbose=false)
    components = format_components(_components)
    params = getparams(components, ["properties/formation.csv"]; userlocations = userlocations, verbose = verbose)
    ΔHf = params["Hf"]
    Sf = params["Sf"]
    packagedparams = ReactiveParams{Float64}(ΔHf,Sf)
    return ReactiveModel(components,reactions,packagedparams,model)
end

reference_chemical_potential_type(model::ReactiveEoSModel) = :pure

ideal_Keq(model::ReactiveEoSModel,T,z) = ideal_Keq(model,T,z,stoichiometric_coefficient(model))

function ideal_Keq(model::ReactiveModel,T,z,ν)
    ΔHf = model.params.ΔHf.values
    Sf = model.params.Sf.values
    ΔGf = ΔHf-T*Sf
    ΔrG = sum(ΔGf.*ν,dims=1)
    Keq = exp.(-ΔrG/Rgas(model.eosmodel)/T)
end

function stoichiometric_coefficient(model::ReactiveEoSModel)
    ν = zeros(length(model.components),length(model.reactions))
    for (i,reaction) in enumerate(model.reactions)
        species = [reaction[k][1] for k in 1:length(reaction)]
        coeff = [reaction[k][2] for k in 1:length(reaction)]
        for (j,component) in enumerate(model.components)
            if component in species
                ν[j,i] = coeff[findfirst(isequal(component),species)]
            end
        end
    end
    return ν
end

function Base.show(io::IO, mime::MIME"text/plain", model::ReactiveModel)
    print(io, typeof(model))
    length(model.eosmodel) == 1 && println(io, " with 1 component:")
    length(model.eosmodel) > 1 && println(io, " with ", length(model.eosmodel), " components:")
    show_pairs(io,model.components)
    println(io, "")

    length(model.reactions) == 1 && println(io, "and 1 reaction:")
    length(model.reactions) > 1 && println(io, "and ", length(model.reactions), " reactions:")
    for reaction in model.reactions
        species = [reaction[k][1] for k in 1:length(reaction)]
        coeff = [reaction[k][2] for k in 1:length(reaction)]
        reactants = ""
        products = ""
        for i in 1:length(species)
            if coeff[i] < 0
                reactants *= string(abs(coeff[i]))*" "*species[i]*" + "
            else
                products *= string(abs(coeff[i]))*" "*species[i]*" + "
            end
        end
        reactants = reactants[1:end-3]
        products = products[1:end-3]
        println(io, reactants, " <=> ", products)
    end
end

export ReactiveModel
