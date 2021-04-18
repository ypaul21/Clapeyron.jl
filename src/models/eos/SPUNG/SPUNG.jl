#=
SPUNG = State Research Program for Utilization of Natural Gas

based on:
Unification of the two-parameter equation of state and the principle ofcorresponding states
Jørgen Mollerup

the idea is to use mollerup shape factors and PropaneRef as the EoS
but any function that provides a_scaled and shape_factors will suffice
=#

struct SPUNG{S<:EoSModel,E<:EoSModel} <: EoSModel
    shape_model::S
    shape_ref::S
    model_ref::E
end

function SPUNG(components::Vector{String},refmodel=PropaneRef(),shapemodel::SHAPE=SRK(components)) where SHAPE<:EoSModel
    refname = component_names(refmodel)
    @show 1
    shape_ref = SHAPE.name.wrapper(refname)
    return SPUNG(shapemodel,shape_ref,refmodel)
end


function eos(model::SPUNG,V,T,z=SA[1.0],phase="unknown")
    f,h = shape_factors(model,V,T,z)
    T0 = T/f
    V0 = V/h
    return eos(model.model_ref,V0,T0)
end

function shape_factors(model::SPUNG{<:ABCubicModel},V,T,z=SA[1.0])
    n = sum(z)
    x = z * (1/n)
    a,b = cubic_ab(model.shape_model,T,x)
    a0,b0 = cubic_ab(model.shape_ref,T,SA[1.0])
    h = b/b0
    fh = n*a/a0
    f = fh/h
    return f,h
end


function Base.show(io::IO,mime::MIME"text/plain",model::SPUNG)
    println(io,"Extended Corresponding States model")
    println(io," reference model: ",model.model_ref)
    print(io," shape model: ",model.shape_model)
end

function Base.show(io::IO,model::SPUNG)
    print(io,"SPUNG(",model.shape_ref,",",model.model_ref,")")
end

function lb_volume(model::SPUNG,z=SA[1.0];phase=:l)
    lb_v0 = lb_volume(model.model_ref)
    T0 = T_scale(model.model_ref)
    f,h = shape_factors(model,lb_v0,T0,z) #h normaly should be independent of temperature
    return lb_v0*h
end

function T_scale(model::SPUNG,z=SA[1.0])
    lb_v0 = lb_volume(model.model_ref)
    T0 = T_scale(model.model_ref)
    f,h = shape_factors(model,lb_v0,T0,z) #h normaly should be independent of temperature
    return T0*f
end

function p_scale(model::SPUNG,z=SA[1.0])
    lb_v0 = lb_volume(model.model_ref)
    T0 = T_scale(model.model_ref)
    p0 = p_scale(model.model_ref)
    f,h = shape_factors(model,lb_v0,T0,z) #h normaly should be independent of temperature
    return p0*f/h 
end

export SPUNG
