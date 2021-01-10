macro comps()
    return :($(esc(:(model.icomponents))))
end

macro groups(component)
    return :($(esc(:(model.igroups[$(component)]))))
end

macro sites(component)
    return :($(esc(:(model.isites[$(component)]))))
end

macro f(func, args...)
    args = [esc(arg) for arg ∈ args]
    return :($(func)($(esc(:model)),$(esc(:V)),$(esc(:T)),$(esc(:z)),$(args...)))
end

macro newmodel(name, parent, paramstype)
    esc(quote struct $name <: $parent
        modelname::String
        components::Array{String,1}
        allcomponentsites::Array{Array{String,1},1}
        allcomponentnsites::Array{Array{Int,1},1}
        params::$paramstype
        icomponents::UnitRange{Int}
        isites::Array{UnitRange{Int},1}
        references::Array{String,1}
        function $name(params::$paramstype, sites::SiteParam; references::Array{String,1}=String[])
            arbitraryparam = getfield(params, first(fieldnames($paramstype)))
            modelname = arbitraryparam.modelname
            components = arbitraryparam.components
            allcomponentsites = sites.allcomponentsites
            allcomponentnsites = sites.allcomponentnsites
            icomponents = 1:length(components)
            isites = [1:length(componentsites) for componentsites ∈ allcomponentsites]
            return new(modelname, components, allcomponentsites, allcomponentnsites, params, icomponents, isites, references)
        end
        function $name(params::$paramstype; references::Array{String,1}=String[])
            arbitraryparam = getfield(params, first(fieldnames($paramstype)))
            sites = SiteParam(arbitraryparam.components, [String[] for _ ∈ 1:length(arbitraryparam.components)], [Int[] for _ ∈ 1:length(arbitraryparam.components)], arbitraryparam.modelname)
            return $name(params, sites; references=references)
        end
    end
    end)
end

macro newmodelgc(name, parent, paramstype)
    esc(quote struct $name <: $parent
        modelname::String
        components::Array{String,1}
        allcomponentgroups::Array{Array{String,1},1}
        allcomponentngroups::Array{Array{Int,1},1}
        flattenedgroups::Array{String,1}
        allcomponentsites::Array{Array{String,1},1}
        allcomponentnsites::Array{Array{Int,1},1}
        params::$paramstype
        icomponents::UnitRange{Int}
        igroups::Array{Array{Int,1},1}
        isites::Array{UnitRange{Int},1}
        references::Array{String,1}
        function $name(params::$paramstype, groups::GCParam, sites::SiteParam; references::Array{String,1}=String[])
            modelname = groups.modelname
            components = groups.components
            allcomponentgroups = groups.allcomponentgroups
            allcomponentngroups = groups.allcomponentngroups
            flattenedgroups = groups.flattenedgroups
            allcomponentsites = sites.allcomponentsites
            allcomponentnsites = sites.allcomponentnsites
            icomponents = 1:length(components)
            igroups = [[findfirst(isequal(group), flattenedgroups) for group ∈ componentgroups] for componentgroups ∈ allcomponentgroups]
            isites = [1:length(componentsites) for componentsites ∈ allcomponentsites]
            return new(modelname, components, allcomponentgroups, allcomponentngroups, flattenedgroups, allcomponentsites, allcomponentnsites, params, icomponents, igroups, isites, references)
        end
        function $name(params::$paramstype, groups::GCParam; references::Array{String,1}=String[])
            sites = SiteParam(groups.components, [String[] for _ ∈ 1:length(groups.components)], [Int[] for _ ∈ 1:length(groups.components)], groups.modelname)
            return $name(params, groups, sites; references=references)
        end
    end
    end)
end
