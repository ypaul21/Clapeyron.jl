function vt_flash(model::EoSModel,V,T,z;kwargs...)
    method = init_preferred_method(vt_flash,model,kwargs)
    return vt_flash(model,V,T,z,method)
end

function init_preferred_method(method::typeof(vt_flash),model::EoSModel,kwargs) 
    GeneralizedXYFlash(;kwargs...)
end

function vt_flash(model,V,T,z,method::FlashMethod)
    check_arraysize(model,z)
    if supports_reduction(method)
        model_r,idx_r = index_reduction(model,z)
        z_r = z[idx_r]
        method_r = index_reduction(method,idx_r)
    else
        model_r,idx_r = model,trues(length(model))
        method_r,z_r = method,z
    end
    if length(model_r) == 1
        if hasfield(typeof(method),:T0)
            T0 == method.T0
        else
            T0 = nothing
        end
        result1 = vt_flash_pure(model_r,V,T,z_r,T0)
        return index_expansion(result1,idx_r)
    end
    
    result = vt_flash_impl(model_r,V,T,z_r,method_r)
    if !issorted(result.volumes)
        #this is in case we catch a bad result.
        result = FlashResult(result)
    end
    ∑β = sum(result.fractions)
    result.fractions ./= ∑β
    result.fractions .*= sum(z)
    return index_expansion(result,idx_r)
end

function vt_flash_impl(model,V,T,z,method::GeneralizedXYFlash)
    flash0 = tx_flash_x0(model,V,T,z,volume,method)
    isone(numphases(flash0)) && return flash0
    spec = FlashSpecifications(volume,V,temperature,T)
    return xy_flash(model,spec,z,flash0,method)
end

function vt_flash_pure(model,V,T,z,T0 = nothing)
    tx_flash_pure(model,V,T,z,volume,T0)
end