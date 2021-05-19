#=
function x0_volume(model::EoS,z; phase = "unknown")
    if phase == "unknown" || is_liquid(phase)
        if typeof(model)<:SAFTgammaMie
            x0 = [log10(π/6*N_A*sum(z[i]*sum(model.group_multiplicities[i][k]*model.params.segment[k]*model.params.shapefactor[k]*model.params.sigma[k]^3 for k in @groups(i)) for i in @comps)/0.8)]
        elseif typeof(model)<:SAFT
            x0 = [log10(π/6*N_A*sum(z[i]*model.params.segment[i]*model.params.sigma[i]^3 for i in model.components)/0.8)]
        elseif typeof(model)<:Cubic
            x0 = [log10(sum(z[i]*z[j]*model.params.b[union(i,j)] for i in model.components for j in model.components)/0.8)]
        end
    elseif is_vapour(phase)
        if typeof(model)<:SAFTgammaMie
            x0 = [log10(π/6*N_A*sum(z[i]*sum(model.group_multiplicities[i][k]*model.params.segment[k]*model.params.shapefactor[k]*model.params.sigma[k]^3 for k in @groups(i)) for i in @comps)/1e-2)]
        elseif typeof(model)<:SAFT
            x0 = [log10(π/6*N_A*sum(z[i]*model.params.segment[i]*model.params.sigma[i]^3 for i in model.components)/1e-2)]
        elseif typeof(model)<:Cubic
            x0 = [log10(sum(z[i]*z[j]*model.params.b[union(i,j)] for i in model.components for j in model.components)/1e-2)]
        end
    elseif is_supercritical(phase)
        if typeof(model)<:SAFTgammaMie
            x0 = [log10(π/6*N_A*sum(z[i]*sum(model.group_multiplicities[i][k]*model.params.segment[k]*model.params.shapefactor[k]*model.params.sigma[k]^3 for k in @groups(i)) for i in @comps)/0.5)]
        elseif typeof(model)<:SAFT
            x0 = [log10(π/6*N_A*sum(z[i]*model.params.segment[i]*model.params.sigma[i]^3 for i in model.components)/0.5)]
        elseif typeof(model)<:Cubic
            x0 = [log10(sum(z[i]*z[j]*model.params.b[union(i,j)] for i in model.components for j in model.components)/0.5)]
        end
    end

    return x0
end
=#

#=x0_volume=#
# function x0_volume(model::SAFTgammaMie,z; phase = "unknown")
#     if phase == "unknown" || is_liquid(phase)
#             x0 = [log10(π/6*N_A*sum(z[i]*sum(model.group_multiplicities[i][k]*model.params.segment[k]*model.params.shapefactor[k]*model.params.sigma[k]^3 for k in @groups(i)) for i in @comps)/0.8)]
#     elseif is_vapour(phase)
#             x0 = [log10(π/6*N_A*sum(z[i]*sum(model.group_multiplicities[i][k]*model.params.segment[k]*model.params.shapefactor[k]*model.params.sigma[k]^3 for k in @groups(i)) for i in @comps)/1e-2)]
#     elseif is_supercritical(phase)
#             x0 = [/0.5)]
#     end
#     return x0
# end


# function x0_volume(model::LJSAFT,z; phase = "unknown")
#     if phase == "unknown" || is_liquid(phase)
#         x0 = [log10(π/6*sum(z[i]*model.params.segment[i]*model.params.b[i] for i in model.components)/0.8)]
#     elseif is_vapour(phase)
#         x0 = [log10(π/6*sum(z[i]*model.params.segment[i]*model.params.b[i] for i in model.components)/1e-2)]
#     elseif is_supercritical(phase)
#         x0 = [/0.5)]
#     end
#     return x0
# end

#=
lb_volume:

SAFTgammaMie:  log10(π/6*N_A*sum(z[i]*sum(model.group_multiplicities[i][k]*model.params.segment[k]*model.params.shapefactor[k]*model.params.sigma[k]^3 for k in @groups(i)) for i in @comps)
LJSAFT: log10(π/6*sum(z[i]*model.params.segment[i]*model.params.b[i] for i in @comps)
=#


function x0_volume_liquid(model,T,z)
    v_lb = lb_volume(model,z)
    return v_lb/0.8
end

function x0_volume_gas(model,p,T,z)
    return volume_virial(model,p,T,z)
end

function x0_volume_sc(model,p,T,z)
    v_sc = lb_volume(model,z)
    return v_sc*2
end

function x0_volume(model::EoSModel,p,T,z; phase = :unknown)
    if phase === :unknown || is_liquid(phase)
        x0val = x0_volume_liquid(model,T,z)
    elseif is_vapour(phase)
        x0val = x0_volume_gas(model,p,T,z)
    elseif is_supercritical(phase)
        x0val = x0_volume_sc(model,p,T,z)
    end
    return x0val
end


#=x0_sat_pure=#

# function x0_sat_pure(model::SAFTVRQMie)
#     x0    = [log10(π/6*N_A*model.params.segment[model.components[1]]*model.params.sigma[model.components[1]]^3/0.2),
#     log10(π/6*N_A*model.params.segment[model.components[1]]*model.params.sigma[model.components[1]]^3/1e-3)]
# end

# function x0_sat_pure(model::LJSAFT)
#     x0    = [log10(π/6*model.params.segment[model.components[1]]*model.params.b[model.components[1]]/0.5),
#     log10(π/6*model.params.segment[model.components[1]]*model.params.b[model.components[1]]/1e-3)]
# end




##=lb_volume=#
#
#lb_volume(model::LJSAFT,z; phase = "unknown") = [log10(π/6*sum(z[i]*model.params.segment[i]*model.params.b[i] for i in model.components)/1)]

function lb_volume(model::SAFTModel, z = SA[1.0]; phase = "unknown")
    seg = model.params.segment.values
    σᵢᵢ = model.params.sigma.diagvalues
    val = π/6*N_A*sum(z[i]*seg[i]*σᵢᵢ[i]^3 for i in @comps)
    return val
end

function lb_volume(model::CubicModel,z = SA[1.0]; phase = "unknown")
    x = z * (1/sum(z))
    b = model.params.b.values
    b̄ = sum(b .* (x * x'))
    return  b̄
end

function lb_volume(model::CPAModel,z = SA[1.0]; phase = "unknown")
    x = z * (1/sum(z))
    b = model.params.b.values
    b̄ = sum(b .* (x * x'))
    return  b̄
end

function lb_volume(model::SAFTgammaMieModel, z = SA[1.0]; phase = "unknown")
    vk  = model.igroups
    seg = model.params.segment.values
    S   = model.params.shapefactor.values
    σᵢᵢ = model.params.sigma.diagvalues


    val = π/6*N_A*sum(z[i]*sum(vk[i][k]*seg[k]*S[k]*σᵢᵢ[k]^3 for k in @groups(i)) for i in @comps)
    return val
end

#=scale_sat_pure=#

# function scale_sat_pure(model::SAFTgammaMie)
#     m̄  = sum(model.group_multiplicities[model.components[1]][k]*model.params.segment[k]*model.params.shapefactor[k] for k in @groups(model.components[1]))
#         σ̄3 = sum(model.group_multiplicities[model.components[1]][k]*model.group_multiplicities[model.components[1]][l]*
#                  model.params.segment[k]*model.params.segment[l]*
#                  model.params.shapefactor[k]*model.params.shapefactor[l]*
#                  model.params.sigma[union(k,l)]^3 for k in @groups(model.components[1]) for l in @groups(model.components[1]))/m̄^2
#         ϵ̄  = T_crit_pure(model)
#         p_scale    = σ̄3*N_A/R̄/ϵ̄
#         μ_scale    = 1/R̄/ϵ̄
#     return p_scale,μ_scale
# end

# function scale_sat_pure(model::LJSAFT)
#     p_scale    = model.params.b[model.components[1]]/N_A/R̄/model.params.T[model.components[1]]
#     μ_scale    = 1/R̄/model.params.T[model.components[1]]
#     return p_scale,μ_scale
# end

function x0_sat_pure(model::SAFTModel,T,z=SA[1.0])
    val = lb_volume(model,z)*one(T)
    x0  = [val/0.5,val/1e-3]
    return log10.(x0)
end

function x0_sat_pure(model::CubicModel,T,z=SA[1.0])
    val = lb_volume(model,z)*one(T)
    x0  = [val/0.9,val/1e-4]
    return log10.(x0)
end

function x0_sat_pure(model::SAFTgammaMieModel,T,z=SA[1.0])
    val = lb_volume(model,z)*one(T)
    #we need to find a better way to scale those values.
    x0  = [val/0.7,val/1e-3]
    return log10.(x0)
end


function x0_sat_pure(model::EoSModel,T,z=SA[1.0])
    val = lb_volume(model,z)*one(T)
    x0  = [val/0.5,val/1e-3]
    return log10.(x0)
end

function scale_sat_pure(model::EoSModel,z=SA[1.0])
    p    = 1/p_scale(model,z)
    μ    = 1/R̄/T_scale(model,z)
    return p,μ
end




#=x0_crit_pure=#

# x0_crit_pure(model::LJSAFT) = [1.5, log10(π/6*model.params.segment[model.components[1]]*model.params.b[model.components[1]]/0.3)]

function x0_crit_pure(model::SAFTModel,z=SA[1.0])
    lb_v = lb_volume(model,z)
    [2, log10(lb_v/0.3)]
end

function x0_crit_pure(model::CubicModel,z=SA[1.0])
    lb_v = lb_volume(model,z)
    [1.0, log10(lb_v/0.3)]
end

function x0_crit_pure(model::CPAModel,z=SA[1.0])
    lb_v = lb_volume(model,z)
    [2.0, log10(lb_v/0.3)]
end

function x0_crit_pure(model::EoSModel,z=SA[1.0])
    lb_v = lb_volume(model,z)
    [1.5, log10(lb_v/0.3)]
end


# x0_crit_pure(model::Cubic) = [1.0, log10(model.params.b[model.components[1]]/0.3)]
# x0_crit_pure(model::CPA) = [2, log10(model.params.b[model.components[1]]/0.3)]



#=T_crit_pure=#
# function T_crit_pure(model::SAFTgammaMie)
#     m̄ = sum(model.group_multiplicities[model.components[1]][k]*model.params.segment[k]*model.params.shapefactor[k] for k in @groups(model.components[1]))
#     return sum(model.group_multiplicities[model.components[1]][k]*model.group_multiplicities[model.components[1]][l]*
#                model.params.segment[k]*model.params.segment[l]*
#                model.params.shapefactor[k]*model.params.shapefactor[l]*
#                model.params.epsilon[union(k,l)] for k in @groups(model.components[1]) for l in @groups(model.components[1]))/m̄^2
# end
# T_crit_pure(model::LJSAFT) = model.params.T[model.components[1]]

T_crit_pure(model::EoSModel,z=SA[1.0]) = T_scale(model,z)

# T_crit_pure(model::Cubic) = model.params.a[model.components[1]]/model.params.b[model.components[1]]/8.314*8/27


#=
 temperature scaling factor,
on critical based EoS, is a function of critical temperature
on SAFT EoS, is a function of ϵ
=#
function T_scale(model::SAFTModel,z=SA[1.0])
    ϵ = model.params.epsilon.diagvalues
    return prod(ϵ)^(1/length(ϵ))
end

#dont use αa, just a, to avoid temperature dependence
function T_scale(model::CubicModel,z=SA[1.0])
    x = z ./ sum(z)
    Ωa,Ωb = ab_consts(model)
    _a = model.params.a.values
    _b = model.params.b.values
    a = dot(x, Symmetric(_a), x)/Ωa
    b = dot(x, Symmetric(_b), x)/Ωb
    return a/b/R̄
end

function T_scale(model::CPAModel,z=SA[1.0])
    x = z ./ sum(z)
    Ωa,Ωb = ab_consts(model)
    _a = model.params.a.values
    _b = model.params.b.values
    a = dot(x, Symmetric(_a), x)/Ωa
    b = dot(x, Symmetric(_b), x)/Ωb
    return a/b/R̄
end

#=
pressure scaling factor
on critical eos, a function of critical pressure
on SAFT, a function of
=#
function p_scale(model::SAFTModel,z=SA[1.0])
    ϵ = model.params.epsilon.diagvalues
    σᵢᵢ = model.params.sigma.diagvalues
    val =  sum(z[i]*σᵢᵢ[i]^3/ϵ[i] for i in 1:length(z))*N_A/R̄
    return 1/val
end

function p_scale(model::SAFTgammaMieModel,z=SA[1.0])
    vk  = model.igroups
    seg = model.params.segment.values
    S   = model.params.shapefactor.values
    σ   = model.params.sigma.values
    m̄  = sum(vk[1][k]*seg[k]*S[k] for k in @groups(model.icomponents[1]))

    σ̄3 = sum(vk[1][k]*vk[1][l]*
                     seg[k]*seg[l]*
                     S[k]*S[l]*
                     σ[k,l]^3 for k in @groups(model.icomponents[1]) for l in @groups(model.icomponents[1]))/m̄^2
    ϵ̄ = T_scale(model,z)
    val    = σ̄3*N_A/R̄/ϵ̄
    return 1/val
end

function p_scale(model::CubicModel,z=SA[1.0])
    sumz = sum(z)
    invsumz = (1/sumz)^2
    Ωa,Ωb = ab_consts(model)
    _a = model.params.a.values
    _b = model.params.b.values
    a = invsumz*dot(z, Symmetric(_a), z)/Ωa
    b = invsumz*dot(z, Symmetric(_b), z)/Ωb
    return a/ (b^2) # Pc mean
end



#=
the following methods are fallbacks,
that require just the definition of T_scale,p_scale and lb_volume
respectively. if possible, each eos should define those
=#
function T_scales(model,z)
    n = length(z)
    x = zeros(n)
    res = zeros(n)
    for i = 1:n
        x[i] = 1.0
        res[i] = T_scale(model,x)
        x[i] = 0.0
    end
    return res
end


function lb_volumes(model,z)
    n = length(z)
    x = zeros(n)
    res = zeros(n)
    for i = 1:n
        x[i] = 1.0
        res[i] = lb_volume(model,x)
        x[i] = 0.0
    end
    return res
end

function p_scales(model,z)
    n = length(z)
    x = zeros(n)
    res = zeros(n)
    for i = 1:n
        x[i] = 1.0
        res[i] = p_scale(model,x)
        x[i] = 0.0
    end
end



function mollerup_k0(tc,pc,p,T)
    return (pc ./ p) .* exp.(5.42 .* (1.0 .- (tc ./ T)))
end

function wilson_k0(ω, tc,pc,p,T)
    return exp.(log.(pc./p).+5.373 .*(1.0 .+ ω).*(1.0 .-tc./T))
end

function bubble_pressure_y0(model::EoSModel,T,x)
    tc = T_scales(model,x)
    pc = p_scales(model,x)
    mollerup(p,t) = mollerup_k0(tc,pc,p,t)
end
