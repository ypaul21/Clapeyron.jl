struct ABCubicParam <: EoSParam
    a::PairParam{Float64}
    b::PairParam{Float64}
    Tc::SingleParam{Float64}
    Pc::SingleParam{Float64}
    Mw::SingleParam{Float64}
end

struct ABCCubicParam <: EoSParam
    a::PairParam{Float64}
    b::PairParam{Float64}
    c::SingleParam{Float64}
    Tc::SingleParam{Float64}
    Pc::SingleParam{Float64}
    Vc::SingleParam{Float64}
    Mw::SingleParam{Float64}
end

"""
    ab_premixing(::Type{T},mixing,Tc,pc,kij) where T <: ABCubicModel

given `Tc::SingleParam`, `pc::SingleParam`, `kij::PairParam` and `mixing <: MixingRule`, it will return 
`PairParam`s `a` and `b`, containing values aᵢⱼ and bᵢⱼ. by default, it performs the van der Wals One-Fluid mixing rule. that is:
```
aᵢⱼ = sqrt(aᵢ*aⱼ)*(1-kᵢⱼ)
bᵢⱼ = (bᵢ + bⱼ)/2
```
"""
function ab_premixing end

function ab_premixing(model,mixing,Tc,pc,kij) 
    Ωa, Ωb = ab_consts(model)
    _Tc = Tc.values
    _pc = pc.values
    components = Tc.components
    a = epsilon_LorentzBerthelot(SingleParam("a",components, @. Ωa*R̄^2*_Tc^2/_pc),kij)
    b = sigma_LorentzBerthelot(SingleParam("b",components, @. Ωb*R̄*_Tc/_pc))
    return a,b
end

ab_premixing(model,mixing,Tc,pc,vc,kij) = ab_premixing(model,mixing,Tc,pc,kij) #ignores the Vc unless dispatch

function c_premixing end

function cubic_ab(model::ABCubicModel,V,T,z=SA[1.0],n=sum(z))
    a = model.params.a.values
    b = model.params.b.values
    T = T * float(one(T))
    α = @f(α_function, model.alpha)
    c = @f(translation, model.translation)
    if length(z) > 1
        ā, b̄, c̄ = @f(mixing_rule, model.mixing, α, a, b, c)
    else
        ā = a[1, 1] * α[1]
        b̄ = b[1, 1]
        c̄ = c[1]
    end
    return ā, b̄, c̄
end

function data(model::ABCubicModel, V, T, z)
    n = sum(z)
    ā, b̄, c̄ = cubic_ab(model, V, T, z, n)
    return n, ā, b̄, c̄
end

function a_res(model::ABCubicModel, V, T, z,_data = data(model,V,T,z))
    n,ā,b̄,c̄ = _data
    Δ1,Δ2 = cubic_Δ(model,z)
    ΔΔ = Δ1 - Δ2
    RT⁻¹ = 1/(R̄*T)
    ρt = (V/n+c̄)^(-1) # translated density
    ρ  = n/V
    b̄ρt = b̄*ρt
    a₁ = -log1p((c̄-b̄)*ρ)
    if Δ1 == Δ2
        return a₁ - ā*ρt*RT⁻¹
    else
        l1 = log1p(Δ1*b̄ρt)
        l2 = log1p(Δ2*b̄ρt)
        return a₁ - ā*RT⁻¹*(l1-l2)/(ΔΔ*b̄) 
    end
end

function cubic_poly(model::ABCubicModel,p,T,z)
    a,b,c = cubic_ab(model,p,T,z)
    RT⁻¹ = 1/(R̄*T)
    A = a*p*RT⁻¹*RT⁻¹
    B = b*p*RT⁻¹
    Δ1,Δ2 = cubic_Δ(model,z)
    ∑Δ = Δ1 + Δ2
    Δ1Δ2 = Δ1*Δ2
    k₀ = -B*evalpoly(B,(A,Δ1Δ2,Δ1Δ2))
    k₁ = evalpoly(B,(A,-∑Δ,Δ1Δ2-∑Δ))
    k₂ = (∑Δ - 1)*B - 1
    k₃ = one(A) # important to enable autodiff
    return (k₀,k₁,k₂,k₃),c
end

function cubic_abp(model::ABCubicModel, V, T, z)
    Δ1,Δ2 = cubic_Δ(model,z)
    n = ∑(z)
    a,b,c = cubic_ab(model,V,T,z,n)
    v = V/n+c
    p = R̄*T/(v-b) - a/((v-Δ1*b)*(v-Δ2*b))
    return a,b,p
end

function pure_cubic_zc(model::ABCubicModel)
    Δ1,Δ2 = cubic_Δ(model,SA[1.0])
    _,Ωb = ab_consts(model)
    Ωb = only(Ωb)
    return (1 - (Δ1+Δ2-1)*Ωb)/3
end

function second_virial_coefficient(model::ABCubicModel,T::Real,z = SA[1.0])
    a,b,c = cubic_ab(model,1/sqrt(eps(float(T))),T,z)
    return b-a/(R̄*T)
end

function lb_volume(model::CubicModel, z=SA[1.0])
    V = 1e-5
    T = 0.0
    n = sum(z)
    invn = one(n) / n
    b = model.params.b.values
    c = @f(translation, model.translation)
    b̄ = dot(z, Symmetric(b), z) * invn #b has m3/mol units, result should have m3 units
    c̄ = dot(z, c) * invn
    return b̄ - c̄
end
#dont use αa, just a, to avoid temperature dependence
function T_scale(model::CubicModel, z=SA[1.0])
    n = sum(z)
    invn2 = one(n) / (n * n)
    Ωa, Ωb = ab_consts(model)
    _a = model.params.a.values
    _b = model.params.b.values
    a = dot(z, Symmetric(_a), z) * invn2 / Ωa
    b = dot(z, Symmetric(_b), z) * invn2 / Ωb
    return a / b / R̄
end

function p_scale(model::CubicModel, z=SA[1.0])
    n = sum(z)
    invn2 = (1 / n)^2
    Ωa, Ωb = ab_consts(model)
    _a = model.params.a.values
    _b = model.params.b.values
    a = invn2 * dot(z, Symmetric(_a), z) / Ωa
    b = invn2 * dot(z, Symmetric(_b), z) / Ωb
    return a / (b^2) # Pc mean
end

function x0_crit_pure(model::CubicModel)
    lb_v = lb_volume(model)
    (1.0, log10(lb_v / 0.3))
end

#works with models with a fixed (Tc,Pc) coordinate
function crit_pure_tp(model)
    Tc = model.params.Tc.values[1]
    Pc = model.params.Pc.values[1]
    Vc = volume(model,Pc,Tc,SA[1.])
    return (Tc,Pc,Vc)
end

function volume_impl(model::ABCubicModel,p,T,z=SA[1.0],phase=:unknown,threaded=false)
    lb_v   =lb_volume(model,z)
    RTp = R̄*T/p
    _poly,c̄ = cubic_poly(model,p,T,z)
    sols = Solvers.roots3(_poly)
    function imagfilter(x)
        absx = abs(imag(x))
        return absx < 8 * eps(typeof(absx))
    end
    x1, x2, x3 = sols
    sols = (x1, x2, x3)
    xx = (x1, x2, x3)
    isreal = imagfilter.(xx)
    vvv = extrema(real.(xx))
    
    zl,zg = vvv
    vvl,vvg = RTp*zl,RTp*zg
    @show vvl,vvg
    err() = @error("model $model Failed to converge to a volume root at pressure p = $p [Pa], T = $T [K] and compositions = $z")
    if sum(isreal) == 3 #3 roots
        vg = vvg
        _vl = vvl
        vl = ifelse(_vl > lb_v, _vl, vg) #catch case where solution is unphysical
    elseif sum(isreal) == 1
        i = findfirst(imagfilter, sols)
        vl = real(sols[i]) * nRTp
        vg = real(sols[i]) * nRTp
    elseif sum(isreal) == 0

        V0 = x0_volume(model, p, T, z; phase)
        v = _volume_compress(model, p, T, z, V0)
        isnan(v) && err()
        return v
    end

    function gibbs(v)
        _df, f = ∂f(model, v, T, z)
        dv, dt = _df
        if abs((p + dv) / p) > 0.03
            return one(dv) / zero(dv)
        else
            return f + p * v
        end
    end
    #this catches the supercritical phase as well
    if vl ≈ vg
        return vl - c̄
    end

    if is_liquid(phase)
        return vl - c̄
    elseif is_vapour(phase)
        return vg - c̄
    else
        gg = gibbs(vg - c̄)
        gl = gibbs(vl - c̄)
        return ifelse(gg < gl, vg - c̄, vl - c̄)
    end
end

function ab_consts(model::CubicModel)
    return ab_consts(typeof(model))
end

function x0_sat_pure(model::ABCubicModel, T)
    a, b, c = cubic_ab(model, 1 / sqrt(eps(float(T))), T)
    Tc = model.params.Tc.values[1]
    pc = model.params.Pc.values[1]
    zc = pure_cubic_zc(model)
    vc = zc*R̄*Tc/pc - c
    if Tc < T
        nan = zero(T) / zero(T)
        return (nan, nan)
    end
    B = b-a/(R̄*T)
    pv0 = -0.25*R̄*T/B
    Δ1,Δ2 = cubic_Δ(model,SA[1.0])
    k⁻¹ = (1 + Δ1)*(1 + Δ2)/2
    vl = b + sqrt(k⁻¹*R̄*T*b^3/a) - c
    pc = model.params.Pc.values[1]
    p_vl = pressure(model, vl, T)
    p_low = min(p_vl, pc)
    pl0 = max(zero(b), p_low)
    p0 = 0.5 * (pl0 + pv0)
    vv = volume_virial(B, p0, T) - c
    if p_vl > pc #improves predictions around critical point
        vlc, vvc = vdw_x0_xat_pure(T, Tc, pc, vc)
        vl = 0.5 * (vl + vlc)
        vv = 0.5 * (vv + vvc)
    end
    return (log10(vl), log10(vv))
end
#=
#on the dpdv limit:
dp/dv = 0
p = RT/v-b - a/pol(v)
dpdv = -RT/(v-b)^2 + a/pol^2 * dpol = a*k -RT/(v-b)^2

vdw: pol = v2 -> pol(b) = b^2, dpol(b) = 2b
pr: pol = v2 + 2bv - b2 -> pol(b) = 2b^2, dpol(b) = 2v + 2b = 4b
rk: pol = v*(v+b) -> pol(b) = 2b2, dpol(b) = 2v + b = 3b

in general:
pol(b) = (b + Δ1b)(b + Δ2b) = b^2(1 + Δ1)(1 + Δ2)
dpol(b) = 2*(1 + Δ1)(1 + Δ2)*b

k = dpol(b)/pol(b)^2
k = 2*(1 + Δ1)(1 + Δ2)/((1 + Δ1)(1 + Δ2))^2 b3
k⁻¹ = (1 + Δ1)(1 + Δ2)/2

vdw: k⁻¹ = 0.50b3
pr:  k⁻¹ = 1.00b3
rk:  k⁻¹ = 1.33b3

solving for dpdv = 0
0 = a*k -RT/(v-b)^2
(v-b)^2 = RT/ak
v2 - 2vb + b2 - RT/ak = 0
v = b ± sqrt(b2 +  RT/ak - b2) #v > b
v = b + sqrt(k⁻¹b3RT/a)
on models with translation:
vl = b + sqrt(k⁻¹RTb3/2a) - c

if k⁻¹ not available, use k⁻¹ = 0.5 (vdw, lowest)
=#

function wilson_k_values(model::ABCubicModel, p, T)
    Pc = model.params.Pc.values
    Tc = model.params.Tc.values

    if hasfield(typeof(model.alpha.params), :acentricfactor)
        ω = model.alpha.params.acentricfactor.values
    else
        pure = split_model(model)
        ω = zero(Tc)
        for i in 1:length(Tc)
            ps = first(saturation_pressure(pure[i], 0.7 * Tc[i]))
            ω[i] = -log10(ps / Pc[i]) - 1.0
        end
    end

    return @. Pc / p * exp(5.373 * (1 + ω) * (1 - Tc / T))
end

function vdw_tv_mix(Tc,Vc,z)
    Tm = zero(first(Tc)+first(Vc))
    Vm = zero(eltype(Vc))
    n = sum(z)
    invn2 = (1/n)^2
    for i in 1:length(z)
        zi = z[i]
        Vi = Vc[i]
        Ti = Tc[i]
        zii = zi*zi
        Vm += zii*Vi
        Tm += zii*Ti*Vi
        for j in 1:i-1
            zj = z[j]
            Vj = Vc[j]
            Tj = Tc[j]
            Tij = sqrt(Ti*Tj)
            Vij = 0.5*(Vi+Vj)
            zij = zj*zi
            Vm += 2zij*Vij
            Tm += zij*Tij*Viij
        end
    end
    Tcm = Tm/Vm
    Vcm = Vm*invn2
    return (Tcm,Vcm)
end