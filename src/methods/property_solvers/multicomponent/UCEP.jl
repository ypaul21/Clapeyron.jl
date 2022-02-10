


"""
    UCEP_mix(model::EoSModel;v0=x0_UCEP_mix(model))

Calculates the Upper Critical End Point of a binary mixture.

returns:
- UCEP Temperature [`K`]
- UCEP Pressure [`Pa`]
- liquid volume at UCEP Point [`m^3`]
- vapour volume at UCEP Point [`m^3`]
- liquid molar composition at UCEP Point
- vapour molar composition at UCEP Point

"""
function UCEP_mix(model::EoSModel;v0=nothing)
    if v0 === nothing
        v0 = x0_UCEP_mix(model)
    end  
    ts = T_scales(model,[0.5,0.5])
    pmix = p_scale(model,[0.5,0.5])
    f! = (F,x) -> Obj_UCEP_mix(model, F, x[1], x[2], exp10(x[3]), exp10(x[4]), x[5],ts,pmix)
    r  = Solvers.nlsolve(f!,v0[1:end],LineSearch(Newton()))
    sol = Solvers.x_sol(r)
    x = FractionVector(sol[1])
    y = FractionVector(sol[2])
    V_l = exp10(sol[3])
    V_v = exp10(sol[4])
    T = sol[5]
    p = pressure(model, V_l, T, x)
    return (T, p, V_l, V_v, x, y)
end



function Obj_UCEP_mix(model::EoSModel,F,x,y,V_l,V_v,T,ts,ps)
    y    = FractionVector(y)
    x    = FractionVector(x)
    μ_l = VT_chemical_potential(model,V_l,T,x)
    μ_v = VT_chemical_potential(model,V_v,T,y)
    p_l = pressure(model,V_l,T,x)
    p_v = pressure(model,V_v,T,y)
    L,detM = LdetM(model,V_l,T,x) 
    for i in 1:length(x)
        F[i] = (μ_l[i]-μ_v[i])/(R̄*ts[i])
    end
    F[3] = (p_l-p_v)/ps
    F[4] = L
    F[5] = detM
    return F
end
"""
    x0_UCEP_mix(model::EoSModel)

Initial point for `UCEP_mix(model)`.

Returns a tuple, containing:
- Initial guess for liquid composition
- Initial guess for vapour composition
- Initial guess for liquid volume `[m³]`
- Initial guess for vapour volume `[m³]`
- Initial guess for UCEP Temperature `[K]`

"""
function x0_UCEP_mix(model::EoSModel)
    T0 = T_scale(model,[0.5,0.5])*1.5
    x0 = 0.5
    y0 = 0.75
    v0 = x0_bubble_pressure(model,T0,[x0,1-x0])
    return [x0,y0,v0[1],v0[2],T0]
end

