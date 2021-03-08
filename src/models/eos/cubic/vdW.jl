struct vdWParam <: EoSParam
    a::PairParam{Float64}
    b::PairParam{Float64}
    Mw::SingleParam{Float64}
end

abstract type vdWModel <: ABCubicModel end
@newmodel vdW vdWModel vdWParam

export vdW
function vdW(components::Array{String,1}; userlocations::Array{String,1}=String[], verbose=false,idealmodel=BasicIdeal)
    params = getparams(components, ["properties/critical.csv", "properties/molarmass.csv","SAFT/PCSAFT/PCSAFT_unlike.csv"]; userlocations=userlocations, verbose=verbose)
    Mw = params["Mw"]
    k = params["k"]
    pc = params["pc"].values
    Tc = params["Tc"].values

    a = epsilon_LorentzBerthelot(SingleParam(params["pc"], @. 27/64*R̄^2*Tc^2/pc/1e6), k)
    b = sigma_LorentzBerthelot(SingleParam(params["pc"], @. 1/8*R̄*Tc/pc/1e6))

    packagedparams = vdWParam(a, b,Mw)
    return vdW(packagedparams,idealmodel)
end


function cubic_ab(model::vdWModel,T,x)
    a = model.params.a.values
    b = model.params.b.values
    ā = sum(a .* (x * x'))
    b̄ = sum(b .* (x * x'))
    return ā,b̄
end

function cubic_abp(model::vdWModel, V, T, z=@SVector [1.0])
    x = z/sum(z)
    n = sum(z)
    v = V/n
    a,b = cubic_ab(model,T,x)
    p = R̄*T/(v-b) - a/(v*v)
    return a,b,p
end

function cubic_poly(model::vdWModel,p,t,x)
    a,b = cubic_ab(model,t,x)
    RT⁻¹ = 1/(R̄*T)
    A = a*p*RT⁻¹*RT⁻¹
    B = b*p*RT⁻¹
    return (-A*B, A, -B-_1, _1)
end



function a_resx(model::vdWModel, v, T, x)
    a,b = cubic_ab(model,T,x)
    RT⁻¹ = 1/(R̄*T)
    ρ = 1/v
    -log(1-b*ρ) - a*ρ*RT⁻¹
    #return -log(V-n*b̄) - ā*n/(R̄*T*V)
end
