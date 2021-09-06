struct PCSAFTParam <: EoSParam
    Mw::SingleParam{Float64}
    segment::SingleParam{Float64}
    sigma::PairParam{Float64}
    epsilon::PairParam{Float64}
    epsilon_assoc::AssocParam{Float64}
    bondvol::AssocParam{Float64}
end

abstract type PCSAFTModel <: SAFTModel end
@newmodel PCSAFT PCSAFTModel PCSAFTParam

export PCSAFT
function PCSAFT(components; idealmodel=BasicIdeal, userlocations=String[], ideal_userlocations=String[], verbose=false)
    params,sites = getparams(components, ["SAFT/PCSAFT","properties/molarmass.csv"]; userlocations=userlocations, verbose=verbose)
    
    segment = params["m"]
    k = params["k"]
    Mw = params["Mw"]
    params["sigma"].values .*= 1E-10
    sigma = sigma_LorentzBerthelot(params["sigma"])
    epsilon = epsilon_LorentzBerthelot(params["epsilon"], k)
    epsilon_assoc = params["epsilon_assoc"]
    bondvol = params["bondvol"]

    packagedparams = PCSAFTParam(Mw, segment, sigma, epsilon, epsilon_assoc, bondvol)
    references = ["10.1021/ie0003887", "10.1021/ie010954d"]

    model = PCSAFT(packagedparams, sites, idealmodel; ideal_userlocations=ideal_userlocations, references=references, verbose=verbose)
    return model
end

function a_res(model::PCSAFTModel, V, T, z)
    return @f(a_hc) + @f(a_disp) + @f(a_assoc)
end

function a_hc(model::PCSAFTModel, V, T, z)
    Σz = sum(z)
    m = model.params.segment.values
    m̄ = dot(z, m)/Σz
    return m̄*@f(a_hs) - ∑(z[i]*(m[i]-1)*log(@f(g_hs,i,i)) for i ∈ @comps)/Σz
end

function a_disp(model::PCSAFTModel, V, T, z)
    Σz = sum(z)
    m = model.params.segment.values
    m̄ = dot(z, m)/Σz
    return -2*π*N_A*Σz/V*@f(I,1)*@f(m2ϵσ3, 1) - π*m̄*N_A*Σz/V*@f(C1)*@f(I,2)*@f(m2ϵσ3,2)
end

function d(model::PCSAFTModel, V, T, z, i)
    ϵii = model.params.epsilon.diagvalues[i]
    σii = model.params.sigma.diagvalues[i]
    return σii * (1 - 0.12exp(-3ϵii/T))
end

function ζ(model::PCSAFTModel, V, T, z, n)
    ∑z = ∑(z)
    x = z * (one(∑z)/∑z)
    m = model.params.segment.values
    res = N_A*∑z*π/6/V * ∑((x[i]*m[i]*@f(d,i)^n for i ∈ @comps))
end

function g_hs(model::PCSAFTModel, V, T, z, i, j)
    di = @f(d,i)
    dj = @f(d,j)
    ζ2 = @f(ζ,2)
    ζ3 = @f(ζ,3)
    return 1/(1-ζ3) + di*dj/(di+dj)*3ζ2/(1-ζ3)^2 + (di*dj/(di+dj))^2*2ζ2^2/(1-ζ3)^3
end

function a_hs(model::PCSAFTModel, V, T, z)
    ζ0 = @f(ζ,0)
    ζ1 = @f(ζ,1)
    ζ2 = @f(ζ,2)
    ζ3 = @f(ζ,3)
    return 1/ζ0 * (3ζ1*ζ2/(1-ζ3) + ζ2^3/(ζ3*(1-ζ3)^2) + (ζ2^3/ζ3^2-ζ0)*log(1-ζ3))
end

function C1(model::PCSAFTModel, V, T, z)
    m = model.params.segment.values
    m̄ = dot(z, m)/sum(z)
    η = @f(ζ,3)
    return (1 + m̄*(8η-2η^2)/(1-η)^4 + (1-m̄)*(20η-27η^2+12η^3-2η^4)/((1-η)*(2-η))^2)^-1
end

function m2ϵσ3(model::PCSAFTModel, V, T, z, n = 1)
    m = model.params.segment.values
    σ = model.params.sigma.values
    ϵ = model.params.epsilon.values
    return ∑(z[i]*z[j]*m[i]*m[j] * (ϵ[i,j]*(1)/T)^n * σ[i,j]^3 for i ∈ @comps, j ∈ @comps)/(sum(z)^2)
end

function I(model::PCSAFTModel, V, T, z, n)
    m = model.params.segment.values
    m̄ = dot(z, m)/sum(z)
    σ = model.params.sigma.values
    η = @f(ζ,3)
    if n == 1
        corr = PCSAFT_consts.corr1
    elseif n == 2
       
        corr = PCSAFT_consts.corr2
    end

    res = zero(η)
    @inbounds for i = 1:7
        ii = i-1 
        corr1,corr2,corr3 = corr[i]
        ki = corr1 + (m̄-1)/m̄*corr2 + (m̄-1)/m̄*(m̄-2)/m̄*corr3
        res +=ki*η^ii
    end
    return res
    end

function a_assoc(model::PCSAFTModel, V, T, z)
    if iszero(length(model.sites.flattenedsites))
        return zero(V+T+first(z))
    end
    X_ = @f(X)
    n = model.sites.n_sites
    res =  ∑(z[i]*∑(n[i][a] * (log(X_[i][a]) - X_[i][a]/2 + 0.5) for a ∈ @sites(i)) for i ∈ @comps)/sum(z)
    return res
end

function X(model::PCSAFTModel, V, T, z)
    _1 = one(V+T+first(z))
    Σz = ∑(z)
    ρ = N_A* Σz/V
    itermax = 100
    dampingfactor = 0.5
    error = 1.
    tol = model.absolutetolerance
    iter = 1
    X_ = [[_1 for a ∈ @sites(i)] for i ∈ @comps]
    X_old = deepcopy(X_)
    while error > tol
        iter > itermax && error("X has failed to converge after $itermax iterations")
        for i ∈ @comps, a ∈ @sites(i)
            rhs = 1/(1+∑(ρ*z[j]*∑(X_old[j][b]*@f(Δ,i,j,a,b) for b ∈ @sites(j)) for j ∈ @comps)/Σz)
            X_[i][a] = (1-dampingfactor)*X_old[i][a] + dampingfactor*rhs
        end
        error = sqrt(∑(∑((X_[i][a] - X_old[i][a])^2 for a ∈ @sites(i)) for i ∈ @comps))
        for i = 1:length(X_)
            X_old[i] .= X_[i]
        end
        iter += 1
    end
    return X_
end

function Δ(model::PCSAFTModel, V, T, z, i, j, a, b)
    ϵ_assoc = model.params.epsilon_assoc.values
    κ = model.params.bondvol.values
    σ = model.params.sigma.values
    gij = @f(g_hs,i,j)
    return gij*σ[i,j]^3*(exp(ϵ_assoc[i,j][a,b]/T)-1)*κ[i,j][a,b]
end

const PCSAFTconsts = (
    corr1 =
    [0.9105631445 -0.3084016918 -0.0906148351;
    0.6361281449 0.1860531159 0.4527842806;
    2.6861347891 -2.5030047259 0.5962700728;
    -26.547362491 21.419793629 -1.7241829131;
    97.759208784 -65.255885330 -4.1302112531;
    -159.59154087 83.318680481 13.776631870;
    91.297774084 -33.746922930 -8.6728470368],

    corr2 =
    [0.7240946941 -0.5755498075 0.0976883116;
    2.2382791861 0.6995095521 -0.2557574982;
    -4.0025849485 3.8925673390 -9.1558561530;
    -21.003576815 -17.215471648 20.642075974;
    26.855641363 192.67226447 -38.804430052;
    206.55133841 -161.82646165 93.626774077;
    -355.60235612 -165.20769346 -29.666905585]
)

const PCSAFT_consts = (
    corr1 =
    [(0.9105631445,-0.3084016918, -0.0906148351),
    (0.6361281449, 0.1860531159, 0.4527842806),
    (2.6861347891, -2.5030047259, 0.5962700728),
    (-26.547362491, 21.419793629, -1.7241829131),
    (97.759208784, -65.255885330, -4.1302112531),
    (-159.59154087, 83.318680481, 13.776631870),
    (91.297774084, -33.746922930, -8.6728470368)],

    corr2 =
    [(0.7240946941, -0.5755498075, 0.0976883116),
    (2.2382791861, 0.6995095521, -0.2557574982),
    (-4.0025849485, 3.8925673390, -9.1558561530),
    (-21.003576815, -17.215471648, 20.642075974),
    (26.855641363, 192.67226447, -38.804430052),
    (206.55133841, -161.82646165, 93.626774077),
    (-355.60235612, -165.20769346, -29.666905585)]
)
