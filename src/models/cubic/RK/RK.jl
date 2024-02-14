const RKParam = ABCubicParam
abstract type RKModel <: ABCubicModel end

struct RK{T <: IdealModel,α,c,M} <: RKModel
    components::Array{String,1}
    alpha::α
    mixing::M
    translation::c
    params::RKParam
    idealmodel::T
    references::Array{String,1}
end

export RK

"""
    RK(components; 
    idealmodel = BasicIdeal,
    alpha = PRAlpha,
    mixing = vdW1fRule,
    activity = nothing,
    translation = NoTranslation,
    userlocations = String[],
    ideal_userlocations = String[],
    alpha_userlocations = String[],
    mixing_userlocations = String[],
    activity_userlocations = String[],
    translation_userlocations = String[],
    verbose = false)

## Input parameters
- `Tc`: Single Parameter (`Float64`) - Critical Temperature `[K]`
- `Pc`: Single Parameter (`Float64`) - Critical Pressure `[Pa]`
- `Mw`: Single Parameter (`Float64`) - Molecular Weight `[g/mol]`
- `k`: Pair Parameter (`Float64`) (optional)
- `l`: Pair Parameter (`Float64`) (optional)

## Model Parameters
- `Tc`: Single Parameter (`Float64`) - Critical Temperature `[K]`
- `Pc`: Single Parameter (`Float64`) - Critical Pressure `[Pa]`
- `Mw`: Single Parameter (`Float64`) - Molecular Weight `[g/mol]`
- `a`: Pair Parameter (`Float64`)
- `b`: Pair Parameter (`Float64`)

## Input models
- `idealmodel`: Ideal Model
- `alpha`: Alpha model
- `mixing`: Mixing model
- `activity`: Activity Model, used in the creation of the mixing model.
- `translation`: Translation Model

## Description
Redlich-Kwong Equation of state.
```
P = RT/(V-Nb) + a•α(T)/(V(V+Nb))
```

## Model Construction Examples
```julia
# Using the default database
model = RK("water") #single input
model = RK(["water","ethanol"]) #multiple components
model = RK(["water","ethanol"], idealmodel = ReidIdeal) #modifying ideal model
model = RK(["water","ethanol"],alpha = Soave2019) #modifying alpha function
model = RK(["water","ethanol"],translation = RackettTranslation) #modifying translation
model = RK(["water","ethanol"],mixing = KayRule) #using another mixing rule
model = RK(["water","ethanol"],mixing = WSRule, activity = NRTL) #using advanced EoS+gᴱ mixing rule

# Passing a prebuilt model

my_alpha = SoaveAlpha(["ethane","butane"],userlocations = Dict(:acentricfactor => [0.1,0.2]))
model =  RK(["ethane","butane"],alpha = my_alpha) #this is efectively now an SRK model

# User-provided parameters, passing files or folders

# Passing files or folders
model = RK(["neon","hydrogen"]; userlocations = ["path/to/my/db","cubic/my_k_values.csv"])

# User-provided parameters, passing parameters directly

model = RK(["neon","hydrogen"];
        userlocations = (;Tc = [44.492,33.19],
                        Pc = [2679000, 1296400],
                        Mw = [20.17, 2.],
                        acentricfactor = [-0.03,-0.21]
                        k = [0. 0.18; 0.18 0.], #k,l can be ommited in single-component models.
                        l = [0. 0.01; 0.01 0.])
                    )
```

## References
1. Redlich, O., & Kwong, J. N. S. (1949). On the thermodynamics of solutions; an equation of state; fugacities of gaseous solutions. Chemical Reviews, 44(1), 233–244. [doi:10.1021/cr60137a013](https://doi.org/10.1021/cr60137a013)
"""
RK

function RK(components;
    idealmodel = BasicIdeal,
    alpha = RKAlpha,
    mixing = vdW1fRule,
    activity = nothing,
    translation = NoTranslation,
    userlocations = String[],
    ideal_userlocations = String[],
    alpha_userlocations = String[],
    mixing_userlocations = String[],
    activity_userlocations = String[],
    translation_userlocations = String[],
    reference_state = nothing,
    verbose = false)
    formatted_components = format_components(components)
    params = getparams(formatted_components, ["properties/critical.csv", "properties/molarmass.csv","SAFT/PCSAFT/PCSAFT_unlike.csv"];
        userlocations=userlocations,
        verbose=verbose,
        ignore_missing_singleparams = __ignored_crit_params(alpha))

    k = get(params,"k",nothing)
    l = get(params,"l",nothing)
    pc = params["Pc"]
    Mw = params["Mw"]
    Tc = params["Tc"]
    acentricfactor = get(params,"acentricfactor",nothing)
    init_mixing = init_model(mixing,components,activity,mixing_userlocations,activity_userlocations,verbose)
    a = PairParam("a",formatted_components,zeros(length(Tc)))
    b = PairParam("b",formatted_components,zeros(length(Tc)))
    init_idealmodel = init_model(idealmodel,components,ideal_userlocations,verbose,reference_state)
    init_alpha = init_alphamodel(alpha,components,acentricfactor,alpha_userlocations,verbose)
    init_translation = init_model(translation,components,translation_userlocations,verbose)
    packagedparams = ABCubicParam(a,b,Tc,pc,Mw)
    references = String["10.1021/cr60137a013"]
    model = RK(formatted_components,init_alpha,init_mixing,init_translation,packagedparams,init_idealmodel,references)
    recombine_cubic!(model,k,l)
    return model
end

function ab_consts(::Type{<:RKModel})
    Ωa =  1/(9*(2^(1/3)-1))
    Ωb = (2^(1/3)-1)/3
    return Ωa,Ωb
end

function cubic_Δ(model::RKModel,z)
    return (0.0,-1.0)
end

#when either Δ1 or Δ2 is equal to zero, requires defining a_res
function a_res(model::RKModel, V, T, z,_data = data(model,V,T,z))
    n,ā,b̄,c̄ = _data
    ρt = (V/n+c̄)^(-1) # translated density
    ρ  = n/V
    RT⁻¹ = 1/(R̄*T)
    return -log1p((c̄-b̄)*ρ) - ā*RT⁻¹*log(b̄*ρt+1)/b̄
    #return -log(V-n*b̄) - ā/(R̄*T*b̄*√(T/T̄c))*log(1+n*b̄/V)
end

crit_pure(model::RKModel) = crit_pure_tp(model)

const CHEB_COEF_L_RK = ([0.9087646203247287,-0.050607640245558645,-0.0008181930140446145,-3.676466126572436e-05,-1.8153449552124612e-06,-4.2148973304623194e-08,8.967883102783958e-09,5.813281384048352e-10,-3.659740635542086e-10,-4.441783746367278e-11,1.5951594212193498e-11,1.7195689316906737e-12,-9.141229440068344e-13,-1.4432899320127035e-14,5.0730253331465747e-14,-6.987466161234579e-15,-1.8318679906315083e-15,1.5126788710517758e-15,-1.734723475976807e-17],
[0.7992363370107459,-0.059445681918569984,-0.0014502328603496623,-7.138164945971748e-05,-3.2189087307149533e-06,-2.486464586864967e-07,-2.4319133658246006e-08,-1.044768556390796e-09,-1.3343107174712543e-10,-1.679712618996021e-11,2.2661039711380226e-13,-1.6726897644758765e-13,-1.0783041126671833e-14,1.4363510381087963e-15,-4.0245584642661925e-16,-7.28583859910259e-16,-4.85722573273506e-17,7.632783294297951e-16,1.734723475976807e-17],
[0.6649628598972782,-0.0760521559973085,-0.002952854145734922,-0.0002237665443423284,-2.2837440018151633e-05,-2.7674939431024392e-06,-3.5249351108129767e-07,-4.743986068561634e-08,-6.602545343792343e-09,-9.417646246179245e-10,-1.37257025190074e-10,-2.031076695718781e-11,-3.0466879019641624e-12,-4.624911564832246e-13,-7.074202335033419e-14,-1.1532441668293814e-14,-1.7277845820728999e-15,3.885780586188048e-16,-2.42861286636753e-17],
[0.5358548459176075,-0.051616853474821534,-0.0019014483430589563,-0.0001522452155350712,-1.6138121186740073e-05,-1.9310179919518333e-06,-2.4780931559553743e-07,-3.3348723178783235e-08,-4.64304821673478e-09,-6.632207977252946e-10,-9.664995193059411e-11,-1.4312308282971031e-11,-2.1472545963518996e-12,-3.260794412263124e-13,-4.987676938128516e-14,-8.229528170033973e-15,-1.2004286453759505e-15,3.469446951953614e-16,-2.0816681711721685e-17],
[0.4485732937583308,-0.03477837295500611,-0.00128197185101675,-0.00010638749439246392,-1.1344133512683485e-05,-1.3588779068512968e-06,-1.7457610832632975e-07,-2.350823218033593e-08,-3.2744967258524493e-09,-4.678942372438044e-10,-6.820332146273422e-11,-1.0101888076041732e-11,-1.5157701482859665e-12,-2.302880108828731e-13,-3.522182545623309e-14,-5.863365348801608e-15,-8.500145032286355e-16,3.0878077872387166e-16,-1.9081958235744878e-17],
[0.3896887959074917,-0.02352363480040463,-0.0008864377996999455,-7.48085471366268e-05,-7.994798847377887e-06,-9.585173595394458e-07,-1.2320694868853121e-07,-1.6596597669088498e-08,-2.312318735792074e-09,-3.3046667832770105e-10,-4.817745855034694e-11,-7.136409518881948e-12,-1.0708690878491467e-12,-1.627031842588167e-13,-2.490022077417109e-14,-4.173744683200198e-15,-5.93275428784068e-16,3.0184188481996443e-16,-1.214306433183765e-17],
[0.3497350745917645,-0.016032294782450524,-0.0006202549276113981,-5.274770800074019e-05,-5.643402364315098e-06,-6.769284142359167e-07,-8.70356488269497e-08,-1.17262048994482e-08,-1.6339517863728226e-09,-2.3353845593487854e-10,-3.404898932246603e-11,-5.04390973432578e-12,-7.568771998034407e-13,-1.1509196373715724e-13,-1.7600504387260685e-14,-3.084338340286763e-15,-4.440892098500626e-16,2.8449465006019636e-16,-3.469446951953614e-18],
[0.3224124351709759,-0.011013634884520482,-0.00043639529901024973,-3.724427437728703e-05,-3.986980578418037e-06,-4.783591501592022e-07,-6.151332742135307e-08,-8.288353604141152e-09,-1.1549863634940927e-09,-1.6508804739912009e-10,-2.4069971710227733e-11,-3.5656928798477594e-12,-5.350303533546708e-13,-8.137934770502397e-14,-1.2437967322753707e-14,-2.203098814490545e-15,-3.157196726277789e-16,2.5673907444456745e-16,0.0],
[0.30358578390142943,-0.0076188736955520135,-0.0003078315487315582,-2.6316343990171603e-05,-2.8179719795684566e-06,-3.381439433146338e-07,-4.348579207141823e-08,-5.8595724626575585e-09,-8.165597746578879e-10,-1.1671785066225127e-10,-1.701778301610446e-11,-2.5209348497590156e-12,-3.7829808730016623e-13,-5.745404152435185e-14,-8.739536871971154e-15,-1.5508427875232655e-15,-1.9081958235744878e-16,2.671474153004283e-16,-3.642919299551295e-17],
[0.29052982685479933,-0.005300319842511697,-0.00021741216965811758,-1.8601576619783206e-05,-1.99216392636134e-06,-2.390659452899957e-07,-3.0745310054813846e-08,-4.142926245120915e-09,-5.773458308655499e-10,-8.252584726697876e-11,-1.2032597140887447e-11,-1.7825428633155838e-12,-2.674006849279209e-13,-4.087702398791748e-14,-6.238065619612598e-15,-1.2177758801357186e-15,-1.457167719820518e-16,3.5735303605122226e-16,1.214306433183765e-17],
[0.281429493928621,-0.00370348431999927,-0.00015364397898540855,-1.3150855018025098e-05,-1.408515659318682e-06,-1.6903172642135367e-07,-2.1738876799376472e-08,-2.929343496937964e-09,-4.0822773633708564e-10,-5.835245481256024e-11,-8.508027615761193e-12,-1.2602800747441023e-12,-1.8902587828328876e-13,-2.898029038966854e-14,-4.3923198411732756e-15,-6.83481049534862e-16,-7.28583859910259e-17,6.210310043996969e-16,3.642919299551295e-17],
[0.27506157036291856,-0.0025962448144549836,-0.00010861133767400138,-9.298192002037131e-06,-9.959154303774975e-07,-1.1951873122098555e-07,-1.5371233163044562e-08,-2.0713067025446286e-09,-2.886544493013732e-10,-4.1260670227694085e-11,-6.015958564642432e-12,-8.912228594004645e-13,-1.3355809513893036e-13,-1.9872992140790302e-14,-3.0531133177191805e-15,-7.181755190543981e-16,-8.500145032286355e-17,2.96637714392034e-16,2.6020852139652106e-18],
[0.2705926697847878,-0.0018244565862979688,-7.678880274703502e-05,-6.574507803324206e-06,-7.04198907909917e-07,-8.451082581135971e-08,-1.0868935571006766e-08,-1.4646166258958093e-09,-2.0410736989440181e-10,-2.9175961993588295e-11,-4.253862886938187e-12,-6.296716620335374e-13,-9.435334458185451e-14,-1.411544492402328e-14,-2.1076890233118206e-15,1.1622647289044608e-16,-1.3877787807814457e-17,2.0816681711721685e-16,-6.071532165918825e-18],
[0.26744977424244143,-0.0012843645575066134,-5.4294008668638594e-05,-4.6487705078930575e-06,-4.979368752669466e-07,-5.97575837875397e-08,-7.685439108034986e-09,-1.035633939716618e-09,-1.4432471884262554e-10,-2.0630074037963908e-11,-3.0077728502275747e-12,-4.461934294264225e-13,-6.715808464896611e-14,-1.0720591081536668e-14,-1.4190038033490282e-15,-4.440892098500626e-16,1.8735013540549517e-16,8.378714388967978e-16,-1.214306433183765e-16],
[0.2652360092057196,-0.0009053089022688643,-3.839029629034568e-05,-3.2871387625962867e-06,-3.5209208477243736e-07,-4.2254784272249823e-08,-5.434404972901041e-09,-7.32300850167511e-10,-1.0205274299235345e-10,-1.458767828754759e-11,-2.126882003850028e-12,-3.151055805172831e-13,-4.713764101271778e-14,-7.455841499748317e-15,-1.1188966420050406e-15,-4.2500725161431774e-16,-1.9081958235744878e-17,4.0072112295064244e-16,1.734723475976807e-17],
[0.2636749575957985,-0.0006387083278175353,-2.7145557162256556e-05,-2.324344536088624e-06,-2.4896583217942636e-07,-2.9878569729469007e-08,-3.8426970945526495e-09,-5.178131846916845e-10,-7.21618736987395e-11,-1.0313611076284701e-11,-1.5037519840443991e-12,-2.2252685805135286e-13,-3.314015728506092e-14,-3.944761184371259e-15,-7.4593109467002705e-16,-1.061650767297806e-15,-7.28583859910259e-17,3.1051550219984847e-16,-8.673617379884035e-19],
[0.26257329202514273,-0.0004509124480249222,-1.9194637512940163e-05,-1.6435549859191573e-06,-1.7604512139550443e-07,-2.11273140666296e-08,-2.7171948059623308e-09,-3.6615067228185083e-10,-5.1026160727274394e-11,-7.291847681223373e-12,-1.0631721197862376e-12,-1.570323732158485e-13,-2.3507237822961713e-14,-4.498137973207861e-15,-5.672545766444159e-16,3.8163916471489756e-16,-1.1622647289044608e-16,-1.5334955527634975e-15,-1.6046192152785466e-16],
[0.261795381349911,-0.00031848138784485065,-1.3572598294689722e-05,-1.1621671743376055e-06,-1.2448259038670695e-07,-1.493926039992932e-08,-1.921346123867629e-09,-2.5890705378228684e-10,-3.6080990625797504e-11,-5.1575497345135446e-12,-7.518204808709683e-13,-1.1052617154838629e-13,-1.6753959330984003e-14,-5.920611223508843e-15,-8.795048023202412e-16,-2.8033131371785203e-15,-1.5265566588595902e-16,1.177877240188252e-15,8.153200337090993e-17],
[0.2612458584790565,-0.00022501920934304802,-9.597255076118291e-06,-8.217756914414426e-07,-8.802244447041196e-08,-1.0563647579939217e-08,-1.358598145442813e-09,-1.83071792100864e-10,-2.551181141341452e-11,-3.646784263455771e-12,-5.306519113013053e-13,-7.894206122127656e-14,-1.3612375115990005e-14,-1.3461454173580023e-15,7.650130529057719e-16,-1.474514954580286e-16,6.262351748276274e-16,-3.9517000782751666e-15,-1.343543332144037e-15],
[0.2608575589284934,-0.0001590219652277923,-6.786276649781153e-06,-5.8108295704698e-07,-6.22412551117002e-08,-7.469627607420426e-09,-9.606729699934702e-10,-1.2945542207654093e-10,-1.804062108035076e-11,-2.578886756920973e-12,-3.76219888575946e-13,-5.927897062107945e-14,-9.209646933960869e-15,-9.674552825522653e-15,-4.822531263215524e-16,4.628242233906121e-15,7.632783294297951e-16,4.489464355827977e-15,3.608224830031759e-16],
[0.2605831256798172,-0.00011240016208966669,-4.79861958626801e-06,-4.108876166461495e-07,-4.401120835134664e-08,-5.281827013953189e-09,-6.792985635800664e-10,-9.153824573338021e-11,-1.2757267464635902e-11,-1.8317153349656223e-12,-2.673868071401131e-13,-5.0074527857546514e-14,-6.376843497690743e-15,2.373101715136272e-15,4.822531263215524e-16,3.7539416020138106e-15,-1.0755285551056204e-16,-4.803449304979779e-15,-4.3281350725621337e-16],
[0.26026738607765193,-0.00023151436012637007,-4.663534061761544e-05,-2.0242818597892773e-05,-1.1449013493607454e-05,-7.455972011261372e-06,-5.310402070114739e-06,-4.027966733814725e-06,-3.2034581341506413e-06,-2.6450585844037255e-06,-2.2524335433048853e-06,-1.969036098618121e-06,-1.7611299692795596e-06,-1.607700890126143e-06,-1.4952529688184307e-06,-1.414974146275727e-06,-1.3611294727017276e-06,-1.3301312448230768e-06,-6.600040751833167e-07])

chebyshev_coef_l(model::RKModel) = CHEB_COEF_L_RK
chebyshev_Trange_l(model::RKModel) = (0.020267685653535945,0.06586997837399182,0.1114722710944477,0.15707456381490356,0.1798757101751315,0.19127628335524546,0.19697656994530244,0.19982671324033094,0.20125178488784518,0.2019643207116023,0.20232058862348087,0.20249872257942014,0.2025877895573898,0.20263232304637463,0.20265458979086703,0.20266572316311324,0.20267128984923632,0.20267407319229788,0.20267546486382865,0.20267616069959404,0.20267650861747674,0.20267668257641808,0.20267685653535944)


const CHEB_COEF_V_RK = ([1.147791769634607e-11,1.8883400952870644e-11,1.0781377594967363e-11,4.424169417038083e-12,1.340057930622405e-12,3.0392837657011474e-13,5.1611350145284516e-14,6.429591419258869e-15,5.526656713208219e-16,2.6838869757647005e-17,-7.144055741834159e-20,-9.585729650435423e-20,-3.2421103573493004e-21,2.5635377164771786e-22,1.5013435876237105e-23,-8.328903877850174e-25,-5.167748904012326e-26,3.2366962941218114e-27,-1.163964265653603e-27],
[1.4330520293934438e-09,2.0653900997266785e-09,8.861422337329386e-10,2.484390807063457e-10,4.771871094962088e-11,6.3539423424235915e-12,5.695290133656733e-13,3.041636735038608e-14,4.618066754497384e-16,-4.882537910833838e-17,-2.2626929302787274e-18,7.735800358738017e-20,5.494308515884454e-21,-1.9098299796206625e-22,-1.0551803468236254e-23,9.358288473073072e-25,8.243533350687156e-25,1.8130926610110065e-24,6.625926732877158e-25],
[3.216589566610815e-07,4.855615137205529e-07,2.2313237230527845e-07,6.564126429055284e-08,1.2457975717948848e-08,1.4356027589196413e-09,7.456366737354983e-11,-2.6879017209951425e-12,-4.739337426481451e-13,7.722229482544613e-15,2.5616030504147617e-15,-7.662739526444225e-17,-1.1817890172164798e-17,8.184648689324978e-19,2.707129500338779e-20,-5.7904193842360474e-21,1.4315394475995835e-22,8.504450672810816e-23,-2.4996364135586805e-23],
[5.7939913777071644e-05,8.379335832127809e-05,3.377936608738917e-05,7.637987852998414e-06,8.302628510075688e-07,1.4927574593586637e-09,-6.691320992196064e-09,1.2878253547929106e-10,7.094465860674562e-11,-5.0660797908678295e-12,-4.994533593074178e-13,9.317691594005963e-14,-1.768263545468354e-15,-8.453117101310306e-16,9.273541699666695e-17,-5.0145673966434665e-19,-7.953374676921863e-19,8.536503921547246e-20,-1.495542547495874e-20],
[0.0031631846269420103,0.003982211787182989,0.0011297528144810199,0.00012650572507611033,-5.2717473574936e-07,-1.4260061216196563e-07,1.3326423389224526e-07,-8.09481292393274e-09,-4.78981067790718e-10,3.091240719345703e-10,-4.0682128710022676e-11,1.5201595949913826e-12,6.298290797186528e-13,-1.5266145620318645e-13,1.6769319426449512e-14,2.1858532236844475e-17,-3.670075175840268e-16,7.421194464073827e-17,-6.8516495103400354e-18],
[0.02527054526018694,0.01958885465789757,0.0029253705571776014,0.00022445447854577003,2.1127284490167714e-05,2.925789705565077e-06,3.436409961475248e-07,4.7523257336800576e-08,6.65800338960857e-09,9.327933602398461e-10,1.3817232769183438e-10,2.024545819618445e-11,3.049514308607615e-12,4.623933614472664e-13,7.069220426050848e-14,1.0914066959216706e-14,1.6929274822274909e-15,2.8888566885876266e-16,4.003416521902725e-17],
[0.06954212908713689,0.023253617876107757,0.0018817558019952863,0.00015099612130791773,1.6119981407421245e-05,1.9323965214038813e-06,2.477548253195308e-07,3.3350354563480705e-08,4.6430122303299515e-09,6.632208966045328e-10,9.665026374713892e-11,1.4312801378119078e-11,2.1476913129869768e-12,3.25807089640584e-13,4.9888478764748e-14,7.651431571664702e-15,1.1947907940790259e-15,2.40692882291782e-16,2.3635607360183997e-17],
[0.11419041805449423,0.020521580327712844,0.0012740772212567466,0.0001062161198180098,1.1343925939054193e-05,1.3588982222490983e-06,1.7457566693114546e-07,2.3508239255405627e-08,3.274496806517091e-09,4.678940377506047e-10,6.82035469767861e-11,1.0102339971507224e-11,1.5161908187288908e-12,2.3003040444669054e-13,3.523483588230292e-14,5.3507545616504615e-15,8.448103328007051e-16,2.393918396847994e-16,2.0816681711721685e-17],
[0.15165818890499952,0.016369427376503674,0.0008840770152282037,7.478700795372859e-05,7.99480167670754e-06,9.58517782677533e-07,1.2320694370814012e-07,1.659659774975314e-08,2.3123189231422092e-09,3.3046648403867174e-10,4.817769534010141e-11,7.1368752921352474e-12,1.0712420533964817e-12,1.6248200701562965e-13,2.4886342986363275e-14,3.694961003830599e-15,5.924080670460796e-16,2.1337098754514727e-16,1.7780915628762273e-17],
[0.1808715039787546,0.0124478423701576,0.0006196163285778558,5.2745024120972456e-05,5.643402899935257e-06,6.769284248055868e-07,8.703564869511071e-08,1.1726204920264882e-08,1.6339519685187875e-09,2.335382477680614e-10,3.40492009587301e-11,5.044305251278303e-12,7.572241444986361e-13,1.1483695938618865e-13,1.7590096046404824e-14,2.5847379792054426e-15,4.198030811863873e-16,2.3765711620882257e-16,3.469446951953614e-18],
[0.202815039861727,0.00921945969494375,0.00043622961797678585,3.7243939807177565e-05,3.986980621388872e-06,4.783591503864509e-07,6.151332734329051e-08,8.28835361281477e-09,1.154986554313675e-09,1.6508787566149596e-10,2.4070209367343942e-11,3.5661248259932776e-12,5.35387706390722e-13,8.11677114409548e-14,1.2444906216657614e-14,1.7468665403086447e-15,2.8102520310824275e-16,2.3071822230491534e-16,6.938893903907228e-18],
[0.2189498032807097,0.006721284919406446,0.0003077893762559753,2.6316302236593578e-05,2.817971982508813e-06,3.381439431966726e-07,4.3485791989886224e-08,5.859572497352028e-09,8.165599741510876e-10,1.16717652903775e-10,1.701801546905024e-11,2.5213789389688657e-12,3.7863462365450573e-13,5.725107887766256e-14,8.748210489351038e-15,1.0911410663894117e-15,1.682681771697503e-16,2.237793284010081e-16,3.209238430557093e-17],
[0.23065921854671798,0.004851398420450417,0.00021740153263174328,1.8601571405155864e-05,1.9921639265139957e-06,2.3906594516856505e-07,3.074531000103742e-08,4.142926290223725e-09,5.773459800517688e-10,8.252565471267292e-11,1.2032836532727131e-11,1.7829522580559143e-12,2.677597726874481e-13,4.066538772384831e-14,6.241535066564552e-15,7.580741590018647e-16,1.3530843112619095e-16,1.5612511283791264e-16,-9.540979117872439e-18],
[0.23908612939229767,0.0034789916333673055,0.00015364130799821497,1.3150854366506329e-05,1.4085156592926612e-06,1.6903172630859664e-07,2.173887675774311e-08,2.9293435507143917e-09,4.0822789593164543e-10,5.835226572770136e-11,8.508240986748739e-12,1.2606964083783367e-12,1.8938843548976791e-13,2.8749572167363624e-14,4.397524011601206e-15,2.5847379792054426e-16,4.5102810375396984e-17,-9.540979117872439e-17,-2.862293735361732e-17],
[0.24511730390452105,0.002483990450115506,0.0001086106684613404,9.298191920630028e-06,9.959154303011697e-07,1.1951873112384104e-07,1.537123309712507e-08,2.0713067459127155e-09,2.8865461930427383e-10,4.126050022479344e-11,6.016184078694309e-12,8.916461319286029e-13,1.339206523454095e-13,1.9670029494101016e-14,3.0687258290029718e-15,2.6020852139652106e-16,7.112366251504909e-17,2.0296264668928643e-16,-8.673617379884035e-19],
[0.24941782042567914,0.0017683273954829087,7.678863526062486e-05,6.57450779313444e-06,7.041989078648142e-07,8.45108257072763e-08,1.086893550855672e-08,1.4646166571208319e-09,2.0410752775423813e-10,2.9175781582346794e-11,4.254084931543112e-12,6.301330984781472e-13,9.470722817095378e-14,1.391421700080997e-14,2.114627917215728e-15,-5.585809592645319e-16,-1.214306433183765e-17,2.8622937353617317e-16,3.469446951953614e-18],
[0.2524765215536145,0.0012562994595141843,5.429396677408306e-05,4.648770506642322e-06,4.979368752443952e-07,5.975758369559936e-08,7.685439031707153e-09,1.0356339657374702e-09,1.4432488711080271e-10,2.0629898830892834e-11,3.0079948948324997e-12,4.4662190612498875e-13,6.752064185544526e-14,1.0500281200087613e-14,1.4103301859691442e-15,-1.734723475976807e-17,-1.8908485888147197e-16,-3.4867941867133823e-16,1.2663481374630692e-16],
[0.2546481887862289,0.0008912762275735033,3.839028581376948e-05,3.287138762469652e-06,3.5209208471866094e-07,4.225478415602335e-08,5.434404913920443e-09,7.32300886596704e-10,1.02052919934148e-10,1.4587508284646944e-11,2.1270849664967173e-12,3.1551150581066167e-13,4.750713711310084e-14,7.212980213111564e-15,1.1084883011491797e-15,-6.591949208711867e-17,1.734723475976807e-18,9.367506770274758e-17,-9.540979117872439e-18],
[0.2561881913448183,0.0006316919590384738,2.7145554542667982e-05,2.3243445360920933e-06,2.4896583211350687e-07,2.9878569618446704e-08,3.842697030367881e-09,5.178132280597714e-10,7.216204023219319e-11,1.0313451481724911e-11,1.5039688244788962e-12,2.2294839585601522e-13,3.350444921501605e-14,3.740063814205996e-15,7.407269242420966e-16,5.967448757360216e-16,6.591949208711867e-17,1.8908485888147197e-16,-7.806255641895632e-18],
[0.2572793323524821,0.00044740425577668823,1.919463685796409e-05,1.6435549859573212e-06,1.7604512135387107e-07,2.1127313966015637e-08,2.717194729634498e-09,3.661506878943621e-10,5.1026320321834184e-11,7.29169849500444e-12,1.0634045727320185e-12,1.5748687076655443e-13,2.3859386688585005e-14,4.270889197854899e-15,5.672545766444159e-16,-8.448103328007051e-16,9.71445146547012e-17,2.0105445086571194e-15,1.448494102440634e-16],
[0.2580519807368856,0.0003167272897560327,1.3572598130872846e-05,1.1621671743965861e-06,1.2448259037282916e-07,1.4939260268090337e-08,1.921346066621754e-09,2.589070850073094e-10,3.6081145016186866e-11,5.157388405230279e-12,7.520546685402252e-13,1.1097546492866428e-13,1.7126924878319016e-14,5.67774993687209e-15,8.673617379884035e-16,2.3071822230491534e-15,1.43982048506075e-16,-6.661338147750939e-16,-7.892991815694472e-17],
[0.2585988724599927,0.00022414215980744862,9.597255035097285e-06,8.217756915021579e-07,8.802244440969664e-08,1.0563647474121085e-08,1.358598076053874e-09,1.8307183026478047e-10,2.5511984885762118e-11,3.646605586937746e-12,5.308739559062303e-13,7.936012957898697e-14,1.395758508770939e-14,1.1587952819525071e-15,-7.546047120499111e-16,-3.0704605524789486e-16,-6.401129626354418e-16,4.4374226515486725e-15,1.3426759704060487e-15],
[0.25898585643609884,0.00015858344033724617,6.786276639469957e-06,5.810829570851439e-07,6.224125504578071e-08,7.46962748945923e-09,9.60672904073978e-10,1.2945542901543483e-10,1.8040785879080978e-11,2.578727162361183e-12,3.764540762452029e-13,5.970050842574182e-14,9.553122182204277e-15,9.42648736845797e-15,4.85722573273506e-16,-5.0844745080880216e-15,-7.73686670285656e-16,-3.9881292712706795e-15,-3.625572064791527e-16],
[0.25925963189740076,0.00011218089961373932,4.798619583632965e-06,4.1088761668951757e-07,4.401120831491745e-08,5.281826880379481e-09,6.792984976605743e-10,9.153828042784973e-11,1.2757453080047831e-11,1.8315314542771688e-12,2.676036475746102e-13,5.0482187874401063e-14,6.720318745934151e-15,-2.589942149633373e-15,-4.805184028455756e-16,-4.217112770099618e-15,1.0061396160665481e-16,5.290906601729262e-15,4.3281350725621337e-16],
[0.25957493297459344,0.0002312950976299348,4.6635340614997745e-05,2.024281859794655e-05,1.1449013493564086e-05,7.45597201117984e-06,5.310402070080045e-06,4.027966733847685e-06,3.2034581343189095e-06,2.6450585842215796e-06,2.2524335435182563e-06,1.969036099069149e-06,1.7611299696455862e-06,1.6077008899179762e-06,1.495252968823635e-06,1.4149741458021475e-06,1.361129472705197e-06,1.330131245348698e-06,6.600040751911229e-07])

chebyshev_coef_v(model::RKModel) = CHEB_COEF_V_RK
chebyshev_Trange_v(model::RKModel) = (0.020267685653535945,0.02596797224359293,0.03166825883364991,0.04306883201376388,0.06586997837399182,0.1114722710944477,0.15707456381490356,0.1798757101751315,0.19127628335524546,0.19697656994530244,0.19982671324033094,0.20125178488784518,0.2019643207116023,0.20232058862348087,0.20249872257942014,0.2025877895573898,0.20263232304637463,0.20265458979086703,0.20266572316311324,0.20267128984923632,0.20267407319229788,0.20267546486382865,0.20267616069959404,0.20267650861747674,0.20267668257641808,0.20267685653535944)


const CHEB_COEF_P_RK = ([2.9225473692726744e-13,4.846211782171669e-13,2.824570002653489e-13,1.1955110559200177e-13,3.771710795961595e-14,9.009392578047045e-15,1.6354246747926244e-15,2.2297564454973646e-16,2.1977303354380907e-17,1.4079420760456038e-18,3.645901247485968e-20,-2.3224668904548756e-21,-2.1119277731077838e-22,1.3269831236005753e-24,7.111685153229056e-25,2.0749322108366765e-27,-2.355803670950289e-27,1.89511506527704e-29,-1.203321029253145e-29],
[4.42411833611544e-11,6.486786260066829e-11,2.8834316573742172e-11,8.49036232570771e-12,1.7382606171711586e-12,2.519225194793423e-13,2.5510836572731127e-14,1.688802065855514e-15,5.658115744876586e-17,-7.526264008617426e-19,-1.3472892665001768e-19,-9.922698129701787e-22,2.6795301476514526e-22,2.2940832430296035e-24,-5.789860391087721e-25,1.1649897848303902e-26,2.589988124021626e-26,5.472801416061292e-26,2.268132874691453e-26],
[1.3403609277725988e-08,2.0613817317201788e-08,9.908855672700315e-09,3.1242347660239074e-09,6.566525761190694e-10,8.934527720041564e-11,6.864947922093292e-12,1.0961157691701013e-13,-2.552558822185046e-14,-1.0755067059569204e-15,1.1595758568473605e-16,4.364935105515505e-18,-6.518267188981206e-19,-2.420947626615306e-21,3.3222352209999887e-21,-1.4160136552428753e-22,-1.0359144702519558e-23,3.0017608947725106e-24,-1.1502980393315541e-24],
[3.6266120105336493e-06,5.40530078063072e-06,2.3535707621866956e-06,6.09797026711632e-07,8.755959081374292e-08,4.4744907993962915e-09,-4.067655714209484e-10,-3.591728379118743e-11,4.465366516258063e-12,1.6493814952480902e-13,-5.3226243319885286e-14,1.7922687946152793e-15,4.005718486965302e-16,-4.988684915276682e-17,4.789814152406683e-19,4.147268524381617e-19,-4.591183271914325e-20,4.350970022028926e-22,-1.0629270871306406e-21],
[0.0003094297409003984,0.0004117661948641284,0.0001333784156235523,1.953127675250332e-05,4.4650824621707415e-07,-1.245752753612688e-07,6.1713045177982556e-09,8.457787641085228e-10,-2.019178315130611e-10,1.671511339141844e-11,7.461694230598466e-13,-4.260132318823111e-13,6.171897036311516e-14,-3.212145765269425e-15,-5.836180885150665e-16,1.719715110942645e-16,-2.147440279526465e-17,8.273394312306379e-19,2.34734005507848e-19],
[0.0074771598568592074,0.008244590094104293,0.0017465545422224592,9.786738241774412e-05,-6.516977077872878e-06,3.2958650511163794e-07,1.5363696321735632e-08,-7.20353161771562e-09,1.2861343701874766e-09,-1.6587511786082245e-10,1.5108394000934097e-11,-3.6747905743763146e-13,-2.33340277167645e-13,7.039936125372015e-14,-1.3651690997269761e-14,2.071300487897776e-15,-2.45029690981724e-16,2.15485181781494e-17,1.4026865606531214e-18])

chebyshev_coef_p(model::RKModel) = CHEB_COEF_P_RK
chebyshev_Trange_p(model::RKModel) = (0.020267685653535945,0.02596797224359293,0.03166825883364991,0.04306883201376388,0.06586997837399182,0.1114722710944477,0.20267685653535944)

#for saturation_temperature
chebyshev_prange_T(model::RKModel) = (6.8022572620412635e-16, 1.2274923410316115e-12, 1.484511633296413e-10, 4.7803462965047013e-8, 1.2086876558042846e-5, 0.000874434393370233, 0.01755999378002107)
