module Clapeyron
using LinearAlgebra
using SparseArrays
#for the assoc solver and the sparse packed VofV
import PackedVectorsOfVectors
const PackedVofV = PackedVectorsOfVectors.PackedVectorOfVectors

#for non allocating vectors of zeros and ones
using Roots: Roots

using Scratch 
import LogExpFunctions
using FillArrays: FillArrays
import BlackBoxOptim
using StaticArrays

using DiffResults, ForwardDiff
using Downloads #for bibtex
using StableTasks #for multithreaded volume

#compatibility and raw julia utilities
include("utils/core_utils.jl")

#misc functions
include("utils/misc.jl")

include("modules/solvers/Solvers.jl")
using .Solvers
using .Solvers: log, sqrt, log1p, ^, dnorm, primalval, Newton2
using NLSolvers
using NLSolvers: NEqOptions
#misc functions, useful for EoS, don't depend on models
include("modules/eosfunctions/EoSFunctions.jl")
using .EoSFunctions

include("modules/fractions/Fractions.jl")
import .Fractions
using .Fractions: FractionVector

#Gas constant, Boltzmann Constant
include("base/constants.jl") 

#The Base of Clapeyron: EoSModel and eos(model,V,T,z)
include("base/EoSModel.jl")

#error handlers
include("base/errors.jl")

#type hierarchy
include("models/types.jl")

#show(model<:EoSModel)
include("base/eosshow.jl")


#EoSParam, ClapeyronParam, All Params
include("database/ClapeyronParam.jl")

#recombine options
include("utils/recombine.jl")

#Combining Rules for Clapeyron Params.
include("database/combiningrules.jl")

#general DB files
using Tables,CSV

#used for reading multiparameter json files
using JSON3

#getparams options
include("database/ParamOptions.jl") 
#getparams definition
include("database/database.jl")
#transform Tables.jl tables to Clapeyron csv files
include("database/UserReader.jl")



#macros, used for defining models
include("utils/macros.jl")

#index reduction
include("utils/index_reduction.jl")

#splitting models, useful for methods.
include("utils/split_model.jl")

#exportting models, useful for parameter estimation.
include("utils/export_model.jl")

# Gustavo: acceleration for successive substitution
include("utils/acceleration_ss.jl")

#Clapeyron methods (AD, property solvers, etc)
include("methods/methods.jl")


#=
the dependency chain is the following:
base --> database(params)  -|-> split_model --> methods -|-> models
                            |-> macros ------------------|
=#

#Clapeyron EoS collection


include("models/ideal/ideal.jl")
include("models/ideal/BasicIdeal.jl")
include("models/ideal/MonomerIdeal.jl")
include("models/ideal/ReidIdeal.jl")
include("models/ideal/WalkerIdeal.jl")
include("models/ideal/JobackIdeal.jl")
include("models/ideal/CPLNGEstIdeal.jl")
include("models/ideal/ShomateIdeal.jl")
include("models/ideal/PPDSIdeal.jl")

#AlyLee Ideal uses gerg 2008 terms
include("models/EmpiricHelmholtz/term_functions.jl")
include("models/ideal/AlyLeeIdeal.jl")

#Basic utility EoS
include("models/utility/EoSVectorParam.jl")
include("models/utility/ZeroResidual.jl")
include("models/utility/TPFlashWrapper.jl")

#Empiric Models uses CompositeModel
include("models/CompositeModel/CompositeModel.jl")

#softSAFT2016 uses LJRef. softSAFT uses x0_sat_pure with LJ correlations (from LJRef)
include("models/EmpiricHelmholtz/SingleFluid/SingleFluid.jl")
include("models/EmpiricHelmholtz/SingleFluid/variants/IAPWS95.jl")
include("models/EmpiricHelmholtz/SingleFluid/variants/PropaneRef.jl")
include("models/EmpiricHelmholtz/SingleFluid/variants/Ammonia2023.jl")
include("models/EmpiricHelmholtz/SingleFluid/variants/TholLJ.jl")
include("models/EmpiricHelmholtz/SingleFluid/variants/XiangDeiters.jl")
include("models/EmpiricHelmholtz/LJRef/LJRef.jl")

#multifluid models
include("models/EmpiricHelmholtz/MultiFluid/multifluid.jl")
include("models/EmpiricHelmholtz/MultiFluid/mixing/mixing.jl")
include("models/EmpiricHelmholtz/MultiFluid/departure/departure.jl")
include("models/EmpiricHelmholtz/MultiFluid/variants/GERG2008.jl")
include("models/EmpiricHelmholtz/MultiFluid/variants/EOS_LNG.jl")
include("models/EmpiricHelmholtz/MultiFluid/variants/EOS_CG.jl")
include("models/EmpiricHelmholtz/MultiFluid/variants/TillnerRothFriend.jl")
include("models/EmpiricHelmholtz/MultiFluid/variants/HelmAct.jl")
include("models/EmpiricHelmholtz/MultiFluid/variants/EmpiricIdeal.jl")
include("models/EmpiricHelmholtz/LKP/LKP.jl")
include("models/EmpiricHelmholtz/LKP/variants/LKPmod.jl")
include("models/EmpiricHelmholtz/LKP/variants/LKPSJT.jl")

#cubic models
include("models/cubic/equations.jl")
include("models/cubic/vdW/vdW.jl")
include("models/cubic/RK/RK.jl")
include("models/cubic/PR/PR.jl")
include("models/cubic/KU/KU.jl")
include("models/cubic/RKPR/RKPR.jl")


#SAFT models
include("models/SAFT/association.jl")
include("models/SAFT/equations.jl")
include("models/SAFT/PCSAFT/PCSAFT.jl")
include("models/SAFT/PCSAFT/variants/sPCSAFT.jl")
include("models/SAFT/PCSAFT/variants/gcsPCSAFT.jl")
include("models/SAFT/PCSAFT/variants/PharmaPCSAFT.jl")
include("models/SAFT/PCSAFT/variants/ADPCSAFT.jl")
include("models/SAFT/PCSAFT/variants/PCPSAFT.jl")
include("models/SAFT/PCSAFT/variants/QPCPSAFT.jl")
include("models/SAFT/PCSAFT/variants/gcPCPSAFT/homogcPCPSAFT.jl")
include("models/SAFT/PCSAFT/variants/gcPCPSAFT/heterogcPCPSAFT.jl")
include("models/SAFT/PCSAFT/variants/gcPCPSAFT/gcPCSAFT.jl")
include("models/SAFT/PCSAFT/variants/CPPCSAFT.jl")
include("models/SAFT/ogSAFT/ogSAFT.jl")
include("models/SAFT/CPA/CPA.jl")
include("models/SAFT/CPA/variants/sCPA.jl")
include("models/SAFT/SAFTVRSW/SAFTVRSW.jl")
include("models/SAFT/LJSAFT/LJSAFT.jl")
include("models/SAFT/softSAFT/softSAFT.jl")
include("models/SAFT/softSAFT/variants/softSAFT2016.jl")
include("models/SAFT/softSAFT/variants/solidsoftSAFT.jl")
include("models/SAFT/SAFTVRMie/SAFTVRMie.jl")
include("models/SAFT/SAFTVRMie/variants/SAFTVRQMie.jl")
include("models/SAFT/SAFTVRMie/variants/SAFTVRMie15.jl")
include("models/SAFT/SAFTVRMie/variants/SAFTVRSMie.jl")
include("models/SAFT/SAFTVRMie/variants/SAFTVRMieGV.jl")
include("models/SAFT/SAFTgammaMie/SAFTgammaMie.jl")
include("models/SAFT/SAFTgammaMie/variants/structSAFTgammaMie.jl")
include("models/SAFT/CKSAFT/CKSAFT.jl")
include("models/SAFT/CKSAFT/variants/sCKSAFT.jl")
include("models/SAFT/BACKSAFT/BACKSAFT.jl")
include("models/SAFT/DAPT/DAPT.jl")
include("models/SAFT/COFFEE/COFFEE.jl")
#Activity models
include("models/Activity/Wilson/Wilson.jl")
include("models/Activity/Wilson/variants/tcPRWilsonRes.jl")
include("models/Activity/NRTL/NRTL.jl")
include("models/Activity/NRTL/variants/aspenNRTL.jl")
include("models/Activity/UNIQUAC/UNIQUAC.jl")
include("models/Activity/UNIFAC/utils.jl")
include("models/Activity/UNIFAC/UNIFAC.jl")
include("models/Activity/UNIFAC/variants/ogUNIFAC.jl")
include("models/Activity/UNIFAC/variants/UNIFACFV.jl")
include("models/Activity/UNIFAC/variants/UNIFACFVPoly.jl")
include("models/Activity/UNIFAC/variants/PSRK.jl")
include("models/Activity/UNIFAC/variants/VTPR.jl")
include("models/Activity/equations.jl")
# include("models/Activity/NRTL/eCPANRTL.jl")

include("models/Activity/COSMOSAC/utils.jl")
include("models/Activity/COSMOSAC/COSMOSAC02.jl")
include("models/Activity/COSMOSAC/COSMOSAC10.jl")
include("models/Activity/COSMOSAC/COSMOSACdsp.jl")

#Cubic variants
include("models/cubic/alphas/alphas.jl")
include("models/cubic/mixing/mixing.jl")
include("models/cubic/translation/translation.jl")
include("models/cubic/vdW/variants/Clausius.jl")
include("models/cubic/vdW/variants/Berthelot.jl")
include("models/cubic/RK/variants/SRK.jl")
include("models/cubic/RK/variants/PSRK.jl")
include("models/cubic/RK/variants/tcRK.jl")
include("models/cubic/PR/variants/PR78.jl")
include("models/cubic/PR/variants/VTPR.jl")
include("models/cubic/PR/variants/UMRPR.jl")
include("models/cubic/PR/variants/QCPR.jl")
include("models/cubic/PR/variants/tcPR.jl")
include("models/cubic/PR/variants/tcPRW.jl")
include("models/cubic/PR/variants/cPR.jl")
include("models/cubic/PR/variants/EPPR78.jl")
include("models/cubic/PatelTeja/PatelTeja.jl")
include("models/cubic/PatelTeja/variants/PatelTejaValderrama.jl")

include("models/SAFT/PCSAFT/variants/GEPCSAFT.jl")

include("models/LatticeFluid/SanchezLacombe/SanchezLacombe.jl")

include("models/Virial/Virial.jl")

include("models/ECS/ECS.jl")
include("models/ECS/variants/SPUNG.jl")
include("models/PeTS/PeTS.jl")

include("models/Electrolytes/equations.jl")
include("models/Electrolytes/base.jl")
include("models/Electrolytes/RSP/ConstRSP.jl")
include("models/Electrolytes/RSP/ZuoFurst.jl")
include("models/Electrolytes/RSP/Schreckenberg.jl")
include("models/Electrolytes/RSP/LinMixRSP.jl")
include("models/Electrolytes/Ion/Born.jl")
include("models/Electrolytes/Ion/DH.jl")
include("models/Electrolytes/Ion/DHBorn.jl")
include("models/Electrolytes/Ion/MSA.jl")
include("models/Electrolytes/Ion/MSABorn.jl")
include("models/Electrolytes/Ion/GCMSABorn.jl")
include("models/SAFT/SAFTgammaMie/variants/SAFTgammaEMie.jl")
include("models/SAFT/SAFTVRMie/variants/SAFTVREMie.jl")
include("models/SAFT/SAFTVRMie/variants/eSAFTVRMie.jl")
include("models/SAFT/PCSAFT/variants/ePCSAFT.jl")
# include("models/Electrolytes/ElectrolyteSAFT/eCPA.jl")
include("methods/property_solvers/electrolytes/electrolytes.jl")
include("methods/property_solvers/multicomponent/tp_flash/electrolyte_flash.jl")
include("models/AnalyticalSLV/AnalyticalSLV.jl")

# Export reactive systems
include("models/Reactive/reactive.jl")
include("models/Reactive/reactive_aq.jl")
include("methods/property_solvers/reactive/reactive.jl")
include("methods/property_solvers/reactive/reactive_aq.jl")

#Unitful support, transition from dependency to ext
if !isdefined(Base,:get_extension)
    include("../ext/ClapeyronUnitfulExt.jl")
end

include("estimation/estimation.jl")

#precompile workload. should be loaded at the end
#include("precompile.jl")
end # module
