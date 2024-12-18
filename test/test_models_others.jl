

@testset "Activity models" begin
    @printline
    let T = 333.15, V = 1e-3,p = 1e5,z = [0.5,0.5],z1 = Clapeyron.SA[1.0],z2 = [0.5,0.5],z3 = [0.333, 0.333,0.333];
    @testset "Wilson" begin
        system = Wilson(["methanol","benzene"])
        @test Clapeyron.activity_coefficient(system,p,T,z)[1] ≈ 1.530046633499114 rtol = 1e-6
        @test Clapeyron.activity_coefficient(system,p,T,z) ≈ Clapeyron.test_activity_coefficient(system,p,T,z)  rtol = 1e-6
        @test Clapeyron.excess_gibbs_free_energy(system,p,T,z) ≈ Clapeyron.test_excess_gibbs_free_energy(system,p,T,z)  rtol = 1e-6
    end

    @testset "NRTL" begin
        system = NRTL(["methanol","benzene"])
        @test Clapeyron.activity_coefficient(system,p,T,z)[1] ≈ 1.5309354738922405 rtol = 1e-6
    end

    @testset "aspen-NRTL" begin
        nrtl_vanilla = NRTL(["methanol","benzene"])
        system = aspenNRTL(["methanol","benzene"])
        system2 = aspenNRTL(nrtl_vanilla)
        @test Clapeyron.activity_coefficient(system,p,T,z)[1] ≈ 1.5309354738922405 rtol = 1e-6
        @test Clapeyron.activity_coefficient(system2,p,T,z)[1] ≈ 1.5309354738922405 rtol = 1e-6
    end

    @testset "UNIQUAC" begin
        system = UNIQUAC(["methanol","benzene"])
        @test Clapeyron.activity_coefficient(system,p,T,z)[1] ≈ 1.3630421218486388 rtol = 1e-6
    end

    @testset "UNIFAC" begin
        system = UNIFAC(["methanol","benzene"])
        @test Clapeyron.activity_coefficient(system,p,T,z)[1] ≈ 1.5322232657797463 rtol = 1e-6
        #when fast UNIFAC works, it should pass this test.
        # system2 = UNIFAC(["methanol","benzene"])
        # prop2 = ()
        # @test Clapeyron.activity_coefficient(system2,1e-4,423.15,[0.,1.])  ≈ [2.0807335111878937,1.0] rtol = 1e-6
    end

    @testset "ogUNIFAC" begin
        system = ogUNIFAC(["methanol","benzene"])
        @test Clapeyron.activity_coefficient(system,p,T,z)[1] ≈ 1.5133696314734384 rtol = 1e-6
    end

    @testset "UNIFAC-FV" begin
        system = UNIFACFV(["benzene","PS(1960)"])
        @test Clapeyron.activity_coefficient(system,p,T,z)[1] ≈ 0.2813003396669342 rtol = 1e-6
    end

    @testset "UNIFAC-FV-poly" begin
        system = system = UNIFACFVPoly(["PMMA(6350)","PS(1390)"])
        @test Clapeyron.activity_coefficient(system,p,T,z)[1] ≈ 2.7045808205365796 rtol = 1e-6
    end

    @testset "COSMOSAC02" begin
        system = COSMOSAC02(["water","ethanol"])
        @test Clapeyron.activity_coefficient(system,p,T,z)[1] ≈ 1.3871817962565904 rtol = 1e-6
        @test Clapeyron.excess_gibbs_free_energy(system,p,T,z) ≈ 610.5706657776052 rtol = 1e-6
    end

    @testset "COSMOSAC10" begin
        system = COSMOSAC10(["water","ethanol"])
        @test Clapeyron.activity_coefficient(system,p,T,z)[1] ≈ 1.4015660588643404 rtol = 1e-6
    end

    @testset "COSMOSACdsp" begin
        system = COSMOSACdsp(["water","ethanol"])
        @test Clapeyron.activity_coefficient(system,p,T,z)[1] ≈ 1.4398951117248127 rtol = 1e-6
    end
    end
    @printline
end

@testset "Ideal models" begin
    @printline
    let T = 298.15, V = 1e-4,p = 1e5,z = Clapeyron.SA[1.0],z1 = Clapeyron.SA[1.0],z2 = [0.5,0.5],z3 = [0.333, 0.333,0.333];
    @testset "Joback" begin
        system = JobackIdeal(["hexane"])
        @test Clapeyron.VT_isobaric_heat_capacity(system,V,298.15) ≈ 143.22076150138616 rtol = 1e-6
        @test Clapeyron.crit_pure(system)[1] ≈ 500.2728274871347 rtol = 1e-6
        @test Clapeyron.a_ideal(system,V,T,z) ≈ 9.210841420941021 rtol = 1e-6
        @test Clapeyron.ideal_consistency(system,V,T,z) ≈ 0.0 atol = 1e-14

        s0 = JobackIdeal("acetone")
        @test Clapeyron.JobackGC.T_c(s0)[1] ≈ 500.5590 rtol = 1e-6
        @test Clapeyron.JobackGC.P_c(s0)[1] ≈ 48.025e5 rtol = 1e-6
        @test Clapeyron.JobackGC.V_c(s0)[1] ≈ 209.5e-6 rtol = 1e-6
        @test Clapeyron.JobackGC.T_b(s0)[1] ≈ 322.1100 rtol = 1e-6
        @test Clapeyron.JobackGC.H_form(s0)[1] ≈ -217830.0 rtol = 1e-6
        @test Clapeyron.JobackGC.G_form(s0)[1] ≈ -154540.0 rtol = 1e-6
        @test Clapeyron.JobackGC.C_p(s0,300)[1] ≈ 75.3264 rtol = 1e-6
        @test Clapeyron.JobackGC.H_fusion(s0)[1] ≈ 5.1250e3 rtol = 1e-6
        @test Clapeyron.JobackGC.H_vap(s0)[1] ≈ 29.0180e3 rtol = 1e-6
        @test Clapeyron.JobackGC.Visc(s0,300)[1] ≈ 0.0002942 rtol = 9e-4
    end

    @testset "Reid" begin
        system = ReidIdeal(["butane"])
        @test Clapeyron.a_ideal(system,V,T,z) ≈ 9.210842104089576 rtol = 1e-6
        @test Clapeyron.ideal_consistency(system,V,T,z) ≈ 0.0 atol = 1e-14
    end

    @testset "Shomate" begin
        system = ShomateIdeal(["water"])
        coeff = system.params.coeffs[1]
        @test Clapeyron.evalcoeff(system,coeff,500) ≈ 35.21836175 rtol = 1e-6
        @test Clapeyron.eval∫coeff(system,coeff,500) ≈ 15979.2447 rtol = 1e-6
        @test Clapeyron.eval∫coeffT(system,coeff,500) ≈ 191.00554 rtol = 1e-6
    end

    @testset "Walker" begin
        system = WalkerIdeal(["hexane"])
        @test Clapeyron.molecular_weight(system)*1000 ≈ 86.21
        @test Clapeyron.a_ideal(system,V,T,z) ≈ 179.51502015696653 rtol = 1e-6
        @test Clapeyron.ideal_consistency(system,V,T,z) ≈ 0.0 atol = 1e-14
    end

    @testset "Monomer" begin
        system = MonomerIdeal(["hexane"])
        @test Clapeyron.a_ideal(system,V,T,z) ≈ -10.00711774776317 rtol = 1e-6
        @test Clapeyron.ideal_consistency(system,V,T,z) ≈ 0.0 atol = 1e-14
    end

    @testset "Empiric" begin
        #Empiric Ideal from JSON
        system = EmpiricIdeal(["water"])
        #
        @test Clapeyron.a_ideal(system,V,T,z) ≈ 7.932205569922042 rtol = 1e-6
        @test Clapeyron.ideal_consistency(system,V,T,z) ≈ 0.0 atol = 1e-14

        #Empiric Ideal from already existing MultiFluid model
        system = Clapeyron.idealmodel(MultiFluid(["water"]))
        @test Clapeyron.a_ideal(system,V,T,z) ≈ 7.932205569922042 rtol = 1e-6
        @test Clapeyron.ideal_consistency(system,V,T,z) ≈ 0.0 atol = 1e-14

        #Empiric Ideal from already existing single fluid model
        system = Clapeyron.idealmodel(system.pures[1])
        @test Clapeyron.a_ideal(system,V,T,z) ≈ 7.932205569922042 rtol = 1e-6
        @test Clapeyron.ideal_consistency(system,V,T,z) ≈ 0.0 atol = 1e-14
    end

    @testset "Aly-Lee" begin
        system = AlyLeeIdeal(["methane"])
        @test_broken Clapeyron.a_ideal(system,V,T,z) ≈ 9.239701647126086 rtol = 1e-6
        @test Clapeyron.ideal_consistency(system,V,T,z) ≈ 0.0 atol = 1e-14

        #we use the default GERG 2008 parameters for methane, test if the Cp is equal
        system_gerg = Clapeyron.idealmodel(GERG2008(["methane"]))
        Cp_system = Clapeyron.VT_isobaric_heat_capacity(system,V,T,z)
        Cp_gerg = Clapeyron.VT_isobaric_heat_capacity(system_gerg,V,T,z)

        @test Cp_system ≈ Cp_gerg rtol = 5e-5
    end

    @testset "Cp - LNG - Estimation" begin
        #Mw to obtain γ₀ = 0.708451
        system = CPLNGEstIdeal(["a1"],userlocations = (;Mw = [20.5200706797]))
        #test at 324.33 K, paper says Cp = 44.232, but the calculations in the paper seem off
        @test Clapeyron.VT_isobaric_heat_capacity(system,0.03,324.33) ≈ 44.231 rtol = 5e-4
    end

    @testset "PPDS" begin
        m1 = PPDSIdeal("krypton")
        @test isobaric_heat_capacity(m1,1,303.15)/Rgas(m1) ≈ 2.5
        mw2 = 32.042 #MonomerIdeal("methanol").params.Mw.values[1]
        m2 = PPDSIdeal("methanol")
        #verification point in ref 1, table A.6
        @test isobaric_heat_capacity(m2,1,303.15)/mw2 ≈ 1.3840 rtol = 1e-4
    end
end
    @printline
end

@testset "Multi-parameter models" begin
    @printline
    let T = 298.15, V = 1e-4,p = 1e5,z = Clapeyron.SA[1.0],z1 = Clapeyron.SA[1.0],z2 = [0.5,0.5],z3 = [0.333, 0.333,0.333];
    #warning, we are in the pseudo maxwell loop, those properties are nonsense, but they evaluate anyway.
    @printline
    @testset "IAPWS95" begin
        system = IAPWS95()
        system_ideal = Clapeyron.idealmodel(system)
        @test Clapeyron.a_ideal(system_ideal, V, T, z) ≈ 7.9322055699220435 rtol = 1e-6
        @test Clapeyron.a_ideal(system, V, T, z) ≈ 7.9322055699220435 rtol = 1e-6
        @test Clapeyron.a_res(system, V, T, z) ≈ -2.1152889226862166e14 rtol = 1e-6
        #because we are in this regime, numerical accuracy suffers. that is why big(V) is used instead.
        @test Clapeyron.ideal_consistency(system,big(V),T,z) ≈ 0.0 atol = 1e-14
    end

    @testset "PropaneRef" begin
        system = PropaneRef()
        @test Clapeyron.a_ideal(system, V, T, z) ≈ 0.6426994942361217 rtol = 1e-6
        @test Clapeyron.a_res(system, V, T, z) ≈ -2.436280448227229 rtol = 1e-6
        @test Clapeyron.ideal_consistency(system,V,T,z) ≈ 0.0 atol = 1e-14
    end

    @testset "GERG2008" begin
        system = GERG2008(["water"])
        @test Clapeyron.a_ideal(system, V, T, z) ≈ 4.500151936577565 rtol = 1e-6
        @test Clapeyron.a_res(system, V, T, z) ≈ -10.122119572808764 rtol = 1e-6
        z4   = [0.25,0.25,0.25,0.25]
        system = GERG2008(["water","carbon dioxide","hydrogen sulfide","argon"])
        @test Clapeyron.a_ideal(system, V, T, z4) ≈ 3.1136322215343917 rtol = 1e-6
        @test Clapeyron.ideal_consistency(system, V, T, z4) ≈ 0.0 rtol = 1e-14
        @test Clapeyron.a_res(system, V, T, z4) ≈ -1.1706377677539772 rtol = 1e-6
        @test Clapeyron.ideal_consistency(system,V,T,z4) ≈ 0.0 atol = 1e-14

        system_R = Clapeyron.GERG2008(["methane","ethane"],Rgas = 8.2)
        @test Clapeyron.Rgas(system_R) == 8.2
    end

    @testset "EOS-LNG" begin
        #LNG paper, table 16
        Tx = 150.0
        Vx = 1/(18002.169)
        zx   = [0.6,0.4]
        system = EOS_LNG(["methane","butane"])
        @test Clapeyron.eos(system,Vx,Tx,zx) ≈ -6020.0044 rtol = 5e-6
    end

    @testset "LKP" begin
        Tx = 150.0
        Vx = 1/(18002.169)
        zx   = [0.6,0.4]
        system = LKP(["methane","butane"])
        #Clapeyron.a_res(EOS_LNG(["methane","butane"]),V,T,z) ≈ -6.56838705236683
        @test Clapeyron.a_res(system,Vx,Tx,zx) ≈ -6.469596957611441 rtol = 5e-6
    end

    @testset "LJRef" begin
        system = LJRef(["methane"])
        Tx = 1.051*Clapeyron.T_scale(system)
        Vx = Clapeyron._v_scale(system)/0.673
        @test Clapeyron.a_ideal(system, Vx, Tx) ≈ 5.704213386278148 rtol = 1e-6
        @test Clapeyron.a_res(system, Vx, Tx) ≈ -2.244730279521925 rtol = 1e-6
        @test Clapeyron.ideal_consistency(system,Vx,Tx,z1) ≈ 0.0 atol = 1e-14
    end

    @testset "Xiang-Deiters" begin
        system = XiangDeiters(["water"])
        #equal to Clapeyron.a_ideal(BasicIdeal(["water"]), V, T, z)
        @test Clapeyron.a_ideal(system, V, T, z1) ≈ -0.33605470137749016 rtol = 1e-6
        @test Clapeyron.a_res(system, V, T)  ≈ -34.16747927719535 rtol = 1e-6
    end
    @printline
    end
end

@testset "SPUNG models" begin
    @printline
    let T = 298.15, V = 1e-4,p = 1e5,z = Clapeyron.SA[1.0],z1 = Clapeyron.SA[1.0],z2 = [0.5,0.5],z3 = [0.333, 0.333,0.333];
    @testset "SRK" begin
        system = SPUNG(["ethane"])
        @test Clapeyron.shape_factors(system, V, T, z)[1] ≈ 0.8246924617474896 rtol = 1e-6
    end


    @testset "PCSAFT" begin
        system = SPUNG(["ethane"],PropaneRef(),PCSAFT(["ethane"]),PCSAFT(["propane"]))
        @test Clapeyron.shape_factors(system, V, T, z)[1] ≈ 0.8090183134644525 rtol = 1e-6
    end
    end
end

@testset "lattice models" begin
    @printline
    let T = 298.15, V = 1e-4,p = 1e5,z = Clapeyron.SA[1.0],z1 = Clapeyron.SA[1.0],z2 = [0.5,0.5],z3 = [0.333, 0.333,0.333];

    @testset "single component" begin
        system = Clapeyron.SanchezLacombe(["carbon dioxide"])
        @test Clapeyron.a_res(system, V, T, z) ≈ -0.9511044462267396 rtol = 1e-6
    end

    @testset "Sanchez-Lacombe,Kij rule" begin
        system = SanchezLacombe(["carbon dioxide","benzoic acid"],mixing = SLKRule)
        @test Clapeyron.a_res(system, V, T, z2) ≈ -6.494291842858994 rtol = 1e-6
    end

    @testset "Sanchez-Lacombe K0-K1-L rule" begin
        system = SanchezLacombe(["carbon dioxide","benzoic acid"],mixing = SLk0k1lMixingRule)
        @test Clapeyron.a_res(system, V, T, z2) ≈ -5.579621796375229 rtol = 1e-6
    end
    end
end


@testset "Correlations" begin
    @testset "DIPPR101Sat" begin
        system = DIPPR101Sat(["water"])
        p0 = saturation_pressure(system,400.01)[1]
        @test p0 ≈ 245338.15099198322 rtol = 1e-6
        @test saturation_temperature(system,p0)[1] ≈ 400.01 rtol = 1e-6
    end

    @testset "LeeKeslerSat" begin
        system = LeeKeslerSat(["water"])
        p0 = saturation_pressure(system,400.01)[1]
        @test p0 ≈ 231731.79240876858 rtol = 1e-6
        @test saturation_temperature(system,p0)[1] ≈ 400.01 rtol = 1e-6
    end

    @testset "COSTALD" begin
        system = COSTALD(["water"])
        @test volume(system,1e5,300.15) ≈ 1.8553472145724288e-5 rtol = 1e-6
        system2 = COSTALD(["water","methanol"])
        @test volume(system2,1e5,300.15,[0.5,0.5]) ≈ 2.834714146558056e-5 rtol = 1e-6
        @test volume(system2,1e5,300.15,[1.,0.]) ≈ 1.8553472145724288e-5 rtol = 1e-6
    end

    @testset "RackettLiquid" begin
        system = RackettLiquid(["water"])
        @test volume(system,1e5,300.15) ≈ 1.6837207241594103e-5 rtol = 1e-6
        system2 = RackettLiquid(["water","methanol"])
        @test volume(system2,1e5,300.15,[0.5,0.5]) ≈ 3.2516352601748416e-5 rtol = 1e-6
        @test volume(system2,1e5,300.15,[1.,0.]) ≈ 1.6837207241594103e-5 rtol = 1e-6
    end


    @testset "AbbottVirial" begin
        system = AbbottVirial(["methane","ethane"])
        @test volume(system,1e5,300,[0.5,0.5]) ≈ 0.024820060368027988 rtol = 1e-6
        #a_res(PR,0.05,300,[0.5,0.5]) == -0.0023705490820905483
        @test Clapeyron.a_res(system,0.05,300,[0.5,0.5]) ≈ -0.0024543543773067693 rtol = 1e-6
    end

    @testset "TsonopoulosVirial" begin
        system = TsonopoulosVirial(["methane","ethane"])
        @test volume(system,1e5,300,[0.5,0.5]) ≈ 0.02485310667780686 rtol = 1e-6
        #a_res(PR,0.05,300,[0.5,0.5]) == -0.0023705490820905483
        @test Clapeyron.a_res(system,0.05,300,[0.5,0.5]) ≈ -0.0017990881811592349 rtol = 1e-6
    end

    @testset "EoSVirial2" begin
        cub = PR(["methane","ethane"])
        system = EoSVirial2(cub)
        #exact equality here, as cubics have an exact second virial coefficient
        @test volume(system,1e5,300,[0.5,0.5]) == Clapeyron.volume_virial(cub,1e5,300,[0.5,0.5])
        #a_res(PR,0.05,300,[0.5,0.5]) == -0.0023705490820905483
        @test Clapeyron.a_res(system,0.05,300,[0.5,0.5]) ≈ -0.0023728381262076137 rtol = 1e-6
    end

    @testset "SolidHfus" begin
        model = SolidHfus(["water"])
        @test chemical_potential(model,1e5,298.15,[1.])[1] ≈ 549.1488193300384 rtol = 1e-6
    end
end
