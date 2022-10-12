function output = ...
    Tnnp06_model(t,X,flag_ode,dt,PGkr,PGks,PGCaL,PGNa,celltype,ex)
    
    %calculate the stimulus current, Istim
    Istim = ex;
    
    %Constants
    R = 8314.472;
    F = 96485.3415;
    T = 310.0;
    RTONF = (R*T) / F;

    Volt = X(1);
    Cai = X(2);
    CaSR = X(3);
    CaSS = X(4);
    Nai = X(5);
    Ki = X(6);
    INa_m = X(7);
    INa_h = X(8);
    INa_j = X(9);
    IKr_xr1 = X(10);
    IKr_xr2 = X(11);
    IKs_xs = X(12);
    Ito_r = X(13);
    Ito_s = X(14);
    ICaL_d = X(15);
    ICaL_f = X(16);
    ICaL_f2 = X(17);
    ICaL_fCaSS = X(18);
    RR = X(19);
    %OO = X(20);

    %External concentrations
    Ko = 5.4;
    Cao = 2.0;
    Nao = 140.0;

    %Intracellular volumes
    Vc = 0.016404;
    Vsr = 0.001094;
    Vss = 0.00005468;

    %Calcium buffering dynamics
    Bufc = 0.2;
    Kbufc = 0.001;
    Bufsr = 10.;
    Kbufsr = 0.3;
    Bufss = 0.4;
    Kbufss = 0.00025;

    %Intracellular calcium flux dynamics
    Vmaxup = 0.006375;
    Kup = 0.00025;
    Vrel = 0.102;%40.8;
    k1_ = 0.15;
    k2_ = 0.045;
    k3 = 0.060;
    k4 = 0.005;%0.000015;
    EC = 1.5;
    maxsr = 2.5;
    minsr = 1.;
    Vleak = 0.00036;
    Vxfer = 0.0038;


    %Cellular capacitance         
    CAPACITANCE = 0.185;

    %Parameters for currents
    %Parameters for IKr
    Gkr = 0.153;
    %Parameters for Iks
    pKNa = 0.03;
    if strcmp('Epi',celltype)
        Gks = 0.392;
    elseif strcmp('Endo', celltype)
        Gks = 0.392;
    else
        Gks = 0.098;
    end
    % //Parameters for Ik1
    GK1 = 5.405;
    % //Parameters for Ito
    if strcmp('Epi',celltype)
        Gto = 0.294;
    elseif strcmp('Endo', celltype)
        Gto = 0.073;
    else
        Gto = 0.294;
    end
    % //Parameters for INa
    GNa = 14.838*PGNa;
    % //Parameters for IbNa
    GbNa = 0.00029;
    % //Parameters for INaK
    KmK = 1.0;
    KmNa = 40.0;
    knak = 2.724;
    % //Parameters for ICaL
    GCaL = 0.00003980;
    % //Parameters for IbCa
    GbCa = 0.000592;
    % //Parameters for INaCa
    knaca = 1000;
    KmNai = 87.5;
    KmCa = 1.38;
    ksat = 0.1;
    n = 0.35;
    % //Parameters for IpCa
    GpCa = 0.1238;
    KpCa = 0.0005;
    % //Parameters for IpK;
    GpK = 0.0146;

    %IKr
    IKr_xr1_inf = 1. / (1. + exp((-26. - Volt) / 7.));
    IKr_xr1_a = 450. / (1. + exp((-45. - Volt) / 10.));
    IKr_xr1_b = 6. / (1. + exp((Volt - (-30.)) / 11.5));
    IKr_xr1_t = IKr_xr1_a * IKr_xr1_b;
    IKr_xr2_inf = 1. / (1. + exp((Volt - (-88.)) / 24.));
    IKr_xr2_a = 3. / (1. + exp((-60. - Volt) / 20.));
    IKr_xr2_b = 1.12 / (1. + exp((Volt - 60.) / 20.));
    IKr_xr2_t = IKr_xr2_a*IKr_xr2_b;
    % dIKr_xr1=(IKr_xr1_inf-IKr_xr1)/IKr_xr1_t;
    % dIKr_xr2=(IKr_xr2_inf-IKr_xr2)/IKr_xr2_t;
    
    Ek = RTONF*(log((Ko / Ki)));
    sxr1 = IKr_xr1;
    sxr2 = IKr_xr2;
    IKr = Gkr*sqrt(Ko / 5.4)*sxr1*sxr2*(Volt - Ek);
    IKr = IKr * PGkr;

    %IKs
    IKs_xs_inf = 1. / (1. + exp((-5. - Volt) / 14.));
	IKs_xs_a = (1400. / (sqrt(1. + exp((5. - Volt) / 6))));
	IKs_xs_b = (1. / (1. + exp((Volt - 35.) / 15.)));
	IKs_xs_t = IKs_xs_a*IKs_xs_b + 80;
    % dIKs_xs = (IKs_xs_inf-IKs_xs)/IKs_xs_t;
    
    Eks = RTONF*(log((Ko + pKNa*Nao) / (Ki + pKNa*Nai)));
    IKs = Gks*IKs_xs*IKs_xs*(Volt - Eks);
    IKs = IKs * PGks;
    
    %IK1
    IK1_a = 0.1 / (1. + exp(0.06*(Volt - Ek - 200)));
	IK1_b = (3.*exp(0.0002*(Volt - Ek + 100))+exp(0.1*(Volt - Ek - 10))) /...
        (1. + exp(-0.5*(Volt - Ek)));
    rec_iK1 = IK1_a / (IK1_a + IK1_b);
    IK1 = GK1*rec_iK1*(Volt - Ek);
    
    %Ito
    if strcmp('Epi',celltype)
        Ito_r_inf = 1. / (1. + exp((20 - Volt) / 6.));
        Ito_s_inf = 1. / (1. + exp((Volt + 20) / 5.));
        Ito_r_t = 9.5*exp(-(Volt + 40.)*(Volt + 40.) / 1800.) + 0.8;
        Ito_s_t = 85.*exp(-(Volt + 45.)*(Volt + 45.) / 320.) + 5. /...
            (1. + exp((Volt - 20.) / 5.)) + 3.;
    elseif strcmp('Endo', celltype)
        Ito_r_inf = 1. / (1. + exp((20 - Volt) / 6.));
        Ito_s_inf = 1. / (1. + exp((Volt + 28) / 5.));
        Ito_r_t = 9.5*exp(-(Volt + 40.)*(Volt + 40.) / 1800.) + 0.8;
        Ito_s_t = 1000.*exp(-(Volt + 67)*(Volt + 67) / 1000.) + 8.;
    else
        Ito_r_inf = 1. / (1. + exp((20 - Volt) / 6.));
        Ito_s_inf = 1. / (1. + exp((Volt + 20) / 5.));
        Ito_r_t = 9.5*exp(-(Volt + 40.)*(Volt + 40.) / 1800.) + 0.8;
        Ito_s_t = 85.*exp(-(Volt + 45.)*(Volt + 45.) / 320.) + 5. / (1. + exp((Volt - 20.) / 5.)) + 3.;
    end
    % dIto_r = (Ito_r_inf-Ito_r)/Ito_r_t;
    % dIto_s = (Ito_s_inf-Ito_s)/Ito_s_t;
    Ito = Gto*Ito_r*Ito_s*(Volt - Ek);
    
    %INa
    Ena = RTONF*(log((Nao / Nai)));
    INa_m_a = 1. / (1. + exp((-60. - Volt) / 5.));
	INa_m_b = 0.1 / (1. + exp((Volt + 35.) / 5.)) + 0.10 / (1. + exp((Volt - 50.) / 200.));
	INa_m_t = INa_m_a*INa_m_b;
	INa_m_inf = 1. / ((1. + exp((-56.86 - Volt) / 9.03))*(1. + exp((-56.86 - Volt) / 9.03)));
    if Volt >= -40.
		INa_h_a = 0.;
		INa_h_b = (0.77 / (0.13*(1. + exp(-(Volt + 10.66) / 11.1))));
		INa_h_t = 1.0 / (INa_h_a + INa_h_b);
    else
		INa_h_a = (0.057*exp(-(Volt + 80.) / 6.8));
		INa_h_b = (2.7*exp(0.079*Volt) + (3.1e5)*exp(0.3485*Volt));
		INa_h_t = 1.0 / (INa_h_a + INa_h_b);
    end
	INa_h_inf = 1. / ((1. + exp((Volt + 71.55) / 7.43))*(1. + exp((Volt + 71.55) / 7.43)));
    if Volt >= -40.
		INa_j_a = 0.;
		INa_j_b = (0.6*exp((0.057)*Volt) / (1. + exp(-0.1*(Volt + 32.))));
		INa_j_t = 1.0 / (INa_j_a + INa_j_b);
    else
		INa_j_a = (((-2.5428e4)*exp(0.2444*Volt)-(6.948e-6)*exp(-0.04391*Volt))...
            *(Volt + 37.78)/(1. + exp(0.311*(Volt + 79.23))));
		INa_j_b = (0.02424*exp(-0.01052*Volt) / (1. + exp(-0.1378*(Volt + 40.14))));
		INa_j_t = 1.0 / (INa_j_a + INa_j_b);
    end
	INa_j_inf = INa_h_inf;
    % dINa_m = (INa_m_inf-INa_m)/INa_m_t;
    % dINa_h = (INa_h_inf-INa_h)/INa_h_t;
    % dINa_j = (INa_j_inf-INa_j)/INa_j_t;
    INa = GNa*INa_m^3*INa_h*INa_j*(Volt - Ena);
    
    %IbNa
    IbNa = GbNa*(Volt - Ena);
    
    %INaK
    rec_iNaK = (1. / (1. + 0.1245*exp(-0.1*Volt*F / (R*T)) + 0.0353*...
        exp(-Volt*F / (R*T))));
    INaK = knak*(Ko / (Ko + KmK))*(Nai / (Nai + KmNa))*rec_iNaK;
    
        
    %ICaL
    ICaL_d_inf = 1. / (1. + exp((-8 - Volt) / 7.5));
	ICaL_d_a = 1.4 / (1. + exp((-35 - Volt) / 13)) + 0.25;
	ICaL_d_b = 1.4 / (1. + exp((Volt + 5) / 5));
	ICaL_d_c = 1. / (1. + exp((50 - Volt) / 20));
	ICaL_d_t = ICaL_d_a*ICaL_d_b + ICaL_d_c;
    
	ICaL_f_inf = 1. / (1. + exp((Volt + 20) / 7));
	ICaL_f_a = 1102.5*exp(-(Volt + 27)*(Volt + 27) / 225);
	ICaL_f_b = 200. / (1 + exp((13 - Volt) / 10.));
	ICaL_f_c = (180. / (1 + exp((Volt + 30) / 10))) + 20;
	ICaL_f_t = ICaL_f_a + ICaL_f_b + ICaL_f_c;
    
	ICaL_f2_inf = 0.67 / (1. + exp((Volt + 35) / 7)) + 0.33;
	ICaL_f2_a = 600 * exp(-(Volt + 25)*(Volt + 25) / 170);
	ICaL_f2_b = 31 / (1. + exp((25 - Volt) / 10));
	ICaL_f2_c = 16 / (1. + exp((Volt + 30) / 10));
	ICaL_f2_t = ICaL_f2_a + ICaL_f2_b + ICaL_f2_c;
    
	ICaL_fCaSS_inf = 0.6 / (1 + (CaSS / 0.05)*(CaSS / 0.05)) + 0.4;
	ICaL_fCaSS_t = 80. / (1 + (CaSS / 0.05)*(CaSS / 0.05)) + 2.;
    
    % dCaL_d = (ICaL_d_inf-ICaL_d)/ICaL_d_t;
    % dCaL_f = (ICaL_f_inf-ICaL_f)/ICaL_f_t;
    % dCaL_f2 = (ICaL_f2_inf-ICaL_f2)/ICaL_f2_t;
    % dCaL_fCaSS = (ICaL_fCaSS_inf-ICaL_fCaSS)/ICaL_fCaSS_t;
    ICaL = GCaL*ICaL_d*ICaL_f*ICaL_f2*ICaL_fCaSS * 4 * (Volt - 15)*...
        (F^2 / (R*T))*(0.25*exp(2 * (Volt - 15)*F /...
        (R*T))*CaSS - Cao) / (exp(2 * (Volt - 15)*F / (R*T)) - 1.);
    ICaL = ICaL * PGCaL;
    
    %IbCa
    ECa = 0.5*RTONF*(log((Cao / Cai)));
    IbCa = GbCa*(Volt - ECa);
    
    %INaCa
    INaCa = knaca*(1. / (KmNai^3 + Nao^3))*(1. / (KmCa + Cao))*...
        (1. / (1 + ksat*exp((n - 1)*Volt*F / (R*T))))*(exp(n*Volt*F / (R*T))...
        *Nai^3*Cao -	exp((n - 1)*Volt*F / (R*T))*Nao^3*Cai*2.5);
    
    %IpCa
    IpCa = GpCa*Cai / (KpCa + Cai);
    
    %IpK
    rec_ipK = 1. / (1. + exp((25 - Volt) / 5.98));
    IpK = GpK*rec_ipK*(Volt - Ek);
    
    %update the membrane voltage
    dv = -(IKr+IKs+IK1+Ito+INa+IbNa+INaK+ICaL+IbCa+INaCa+IpCa+IpK+Istim);
    
    %Ca transient and intracellular concentrations
    inverseVcF2 = 1 / (2 * Vc*F);
	inverseVcF = 1. / (Vc*F);
	inversevssF2 = 1 / (2 * Vss*F);
    
    kCaSR = maxsr - ((maxsr - minsr) / (1 + (EC / CaSR)*(EC / CaSR)));
	k1 = k1_ / kCaSR;
	k2 = k2_ * kCaSR;
	dRR = k4 * (1 - RR) - k2*CaSS*RR;
	OO = k1*CaSS*CaSS*RR / (k3 + k1*CaSS*CaSS);


	Irel = Vrel*OO*(CaSR - CaSS);
	Ileak = Vleak*(CaSR - Cai);
	Iup = Vmaxup / (1. + ((Kup^2) / (Cai^2)));
	Ixfer = Vxfer*(CaSS - Cai);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%    
    CaCSQN = Bufsr*CaSR / (CaSR + Kbufsr);
	dCaSR = (Iup - Irel - Ileak)*dt;
	bjsr = Bufsr - CaCSQN - dCaSR - CaSR + Kbufsr;
	cjsr = Kbufsr*(CaCSQN + dCaSR + CaSR);
	CaSR = (sqrt(bjsr*bjsr + 4 * cjsr) - bjsr) / 2;

	CaSSBuf = Bufss*CaSS / (CaSS + Kbufss);
	dCaSS = (-Ixfer*(Vc / Vss) + Irel*(Vsr / Vss) +...
        (-ICaL*inversevssF2*CAPACITANCE))*dt;
	bcss = Bufss - CaSSBuf - dCaSS - CaSS + Kbufss;
	ccss = Kbufss*(CaSSBuf + dCaSS + CaSS);
	CaSS = (sqrt(bcss*bcss + 4 * ccss) - bcss) / 2;

	CaBuf = Bufc*Cai / (Cai + Kbufc);
	dCai = ((-(IbCa + IpCa - 2 * INaCa)*inverseVcF2*CAPACITANCE) -...
        (Iup - Ileak)*(Vsr / Vc) + Ixfer)*dt;
	bc = Bufc - CaBuf - dCai - Cai + Kbufc;
	cc = Kbufc*(CaBuf + dCai + Cai);
	Cai = (sqrt(bc*bc + 4 * cc) - bc) / 2;

	dNai = -(INa + IbNa + 3 * INaK + 3 * INaCa)*inverseVcF*CAPACITANCE;
	dKi = -(Istim + IK1 + Ito + IKr + IKs - 2 * INaK + IpK)*inverseVcF*CAPACITANCE;
	
	RR = dRR*dt + RR;
    Nai = dNai*dt + Nai;
    Ki = dKi*dt + Ki;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    INa_m = INa_m_inf-(INa_m_inf-INa_m)*exp(-dt/INa_m_t);
    INa_h = INa_h_inf-(INa_h_inf-INa_h)*exp(-dt/INa_h_t);
    INa_j = INa_j_inf-(INa_j_inf-INa_j)*exp(-dt/INa_j_t);
    IKr_xr1 = IKr_xr1_inf-(IKr_xr1_inf-IKr_xr1)*exp(-dt/IKr_xr1_t);
    IKr_xr2 = IKr_xr2_inf-(IKr_xr2_inf-IKr_xr2)*exp(-dt/IKr_xr2_t);
    IKs_xs = IKs_xs_inf-(IKs_xs_inf-IKs_xs)*exp(-dt/IKs_xs_t);
    Ito_s= Ito_s_inf-(Ito_s_inf-Ito_s)*exp(-dt/Ito_s_t);
    Ito_r= Ito_r_inf-(Ito_r_inf-Ito_r)*exp(-dt/Ito_r_t);
    ICaL_d = ICaL_d_inf-(ICaL_d_inf-ICaL_d)*exp(-dt/ICaL_d_t);
    ICaL_f =ICaL_f_inf-(ICaL_f_inf-ICaL_f)*exp(-dt/ICaL_f_t);
    ICaL_f2 =ICaL_f2_inf-(ICaL_f2_inf-ICaL_f2)*exp(-dt/ICaL_f2_t);
    ICaL_fCaSS =ICaL_fCaSS_inf-(ICaL_fCaSS_inf-ICaL_fCaSS)*exp(-dt/ICaL_fCaSS_t);
	
    Volt = dv*dt + Volt;
	
    %return values
    if flag_ode == 1
		X(1) = Volt;X(2) = Cai;X(3) = CaSR;X(4) = CaSS;X(5) = Nai;
		X(6) = Ki;X(7) = INa_m;X(8) = INa_h;X(9) = INa_j;X(10) = IKr_xr1;
		X(11) = IKr_xr2;X(12) = IKs_xs;X(13) = Ito_r;X(14) = Ito_s;
		X(15) = ICaL_d;X(16) = ICaL_f;X(17) = ICaL_f2;X(18) = ICaL_fCaSS;
		X(19) = RR;
        output = X;
    else 
        output = [IKr IKs IK1 Ito INa IbNa INaK ICaL IbCa INaCa IpCa IpK Istim];
    end
end