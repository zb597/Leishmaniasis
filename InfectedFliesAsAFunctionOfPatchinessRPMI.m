%Mice
MF = [1918.327	5368.238	14975.54	3805.637	4451.461	17625.02	27372.53	17554.89	250.9739	11183.28	7324	7664.599	1906.456	5972.289	5060.821	11847.68	7574.186	717.4968	1669.183	10734.92	43123.45	2413.856	731.7327	42438.04];
MJ = [15997.09	3017.408	5034.854	9239.502	2739.726	513.6808	1732.656	1191.046	1036.687	3348.615	1548.505	860.771	7333.507	1555.587	2160.678	3284.401	3297.141	922.8627	4877.168	6981.653	555.0819	763.1017	27997.12	45156.43];
M34 = [79325.4	8075.2	12407.7	1621.6	12520.3	14747.8	6971.3	4854.6	2778.3	1358.6	14703.1	5527.3	3123.6	3629.4	5529.5	4353.3	6067.5	2363.5	4163.9	2053.3	5635.2	3689.9	8057.8	12618.4];
M10 = [4243310	5552395	4087400	4627846	5546594	17248300	17070900	9188658	3691940	25452820	34302760	23473170	7529109	9002134	4844421	3201581	3572419	5662281	9909828	4986509	1713397	193055.3	593888.8	833192.3];
M35 = [32331.17	11941.36	30553.34	23279.64	148774.2	335843	42148.09	38192.13	12640.68	31603.28	103105.3	68400.34	76280.03	174505.1	40626.54	108187.5	23332.35	43019.24	114974.5	98845.32	6164.186	13454.91	26539.47	11808.43];
M4 = [6816.114	9865.926	8074.566	2144.396	3790.682	4670.691	287609.7	17519.6	675.1464	2809.901	4721.435	17965.03	10342.02	15399.22	1696.771	16272.53	3538.723	1732.837	4501.902	630.3101	6376.592	79432.82	1688.43	4939.643];
M7 = [8176.658	4045.592	11200.34	27846.53	7319.354	3848.755	2961.245	10853.75	13593.43	82828.38	5796.723	14738.95	44362.54	2780.33	11222.72	18019.24	31724.3	4964.867	7930.374	8426.884	40682.51	12424.21	6102.104	2921.042];
M9 = [94325.42	62468.18	71573.48	59310.95	74108.88	103655.9	51661.73	193353.3	83582.15	84971.95	27241.44	274852.2	253882.5	50821.75	153677.9	296190.4	44208.68	353988.8	340442.5	128945	28428.64	28618.89	6126.26	7119.881];
M11 = [43841.88	29074.17	13066.45	81767.38	119878.5	33898.7	165113.4	100725.1	73412.82	147603.8	58429.43	247062.3	122996.5	214556.3	62988.18	116733.1	114132	179094.4	68986.91	118917.3	35624.51	3932.968	12931.35	70823.13];

MeanA = mean(MF); MeanB = mean(MJ); MeanC = mean(M34); MeanD = mean(M10); MeanE = mean(M35); MeanF = mean(M4); MeanG = mean(M7); MeanH = mean(M9); MeanI = mean(M11); 
A = round(reshape(MF, [numel(MF),1])); B = round(reshape(MJ, [numel(MJ),1])); C = round(reshape(M34, [numel(M34),1])); D = round(reshape(M10, [numel(M10),1])); E = round(reshape(M35, [numel(M35),1])); F = round(reshape(M4, [numel(M4),1])); G = round(reshape(M7, [numel(M7),1])); H = round(reshape(M9, [numel(M9),1])); I = round(reshape(M11, [numel(M11),1]));

%Fit Negative Binomial Distribution.
MaNBA = fitdist(A, 'Negative Binomial'); MaNBB = fitdist(B, 'Negative Binomial'); MaNBC = fitdist(C, 'Negative Binomial'); MaNBD = fitdist(D, 'Negative Binomial'); MaNBE = fitdist(E, 'Negative Binomial'); MaNBF = fitdist(F, 'Negative Binomial'); MaNBG = fitdist(G, 'Negative Binomial'); MaNBH = fitdist(H, 'Negative Binomial'); MaNBI = fitdist(I, 'Negative Binomial'); 

%Volume of punch biopsy in m^3
PunchVol = 5.03 * 10^-8 ;

%Total BP volume of 1 uL (cuboidal dimensions = 1.4mm x 1.4mm x 500um)
BPVol = 1 * 10^-9;

%Feed volume of fly
FeedVol = 0.7*10^-9;

for ii = 1:1000;
    
    %%%% MICRO-PATCHINESS %%%%
    
    %Number of parasites at macro site
    %Use if macro is patchy
    MaPPPA = nbinrnd(MaNBA.r, MaNBA.p); MaPPPB = nbinrnd(MaNBB.r, MaNBB.p); MaPPPC = nbinrnd(MaNBC.r, MaNBC.p); MaPPPD = nbinrnd(MaNBD.r, MaNBD.p); MaPPPE = nbinrnd(MaNBE.r, MaNBE.p); MaPPPF = nbinrnd(MaNBF.r, MaNBF.p); MaPPPG = nbinrnd(MaNBG.r, MaNBG.p); MaPPPH = nbinrnd(MaNBH.r, MaNBH.p); MaPPPI = nbinrnd(MaNBI.r, MaNBI.p); 
        
    %Parasite density in punch biopsy
    MaPDPPA = MaPPPA/PunchVol; MaPDPPB = MaPPPB/PunchVol; MaPDPPC = MaPPPC/PunchVol; MaPDPPD = MaPPPD/PunchVol; MaPDPPE = MaPPPE/PunchVol; MaPDPPF = MaPPPF/PunchVol; MaPDPPG = MaPPPG/PunchVol; MaPDPPH = MaPPPH/PunchVol; MaPDPPI = MaPPPI/PunchVol;
    
    %Expected number of parasites
    ExpParaPPA = MaPDPPA * BPVol; ExpParaPPB = MaPDPPB * BPVol; ExpParaPPC = MaPDPPC * BPVol; ExpParaPPD = MaPDPPD * BPVol; ExpParaPPE = MaPDPPE * BPVol; ExpParaPPF = MaPDPPF * BPVol; ExpParaPPG = MaPDPPG * BPVol; ExpParaPPH = MaPDPPH * BPVol; ExpParaPPI = MaPDPPI * BPVol; 
     
    %Number of parasites at a micro sites is a random number from the
    %exponential distribution with the mean from the macro data.
    MiPPPA = exprnd(ExpParaPPA); MiPPPB = exprnd(ExpParaPPB); MiPPPC = exprnd(ExpParaPPC); MiPPPD = exprnd(ExpParaPPD); MiPPPE = exprnd(ExpParaPPE); MiPPPF = exprnd(ExpParaPPF); MiPPPG = exprnd(ExpParaPPG); MiPPPH = exprnd(ExpParaPPH); MiPPPI = exprnd(ExpParaPPI);
    
    %Concentration of parasites in BP
    ConcPPA = MiPPPA/BPVol; ConcPPB = MiPPPB/BPVol; ConcPPC = MiPPPC/BPVol; ConcPPD = MiPPPD/BPVol; ConcPPE = MiPPPE/BPVol; ConcPPF = MiPPPF/BPVol; ConcPPG = MiPPPG/BPVol; ConcPPH = MiPPPH/BPVol; ConcPPI = MiPPPI/BPVol;   
    
    %Number of parasites picked up by the fly on feeding on the full 0.7 uL
    %feed volume %BP Scale
    FlyParaPPA(ii) = ConcPPA * FeedVol; FlyParaPPB(ii) = ConcPPB * FeedVol; FlyParaPPC(ii) = ConcPPC * FeedVol; FlyParaPPD(ii) = ConcPPD * FeedVol; FlyParaPPE(ii) = ConcPPE * FeedVol; FlyParaPPF(ii) = ConcPPF * FeedVol; FlyParaPPG(ii) = ConcPPG * FeedVol; FlyParaPPH(ii) = ConcPPH * FeedVol; FlyParaPPI(ii) = ConcPPI * FeedVol; 
      
    %%%INFECTION STATUS
    
    %Number of parasites 
    PoTPPA(ii) = round(FlyParaPPA(ii)); PoTPPB(ii) = round(FlyParaPPB(ii)); PoTPPC(ii) = round(FlyParaPPC(ii)); PoTPPD(ii) = round(FlyParaPPD(ii)); PoTPPE(ii) = round(FlyParaPPE(ii)); PoTPPF(ii) = round(FlyParaPPF(ii)); PoTPPG(ii) = round(FlyParaPPG(ii)); PoTPPH(ii) = round(FlyParaPPH(ii)); PoTPPI(ii) = round(FlyParaPPI(ii));
    
    if PoTPPA(ii) > 500 %This is the infection threshold
        
        pNI = 0.75;
        pInfPPA = 1 - (pNI^PoTPPA(ii));
        InfProbPPA(ii) = pInfPPA;
        
    else
        InfProbPPA(ii) = 0;
        
    end
    
    rnPPA = rand;
    
    if rnPPA < InfProbPPA(ii);
        InfYesNoPPA(ii) = 1;
    else
        InfYesNoPPA(ii) = 0;
    end
    
    if PoTPPB(ii) > 500 %This is the infection threshold
        
        pNI = 0.75;
        pInfPPB = 1 - (pNI^PoTPPB(ii));
        InfProbPPB(ii) = pInfPPB;
        
    else
        InfProbPPB(ii) = 0;
        
    end
    
    rnPPB = rand;
    
    if rnPPB < InfProbPPB(ii);
        InfYesNoPPB(ii) = 1;
    else
        InfYesNoPPB(ii) = 0;
    end
    
    if PoTPPC(ii) > 500 %This is the infection threshold
        
        pNI = 0.75;
        pInfPPC = 1 - (pNI^PoTPPC(ii));
        InfProbPPC(ii) = pInfPPC;
        
    else
        InfProbPPC(ii) = 0;
        
    end
    
    rnPPC = rand;
    
    if rnPPC < InfProbPPC(ii);
        InfYesNoPPC(ii) = 1;
    else
        InfYesNoPPC(ii) = 0;
    end
    
    if PoTPPD(ii) > 500 %This is the infection threshold
        
        pNI = 0.75;
        pInfPPD = 1 - (pNI^PoTPPD(ii));
        InfProbPPD(ii) = pInfPPD;
        
    else
        InfProbPPD(ii) = 0;
        
    end
    
    rnPPD = rand;
    
    if rnPPD < InfProbPPD(ii);
        InfYesNoPPD(ii) = 1;
    else
        InfYesNoPPD(ii) = 0;
    end
    
    if PoTPPE(ii) > 500 %This is the infection threshold
        
        pNI = 0.75;
        pInfPPE = 1 - (pNI^PoTPPE(ii));
        InfProbPPE(ii) = pInfPPE;
        
    else
        InfProbPPE(ii) = 0;
        
    end
    
    rnPPE = rand;
    
    if rnPPE < InfProbPPE(ii);
        InfYesNoPPE(ii) = 1;
    else
        InfYesNoPPE(ii) = 0;
    end
    
    if PoTPPF(ii) > 500 %This is the infection threshold
        
        pNI = 0.75;
        pInfPPF = 1 - (pNI^PoTPPF(ii));
        InfProbPPF(ii) = pInfPPF;
        
    else
        InfProbPPF(ii) = 0;
        
    end
    
    rnPPF = rand;
    
    if rnPPF < InfProbPPF(ii);
        InfYesNoPPF(ii) = 1;
    else
        InfYesNoPPF(ii) = 0;
    end
    
    if PoTPPG(ii) > 500 %This is the infection threshold
        
        pNI = 0.75;
        pInfPPG = 1 - (pNI^PoTPPG(ii));
        InfProbPPG(ii) = pInfPPG;
        
    else
        InfProbPPG(ii) = 0;
        
    end
    
    rnPPG = rand;
    
    if rnPPG < InfProbPPG(ii);
        InfYesNoPPG(ii) = 1;
    else
        InfYesNoPPG(ii) = 0;
    end
    
    if PoTPPH(ii) > 500 %This is the infection threshold
        
        pNI = 0.75;
        pInfPPH = 1 - (pNI^PoTPPH(ii));
        InfProbPPH(ii) = pInfPPH;
        
    else
        InfProbPPH(ii) = 0;
        
    end
    
    rnPPH = rand;
    
    if rnPPH < InfProbPPH(ii);
        InfYesNoPPH(ii) = 1;
    else
        InfYesNoPPH(ii) = 0;
    end
    
    if PoTPPI(ii) > 500 %This is the infection threshold
        
        pNI = 0.75;
        pInfPPI = 1 - (pNI^PoTPPI(ii));
        InfProbPPI(ii) = pInfPPI;
        
    else
        InfProbPPI(ii) = 0;
        
    end
    
    rnPPI = rand;
    
    if rnPPI < InfProbPPI(ii);
        InfYesNoPPI(ii) = 1;
    else
        InfYesNoPPI(ii) = 0;
    end
    
end

YesPPA = nnz(InfYesNoPPA); YesPPB = nnz(InfYesNoPPB); YesPPC = nnz(InfYesNoPPC); YesPPD = nnz(InfYesNoPPE); YesPPE = nnz(InfYesNoPPE); YesPPF = nnz(InfYesNoPPF); YesPPG = nnz(InfYesNoPPG); YesPPH = nnz(InfYesNoPPH); YesPPI = nnz(InfYesNoPPI);

for ii = 1:1000;
    
    %%%% MICRO-PATCHINESS %%%%
    
    %Number of parasites at macro site   
    %Use if macro is not patchy
    MaPNPPA = MeanA; MaPNPPB = MeanB; MaPNPPC = MeanC; MaPNPPD = MeanD; MaPNPPE = MeanE; MaPNPPF = MeanF; MaPNPPG = MeanG; MaPNPPH = MeanH; MaPNPPI = MeanI; 
    
    %Parasite density in punch biopsy
    MaPDNPPA = MaPNPPA/PunchVol; MaPDNPPB = MaPNPPB/PunchVol; MaPDNPPC = MaPNPPC/PunchVol; MaPDNPPD = MaPNPPD/PunchVol; MaPDNPPE = MaPNPPE/PunchVol; MaPDNPPF = MaPNPPF/PunchVol; MaPDNPPG = MaPNPPG/PunchVol; MaPDNPPH = MaPNPPH/PunchVol; MaPDNPPI = MaPNPPI/PunchVol;
    
    %Expected number of parasites
    ExpParaNPPA = MaPDNPPA * BPVol; ExpParaNPPB = MaPDNPPB * BPVol; ExpParaNPPC = MaPDNPPC * BPVol; ExpParaNPPD = MaPDNPPD * BPVol; ExpParaNPPE = MaPDNPPE * BPVol; ExpParaNPPF = MaPDNPPF * BPVol; ExpParaNPPG = MaPDNPPG * BPVol; ExpParaNPPH = MaPDPPH * BPVol; ExpParaNPPI = MaPDNPPI * BPVol;  
     
    %Number of parasites at a micro sites is a random number from the
    %exponential distribution with the mean from the macro data.
    MiPNPPA = exprnd(ExpParaNPPA); MiPNPPB = exprnd(ExpParaNPPB); MiPNPPC = exprnd(ExpParaNPPC); MiPNPPD = exprnd(ExpParaNPPD); MiPNPPE = exprnd(ExpParaNPPE); MiPNPPF = exprnd(ExpParaNPPF); MiPNPPG = exprnd(ExpParaNPPG); MiPNPPH = exprnd(ExpParaNPPH); MiPNPPI = exprnd(ExpParaNPPI); 
    
    %Concentration of parasites in BP
    ConcNPPA = MiPNPPA/BPVol; ConcNPPB = MiPNPPB/BPVol; ConcNPPC = MiPNPPC/BPVol; ConcNPPD = MiPNPPD/BPVol; ConcNPPE = MiPNPPE/BPVol; ConcNPPF = MiPNPPF/BPVol; ConcNPPG = MiPNPPG/BPVol; ConcNPPH = MiPNPPH/BPVol; ConcNPPI = MiPNPPI/BPVol;   
    
    %Number of parasites picked up by the fly on feeding on the full 0.7 uL
    %feed volume %BP Scale
    FlyParaNPPA(ii) = ConcNPPA * FeedVol; FlyParaNPPB(ii) = ConcNPPB * FeedVol; FlyParaNPPC(ii) = ConcNPPC * FeedVol; FlyParaNPPD(ii) = ConcNPPD * FeedVol; FlyParaNPPE(ii) = ConcNPPE * FeedVol; FlyParaNPPF(ii) = ConcNPPF * FeedVol; FlyParaNPPG(ii) = ConcNPPG * FeedVol; FlyParaNPPH(ii) = ConcNPPH * FeedVol; FlyParaNPPI(ii) = ConcNPPI * FeedVol;  
      
    %%%INFECTION STATUS
    
    %Number of parasites 
    PoTNPPA(ii) = round(FlyParaNPPA(ii)); PoTNPPB(ii) = round(FlyParaNPPB(ii)); PoTNPPC(ii) = round(FlyParaNPPC(ii)); PoTNPPD(ii) = round(FlyParaNPPD(ii)); PoTNPPE(ii) = round(FlyParaNPPE(ii)); PoTNPPF(ii) = round(FlyParaNPPF(ii)); PoTNPPG(ii) = round(FlyParaNPPG(ii)); PoTNPPH(ii) = round(FlyParaNPPH(ii)); PoTNPPI(ii) = round(FlyParaNPPI(ii));   
    if PoTNPPA(ii) > 500 %This is the infection threshold
        
        pNI = 0.75;
        pInfNPPA = 1 - (pNI^PoTNPPA(ii));
        InfProbNPPA(ii) = pInfNPPA;
        
    else
        InfProbNPPA(ii) = 0;
        
    end
    
    rnNPPA = rand;
    
    if rnNPPA < InfProbNPPA(ii);
        InfYesNoNPPA(ii) = 1;
    else
        InfYesNoNPPA(ii) = 0;
    end
    
    if PoTNPPB(ii) > 500 %This is the infection threshold
        
        pNI = 0.75;
        pInfNPPB = 1 - (pNI^PoTNPPB(ii));
        InfProbNPPB(ii) = pInfNPPB;
        
    else
        InfProbNPPB(ii) = 0;
        
    end
    
    rnNPPB = rand;
    
    if rnNPPB < InfProbNPPB(ii);
        InfYesNoNPPB(ii) = 1;
    else
        InfYesNoNPPB(ii) = 0;
    end
    
    if PoTNPPC(ii) > 500 %This is the infection threshold
        
        pNI = 0.75;
        pInfNPPC = 1 - (pNI^PoTNPPC(ii));
        InfProbNPPC(ii) = pInfNPPC;
        
    else
        InfProbNPPC(ii) = 0;
        
    end
    
    rnNPPC = rand;
    
    if rnNPPC < InfProbNPPC(ii);
        InfYesNoNPPC(ii) = 1;
    else
        InfYesNoNPPC(ii) = 0;
    end
    
    if PoTNPPD(ii) > 500 %This is the infection threshold
        
        pNI = 0.75;
        pInfNPPD = 1 - (pNI^PoTNPPD(ii));
        InfProbNPPD(ii) = pInfNPPD;
        
    else
        InfProbNPPD(ii) = 0;
        
    end
    
    rnNPPD = rand;
    
    if rnNPPD < InfProbNPPD(ii);
        InfYesNoNPPD(ii) = 1;
    else
        InfYesNoNPPD(ii) = 0;
    end
    
    if PoTNPPE(ii) > 500 %This is the infection threshold
        
        pNI = 0.75;
        pInfNPPE = 1 - (pNI^PoTNPPE(ii));
        InfProbNPPE(ii) = pInfNPPE;
        
    else
        InfProbNPPE(ii) = 0;
        
    end
    
    rnNPPE = rand;
    
    if rnNPPE < InfProbNPPE(ii);
        InfYesNoNPPE(ii) = 1;
    else
        InfYesNoNPPE(ii) = 0;
    end
    
    if PoTNPPF(ii) > 500 %This is the infection threshold
        
        pNI = 0.75;
        pInfNPPF = 1 - (pNI^PoTNPPF(ii));
        InfProbNPPF(ii) = pInfNPPF;
        
    else
        InfProbNPPF(ii) = 0;
        
    end
    
    rnNPPF = rand;
    
    if rnNPPF < InfProbNPPF(ii);
        InfYesNoNPPF(ii) = 1;
    else
        InfYesNoNPPF(ii) = 0;
    end
    
    if PoTNPPG(ii) > 500 %This is the infection threshold
        
        pNI = 0.75;
        pInfNPPG = 1 - (pNI^PoTNPPG(ii));
        InfProbNPPG(ii) = pInfNPPG;
        
    else
        InfProbNPPG(ii) = 0;
        
    end
    
    rnNPPG = rand;
    
    if rnNPPG < InfProbNPPG(ii);
        InfYesNoNPPG(ii) = 1;
    else
        InfYesNoNPPG(ii) = 0;
    end
    
    if PoTNPPH(ii) > 500 %This is the infection threshold
        
        pNI = 0.75;
        pInfNPPH = 1 - (pNI^PoTNPPH(ii));
        InfProbNPPH(ii) = pInfNPPH;
        
    else
        InfProbNPPH(ii) = 0;
        
    end
    
    rnNPPH = rand;
    
    if rnNPPH < InfProbNPPH(ii);
        InfYesNoNPPH(ii) = 1;
    else
        InfYesNoNPPH(ii) = 0;
    end
    
    if PoTNPPI(ii) > 500 %This is the infection threshold
        
        pNI = 0.75;
        pInfNPPI = 1 - (pNI^PoTNPPI(ii));
        InfProbNPPI(ii) = pInfNPPI;
        
    else
        InfProbNPPI(ii) = 0;
        
    end
    
    rnNPPI = rand;
    
    if rnNPPI < InfProbNPPI(ii);
        InfYesNoNPPI(ii) = 1;
    else
        InfYesNoNPPI(ii) = 0;
    end
    
end

YesNPPA = nnz(InfYesNoNPPA); YesNPPB = nnz(InfYesNoNPPB); YesNPPC = nnz(InfYesNoNPPC); YesNPPD = nnz(InfYesNoNPPE); YesNPPE = nnz(InfYesNoNPPE); YesNPPF = nnz(InfYesNoNPPF); YesNPPG = nnz(InfYesNoNPPG); YesNPPH = nnz(InfYesNoNPPH); YesNPPI = nnz(InfYesNoNPPI);

for ii = 1:1000;    
      
    %%%% MICRO NOT PATCHY %%%%
    
    %Number of parasites at macro site
    %Macro not patchy
    MaPNPNPA = MeanA; MaPNPNPB = MeanB; MaPNPNPC = MeanC;MaPNPNPD = MeanD; MaPNPNPE = MeanE; MaPNPNPF = MeanF; MaPNPNPG = MeanG; MaPNPNPH = MeanH; MaPNPNPI = MeanI;
    
    %Parasite density in punch biopsy
    MaPDNPNPA = MaPNPNPA/PunchVol; MaPDNPNPB = MaPNPNPB/PunchVol; MaPDNPNPC = MaPNPNPC/PunchVol; MaPDNPNPD = MaPNPNPD/PunchVol; MaPDNPNPE = MaPNPNPE/PunchVol; MaPDNPNPF = MaPNPNPF/PunchVol; MaPDNPNPG = MaPNPNPG/PunchVol; MaPDNPNPH = MaPNPNPH/PunchVol; MaPDNPNPI = MaPNPNPI/PunchVol;
    
    %Expected number of parasites
    ExpParaNPNPA = MaPDNPNPA * BPVol; ExpParaNPNPB = MaPDNPNPB * BPVol; ExpParaNPNPC = MaPDNPNPC * BPVol; ExpParaNPNPD = MaPDNPNPD * BPVol; ExpParaNPNPE = MaPDNPNPE * BPVol; ExpParaNPNPF = MaPDNPNPF * BPVol; ExpParaNPNPG = MaPDNPNPG * BPVol; ExpParaNPNPH = MaPDNPNPH * BPVol; ExpParaNPNPI = MaPDNPNPI * BPVol; 
    
    %Concentration of parasites in BP
    ConcNPNPA = ExpParaNPNPA/BPVol; ConcNPNPB = ExpParaNPNPB/BPVol; ConcNPNPC = ExpParaNPNPC/BPVol; ConcNPNPD = ExpParaNPNPD/BPVol; ConcNPNPE = ExpParaNPNPE/BPVol; ConcNPNPF = ExpParaNPNPF/BPVol; ConcNPNPG = ExpParaNPNPG/BPVol; ConcNPNPH = ExpParaNPNPH/BPVol; ConcNPNPI = ExpParaNPNPI/BPVol;     
    
    %Number of parasites picked up by the fly on feeding on the full 0.7 uL
    %feed volume
    FlyParaNPNPA(ii) = ConcNPNPA * FeedVol; FlyParaNPNPB(ii) = ConcNPNPB * FeedVol; FlyParaNPNPC(ii) = ConcNPNPC * FeedVol; FlyParaNPNPD(ii) = ConcNPNPD * FeedVol; FlyParaNPNPE(ii) = ConcNPNPE * FeedVol; FlyParaNPNPF(ii) = ConcNPNPF * FeedVol; FlyParaNPNPG(ii) = ConcNPNPG * FeedVol; FlyParaNPNPH(ii) = ConcNPNPH * FeedVol; FlyParaNPNPI(ii) = ConcNPNPI * FeedVol;      
        
    %%%INFECTION STATUS
    
    %Number of parasites 
    PoTNPNPA(ii) = round(FlyParaNPNPA(ii)); PoTNPNPB(ii) = round(FlyParaNPNPB(ii)); PoTNPNPC(ii) = round(FlyParaNPNPC(ii)); PoTNPNPD(ii) = round(FlyParaNPNPD(ii)); PoTNPNPE(ii) = round(FlyParaNPNPE(ii)); PoTNPNPF(ii) = round(FlyParaNPNPF(ii)); PoTNPNPG(ii) = round(FlyParaNPNPG(ii)); PoTNPNPH(ii) = round(FlyParaNPNPH(ii)); PoTNPNPI(ii) = round(FlyParaNPNPI(ii)); 
    
    if PoTNPNPA(ii) > 500 %This is the infection threshold
        
        pNI = 0.75;
        pInfNPNPA = 1 - (pNI^PoTNPNPA(ii));
        InfProbNPNPA(ii) = pInfNPNPA;
        
    else
        InfProbNPNPA(ii) = 0;
        
    end
    
    rnNPNPA = rand;
    
    if rnNPNPA < InfProbNPNPA(ii);
        InfYesNoNPNPA(ii) = 1;
    else
        InfYesNoNPNPA(ii) = 0;
    end
    
    if PoTNPNPB(ii) > 500 %This is the infection threshold
        
        pNI = 0.75;
        pInfNPNPB = 1 - (pNI^PoTNPNPB(ii));
        InfProbNPNPB(ii) = pInfNPNPB;
        
    else
        InfProbNPNPB(ii) = 0;
        
    end
    
    rnNPNPB = rand;
    
    if rnNPNPB < InfProbNPNPB(ii);
        InfYesNoNPNPB(ii) = 1;
    else
        InfYesNoNPNPB(ii) = 0;
    end
    
    if PoTNPNPC(ii) > 500 %This is the infection threshold
        
        pNI = 0.75;
        pInfNPNPC = 1 - (pNI^PoTNPNPC(ii));
        InfProbNPNPC(ii) = pInfNPNPC;
        
    else
        InfProbNPNPC(ii) = 0;
        
    end
    
    rnNPNPC = rand;
    
    if rnNPNPC < InfProbNPNPC(ii);
        InfYesNoNPNPC(ii) = 1;
    else
        InfYesNoNPNPC(ii) = 0;
    end
    
    if PoTNPNPD(ii) > 500 %This is the infection threshold
        
        pNI = 0.75;
        pInfNPNPD = 1 - (pNI^PoTNPNPD(ii));
        InfProbNPNPD(ii) = pInfNPNPD;
        
    else
        InfProbNPNPD(ii) = 0;
        
    end
    
    rnNPNPD = rand;
    
    if rnNPNPD < InfProbNPNPD(ii);
        InfYesNoNPNPD(ii) = 1;
    else
        InfYesNoNPNPD(ii) = 0;
    end
    
    if PoTNPNPE(ii) > 500 %This is the infection threshold
        
        pNI = 0.75;
        pInfNPNPE = 1 - (pNI^PoTNPNPE(ii));
        InfProbNPNPE(ii) = pInfNPNPE;
        
    else
        InfProbNPNPE(ii) = 0;
        
    end
    
    rnNPNPE = rand;
    
    if rnNPNPE < InfProbNPNPE(ii);
        InfYesNoNPNPE(ii) = 1;
    else
        InfYesNoNPNPE(ii) = 0;
    end
    
    if PoTNPNPF(ii) > 500 %This is the infection threshold
        
        pNI = 0.75;
        pInfNPNPF = 1 - (pNI^PoTNPNPF(ii));
        InfProbNPNPF(ii) = pInfNPNPF;
        
    else
        InfProbNPNPF(ii) = 0;
        
    end
    
    rnNPNPF = rand;
    
    if rnNPNPF < InfProbNPNPF(ii);
        InfYesNoNPNPF(ii) = 1;
    else
        InfYesNoNPNPF(ii) = 0;
    end
    
    if PoTNPNPG(ii) > 500 %This is the infection threshold
        
        pNI = 0.75;
        pInfNPNPG = 1 - (pNI^PoTNPNPG(ii));
        InfProbNPNPG(ii) = pInfNPNPG;
        
    else
        InfProbNPNPG(ii) = 0;
        
    end
    
    rnNPNPG = rand;
    
    if rnNPNPG < InfProbNPNPG(ii);
        InfYesNoNPNPG(ii) = 1;
    else
        InfYesNoNPNPG(ii) = 0;
    end
    
    if PoTNPNPH(ii) > 500 %This is the infection threshold
        
        pNI = 0.75;
        pInfNPNPH = 1 - (pNI^PoTNPNPH(ii));
        InfProbNPNPH(ii) = pInfNPNPH;
        
    else
        InfProbNPNPH(ii) = 0;
        
    end
    
    rnNPNPH = rand;
    
    if rnNPNPH < InfProbNPNPH(ii);
        InfYesNoNPNPH(ii) = 1;
    else
        InfYesNoNPNPH(ii) = 0;
    end
    
    if PoTNPNPI(ii) > 500 %This is the infection threshold
        
        pNI = 0.75;
        pInfNPNPI = 1 - (pNI^PoTNPNPI(ii));
        InfProbNPNPI(ii) = pInfNPNPI;
        
    else
        InfProbNPNPI(ii) = 0;
        
    end
    
    rnNPNPI = rand;
    
    if rnNPNPI < InfProbNPNPI(ii);
        InfYesNoNPNPI(ii) = 1;
    else
        InfYesNoNPNPI(ii) = 0;
    end
    
end

YesNPNPA = nnz(InfYesNoNPNPA); YesNPNPB = nnz(InfYesNoNPNPB); YesNPNPC = nnz(InfYesNoNPNPC); YesNPNPD = nnz(InfYesNoNPNPE); YesNPNPE = nnz(InfYesNoNPNPE); YesNPNPF = nnz(InfYesNoNPNPF); YesNPNPG = nnz(InfYesNoNPNPG); YesNPNPH = nnz(InfYesNoNPNPH); YesNPNPI = nnz(InfYesNoNPNPI);

for ii = 1:1000;    
      
    %%%% MICRO NOT PATCHY %%%%
    
    %Number of parasites at macro site
    %Macro patchy
    MaPPNPA = nbinrnd(MaNBA.r, MaNBA.p); MaPPNPB = nbinrnd(MaNBB.r, MaNBB.p); MaPPNPC = nbinrnd(MaNBC.r, MaNBC.p); MaPPNPD = nbinrnd(MaNBD.r, MaNBD.p); MaPPNPE = nbinrnd(MaNBE.r, MaNBE.p); MaPPNPF = nbinrnd(MaNBF.r, MaNBF.p); MaPPNPG = nbinrnd(MaNBG.r, MaNBG.p); MaPPNPH = nbinrnd(MaNBH.r, MaNBH.p); MaPPNPI = nbinrnd(MaNBI.r, MaNBI.p); 
    
    %Parasite density in punch biopsy
    MaPDPNPA = MaPPNPA/PunchVol; MaPDPNPB = MaPPNPB/PunchVol; MaPDPNPC = MaPPNPC/PunchVol; MaPDPNPD = MaPPNPD/PunchVol; MaPDPNPE = MaPPNPE/PunchVol; MaPDPNPF = MaPPNPF/PunchVol; MaPDPNPG = MaPPNPG/PunchVol; MaPDPNPH = MaPPNPH/PunchVol; MaPDPNPI = MaPPNPI/PunchVol;  
    
    %Expected number of parasites
    ExpParaPNPA = MaPDPNPA * BPVol; ExpParaPNPB = MaPDPNPB * BPVol; ExpParaPNPC = MaPDPNPC * BPVol; ExpParaPNPD = MaPDPNPD * BPVol; ExpParaPNPE = MaPDPNPE * BPVol; ExpParaPNPF = MaPDPNPF * BPVol; ExpParaPNPG = MaPDPNPG * BPVol; ExpParaPNPH = MaPDPNPH * BPVol; ExpParaPNPI = MaPDPNPI * BPVol; 
    
    %Concentration of parasites in BP
    ConcPNPA = ExpParaPNPA/BPVol; ConcPNPB = ExpParaPNPB/BPVol; ConcPNPC = ExpParaPNPC/BPVol; ConcPNPD = ExpParaPNPD/BPVol; ConcPNPE = ExpParaPNPE/BPVol; ConcPNPF = ExpParaPNPF/BPVol; ConcPNPG = ExpParaPNPG/BPVol; ConcPNPH = ExpParaPNPH/BPVol; ConcPNPI = ExpParaPNPI/BPVol;    
    
    %Number of parasites picked up by the fly on feeding on the full 0.7 uL
    %feed volume
    FlyParaPNPA(ii) = ConcPNPA * FeedVol; FlyParaPNPB(ii) = ConcPNPB * FeedVol; FlyParaPNPC(ii) = ConcPNPC * FeedVol; FlyParaPNPD(ii) = ConcPNPD * FeedVol; FlyParaPNPE(ii) = ConcPNPE * FeedVol; FlyParaPNPF(ii) = ConcPNPF * FeedVol; FlyParaPNPG(ii) = ConcPNPG * FeedVol; FlyParaPNPH(ii) = ConcPNPH * FeedVol; FlyParaPNPI(ii) = ConcPNPI * FeedVol;      
        
    %%%INFECTION STATUS
    
    %Number of parasites 
    PoTPNPA(ii) = round(FlyParaPNPA(ii)); PoTPNPB(ii) = round(FlyParaPNPB(ii)); PoTPNPC(ii) = round(FlyParaPNPC(ii)); PoTPNPD(ii) = round(FlyParaPNPD(ii)); PoTPNPE(ii) = round(FlyParaPNPE(ii)); PoTPNPF(ii) = round(FlyParaPNPF(ii)); PoTPNPG(ii) = round(FlyParaPNPG(ii)); PoTPNPH(ii) = round(FlyParaPNPH(ii)); PoTPNPI(ii) = round(FlyParaPNPI(ii)); 
    
    if PoTPNPA(ii) > 500 %This is the infection threshold
        
        pNI = 0.75;
        pInfPNPA = 1 - (pNI^PoTPNPA(ii));
        InfProbPNPA(ii) = pInfPNPA;
        
    else
        InfProbPNPA(ii) = 0;
        
    end
    
    rnPNPA = rand;
    
    if rnPNPA < InfProbPNPA(ii);
        InfYesNoPNPA(ii) = 1;
    else
        InfYesNoPNPA(ii) = 0;
    end
    
    if PoTPNPB(ii) > 500 %This is the infection threshold
        
        pNI = 0.75;
        pInfPNPB = 1 - (pNI^PoTPNPB(ii));
        InfProbPNPB(ii) = pInfPNPB;
        
    else
        InfProbPNPB(ii) = 0;
        
    end
    
    rnPNPB = rand;
    
    if rnPNPB < InfProbPNPB(ii);
        InfYesNoPNPB(ii) = 1;
    else
        InfYesNoPNPB(ii) = 0;
    end
    
    if PoTPNPC(ii) > 500 %This is the infection threshold
        
        pNI = 0.75;
        pInfPNPC = 1 - (pNI^PoTPNPC(ii));
        InfProbPNPC(ii) = pInfPNPC;
        
    else
        InfProbPNPC(ii) = 0;
        
    end
    
    rnPNPC = rand;
    
    if rnPNPC < InfProbPNPC(ii);
        InfYesNoPNPC(ii) = 1;
    else
        InfYesNoPNPC(ii) = 0;
    end
    
    if PoTPNPD(ii) > 500 %This is the infection threshold
        
        pNI = 0.75;
        pInfPNPD = 1 - (pNI^PoTPNPD(ii));
        InfProbPNPD(ii) = pInfPNPD;
        
    else
        InfProbPNPD(ii) = 0;
        
    end
    
    rnPNPD = rand;
    
    if rnPNPD < InfProbPNPD(ii);
        InfYesNoPNPD(ii) = 1;
    else
        InfYesNoPNPD(ii) = 0;
    end
    
    if PoTPNPE(ii) > 500 %This is the infection threshold
        
        pNI = 0.75;
        pInfPNPE = 1 - (pNI^PoTPNPE(ii));
        InfProbPNPE(ii) = pInfPNPE;
        
    else
        InfProbPNPE(ii) = 0;
        
    end
    
    rnPNPE = rand;
    
    if rnPNPE < InfProbPNPE(ii);
        InfYesNoPNPE(ii) = 1;
    else
        InfYesNoPNPE(ii) = 0;
    end
    
    if PoTPNPF(ii) > 500 %This is the infection threshold
        
        pNI = 0.75;
        pInfPNPF = 1 - (pNI^PoTPNPF(ii));
        InfProbPNPF(ii) = pInfPNPF;
        
    else
        InfProbPNPF(ii) = 0;
        
    end
    
    rnPNPF = rand;
    
    if rnPNPF < InfProbPNPF(ii);
        InfYesNoPNPF(ii) = 1;
    else
        InfYesNoPNPF(ii) = 0;
    end
    
    if PoTPNPG(ii) > 500 %This is the infection threshold
        
        pNI = 0.75;
        pInfPNPG = 1 - (pNI^PoTPNPG(ii));
        InfProbPNPG(ii) = pInfPNPG;
        
    else
        InfProbPNPG(ii) = 0;
        
    end
    
    rnPNPG = rand;
    
    if rnPNPG < InfProbPNPG(ii);
        InfYesNoPNPG(ii) = 1;
    else
        InfYesNoPNPG(ii) = 0;
    end
    
    if PoTPNPH(ii) > 500 %This is the infection threshold
        
        pNI = 0.75;
        pInfPNPH = 1 - (pNI^PoTPNPH(ii));
        InfProbPNPH(ii) = pInfPNPH;
        
    else
        InfProbPNPH(ii) = 0;
        
    end
    
    rnPNPH = rand;
    
    if rnPNPH < InfProbPNPH(ii);
        InfYesNoPNPH(ii) = 1;
    else
        InfYesNoPNPH(ii) = 0;
    end
    
    if PoTPNPI(ii) > 500 %This is the infection threshold
        
        pNI = 0.75;
        pInfPNPI = 1 - (pNI^PoTPNPI(ii));
        InfProbPNPI(ii) = pInfPNPI;
        
    else
        InfProbPNPI(ii) = 0;
        
    end
    
    rnPNPI = rand;
    
    if rnPNPI < InfProbPNPI(ii);
        InfYesNoPNPI(ii) = 1;
    else
        InfYesNoPNPI(ii) = 0;
    end
    
end

YesPNPA = nnz(InfYesNoPNPA); YesPNPB = nnz(InfYesNoPNPB); YesPNPC = nnz(InfYesNoPNPC); YesPNPD = nnz(InfYesNoPNPE); YesPNPE = nnz(InfYesNoPNPE); YesPNPF = nnz(InfYesNoPNPF); YesPNPG = nnz(InfYesNoPNPG); YesPNPH = nnz(InfYesNoPNPH); YesPNPI = nnz(InfYesNoPNPI); 

Patchiness = {'MaP MiP';'MaNP MiP';'MaNP MiNP';'MaP MiNP'};
data = [YesPPA YesPPB YesPPC YesPPD YesPPE YesPPF YesPPG YesPPH YesPPI; YesNPPA YesNPPB YesNPPC YesNPPD YesNPPE YesNPPF YesNPPG YesNPPH YesNPPI; YesNPNPA YesNPNPB YesNPNPC YesNPNPD YesNPNPE YesNPNPF YesNPNPG YesNPNPH YesNPNPI; YesPNPA YesPNPB YesPNPC YesPNPD YesPNPE YesPNPF YesPNPG YesPNPH YesPNPI]
bar(data,'grouped')
set(gca,'xticklabel',Patchiness)
ylabel('Number of Infected Flies')
title('Number of Infected Flies on the BP Scale as a Function of Patchiness - RPMI Group')
legend('MouseF Exp2','MouseJ Exp2','Rag3 Exp4','Rag10 Exp4','Rag3 Exp5','Rag4 Exp5','Rag7 Exp5','Rag9 Exp5','Rag11 Exp5')