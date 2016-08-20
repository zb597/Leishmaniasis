%Mice
Rag1 = [1734.9	4606.1	4255.6	5866.0	8419.2	1115.6	915.1	813.4	1905.9	7063.1	915.2	1245.8	1602.3	8424.4	1362.0	2523.7	1140.7	8686.8	12205.9	2051.1	3327.6	1340.7	980.6	173.5];
Rag2 = [5581.7	1476.5	756.9	2162.6	16668.0	422.8	1995.1	11084.2	2355.6	6956.3	221.8	2708.2	1839.7	1675.0	2105.4	737.8	249.1	2177.8	470.5	167.1	3667.1	1235.0	656.6	3548.2];
Rag3 = [505.5	173.2	1884.3	5646.4	1433.9	1635.7	800.3	1670.8	16082.9	3466.9	1376.7	524.0	1751.0	4847.2	3003.9	1316.0	1356.7	2526.7	2671.0	5941.0	7748.1	6120.4	3849.7	714.5];
Rag4 = [1082.3	841.0	787.7	680.5	1473.9	66.4	578.4	521.7	755.7	529.4	781.4	781.1	1660.2	1170.1	1677.1	567.9	1291.2	1171.7	1797.6	700.2	2031.9	1318.6	1630.4	1447.0];
Rag5 = [868.3	303.7	2543.6	248.4	12076.7	98.8	246.4	413.2	1042.1	685.7	792.8	1656.2	1659.3	717.2	2846.7	2748.0	623.6	564.2	917.4	2703.1	1084.0	1620.8	402.9	129.7];
Rag6 = [18815.3	61973.8	37314.1	2225.0	144193.4	28475.9	122334.3	15205.1	6420.3	4582.2	11234.9	47808.0	15934.2	3369.7	22128.1	5687.2	7504.0	1921.5	3141.7	1870.7	34220.9	5318.9	899.5	48252.0];
Rag7 = [81094.5	397965.5	95943.1	189913.2	54714.1	53272.7	820945.0	1489117.0	76049.7	133271.6	374300.6	943975.9	19880.9	78146.6	369051.2	237573.3	23747.6	262029.2	433661.6	1001859.9	1429.9	15965.9	15016.3	69149.0];
Rag8 = [266994.7	612286.0	587857.0	578527.3	37708.9	192540.2	204166.2	349983.6	4876.1	61548.5	228418.5	54736.7	77599.3	107644.2	159764.8	34644.1	37963.3	15894.1	154566.9	132574.3	5059.9	2870.5	157096.4	98691.8];
Rag9 = [12413.2	10824.3	3444.4	8597.4	703.2	20057.6	82067.7	12291.2	16895.1	22202.2	26330.9	2902.4	63653.7	41884.9	77456.6	6736.4	5400.5	24222.7	35832.7	16875.7	318.0	572.1	996.2	2235.6];

MeanA = mean(Rag1); MeanB = mean(Rag2); MeanC = mean(Rag3); MeanD = mean(Rag4); MeanE = mean(Rag5); MeanF = mean(Rag6); MeanG = mean(Rag7); MeanH = mean(Rag8); MeanI = mean(Rag9); 
A = round(reshape(Rag1, [numel(Rag1),1])); B = round(reshape(Rag2, [numel(Rag2),1])); C = round(reshape(Rag3, [numel(Rag3),1])); D = round(reshape(Rag4, [numel(Rag4),1])); E = round(reshape(Rag5, [numel(Rag5),1])); F = round(reshape(Rag6, [numel(Rag6),1])); G = round(reshape(Rag7, [numel(Rag7),1])); H = round(reshape(Rag8, [numel(Rag8),1])); I = round(reshape(Rag9, [numel(Rag9),1]));

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

Patchiness = {'MaP MiP';'MaNP MiP'};
data = [YesPPA YesPPB YesPPC YesPPD YesPPE YesPPF YesPPG YesPPH YesPPI; YesNPPA YesNPPB YesNPPC YesNPPD YesNPPE YesNPPF YesNPPG YesNPPH YesNPPI]
bar(data,'grouped')
set(gca,'xticklabel',Patchiness)
ylabel('Number of Infected Flies')
title('Number of Infected Flies on the BP Scale as a Function of Patchiness - Untreated Group')
legend('Rag1','Rag2','Rag3','Rag4','Rag5','Rag6','Rag7','Rag8','Rag9')