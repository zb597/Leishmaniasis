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

Patchiness = {'MaNP';'MaP'};
data = [YesNPNPA YesNPNPB YesNPNPC YesNPNPD YesNPNPE YesNPNPF YesNPNPG YesNPNPH YesNPNPI; YesPNPA YesPNPB YesPNPC YesPNPD YesPNPE YesPNPF YesPNPG YesPNPH YesPNPI]
bar(data,'grouped')
set(gca,'xticklabel',Patchiness)
ylabel('Number of Infected Flies')
title('Number of Infected Flies on the BP Scale as a Function of Patchiness - Untreated Group')
legend('Rag1','Rag2','Rag3','Rag4','Rag5','Rag6','Rag7','Rag8','Rag9')
