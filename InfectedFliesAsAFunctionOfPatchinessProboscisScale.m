%Mice
BiopA = [6242.308594 40213.92188 31013.05078 208242.6563 60042.10156 98883.60938 60476.59375 249944.2188 10892.03223 6057.03125	112096.5313	35391.50781	86337.54688	45613.80078	46020.12109	64171.25 56225.69141 20209.32813 55936.78906 42591.23438 88499.5 5941.76123	5517.262695	139638.9688];
BiopB = [113652.8672 23657.07813 7039.90332	67137.04688	64529.57031	131615.2969	518365.5938	541672.125	1800.070679	14785.81934	1394753	706560.6875	130131.8047	23352.04688	208169.5938	608231.375	106441.1953	32235.32422	44494.23438	94896.08594	6424.298828	5645.563477	27523.33594	32097.01172];
BiopC = [86988.09375 2483.02417	14647.36816	48430.07813	10765.91016	9871.273438	9871.273438	28048.60938	201290.25	13951.59961	42517.27734	971090.375	2863765.75	36039.0625	42434.35156	179498.3125	220141.0469	10395.6875	40059.78906	52496.25781	44046.84375	6026.121582	2963.452393	2269.223633];
BiopD = [3636.382568 9097.695313 3001.909668 23878.79297 918.046875	4182.071777	1520.401367	1117.206055	2155.180176	5595.969727	1236.572754	17303.19531	290.6083984	1136.777832	13201.63477	1235.423096	1470.355835	4630.273926	6996.893555	2384.648926	22086.73828	3055.782959	370.6409302	16257.79492];
BiopE = [15409.66016 38189.30078 7407.382813 12525.3877	14765.42969	48938.80859	11563.5166	16625.5625	6094.720703	3692.590332	41796.57813	14807.84375	9591.882813	3691.366699	39645.85156	9301.40625	29549.90625	12585.88477	3411.382324	41044.53906	2181.679932	5054.815918	25538.76953	3942.522705];
BiopF = [1918.327393 5368.237793 14975.54102 3805.636719 4451.460938 17625.02344 27372.53125 17554.89453 250.9739227 11183.27734 7323.999512 7664.599121 1906.455933 5972.289063 5060.821289 11847.68164 7574.185547 717.4968262 1669.182861 10734.91504 43123.44922 2413.856445 731.732666	42438.03906];
BiopG = [5879.079102 22170.35547 6642.881836 175939.375	29015.39844	8576.245117	15615.9082	20745.52734	1759.902832	2210.585449	8686.671875	7552.5625	4470.313965	4072.462891	10441.4375	5817.073242	12959.82227	4377.237305	11262.03516	14850.03125	1024.440186	2940.30249	6098.556641	6186.19043];
BiopH = [7838.777344 4216.459473 4139.382813 13826.62305 1222.324097 14711.82324 14711.82324 27832.98828 35677 2495.137695 5517.897461 13922.68457 9070.527344 2650.269531 8553.928711 17811.33398 7258.19043 15696.14648 2292.208496 7610.755859 50638.11719 11762.1748 1077.937134 32804.46484];
BiopI = [1637.440186 1071.902344 1039.78186	1847.271973	1918.505127	696.3065186	3563.663086	2304.054688	2853.534668	674.1523438	3040.928223	322.5037231	2841.975586	509.4277344	514.2527466	1855.881958	816.8228149	1061.693604	1126.812866	1587.541748	695.0185547	215.0773773	571.5085449	7049.145996];
BiopJ = [15997.08789 3017.408203 5034.854004 9239.501953 2739.726074 513.6808472 1732.656006 1191.046387 1036.687012 3348.614746 1548.505249 860.7709961	7333.506836	1555.587158	2160.678223	3284.401367	3297.141357	922.8627319	4877.167969	6981.65332	555.0819092	763.1016846	27997.11719	45156.42969];

MeanA = mean(BiopA); MeanB = mean(BiopB); MeanC = mean(BiopC); MeanD = mean(BiopD); MeanE = mean(BiopE); MeanF = mean(BiopF); MeanG = mean(BiopG); MeanH = mean(BiopH); MeanI = mean(BiopI); MeanJ = mean(BiopJ);
A = round(reshape(BiopA, [numel(BiopA),1])); B = round(reshape(BiopB, [numel(BiopB),1])); C = round(reshape(BiopC, [numel(BiopC),1])); D = round(reshape(BiopD, [numel(BiopD),1])); E = round(reshape(BiopE, [numel(BiopE),1])); F = round(reshape(BiopF, [numel(BiopF),1])); G = round(reshape(BiopG, [numel(BiopG),1])); H = round(reshape(BiopH, [numel(BiopH),1])); I = round(reshape(BiopI, [numel(BiopI),1])); J = round(reshape(BiopJ, [numel(BiopJ),1]));

%Fit Negative Binomial Distribution.
MaNBA = fitdist(A, 'Negative Binomial'); MaNBB = fitdist(B, 'Negative Binomial'); MaNBC = fitdist(C, 'Negative Binomial'); MaNBD = fitdist(D, 'Negative Binomial'); MaNBE = fitdist(E, 'Negative Binomial'); MaNBF = fitdist(F, 'Negative Binomial'); MaNBG = fitdist(G, 'Negative Binomial'); MaNBH = fitdist(H, 'Negative Binomial'); MaNBI = fitdist(I, 'Negative Binomial'); MaNBJ = fitdist(J, 'Negative Binomial');

%Volume of punch biopsy in m^3
PunchVol = 5.03 * 10^-8 ;

%Proboscis Volume. Cylindrical dimensions = 75 um diameter, 350 um length).
%Cuboidal dimensions = 6.65 * 10^-5 x 6.65 * 10^-5 x 350 * 10^-6 m
ProbVol = 1.55 * 10^-12;

%Feed volume of fly
FeedVol = 0.7*10^-9;

for ii = 1:1000;
    
    %%%% MICRO-PATCHINESS %%%%
    
    %Number of parasites at macro site
    %Use if macro is patchy
    MaPPPA = nbinrnd(MaNBA.r, MaNBA.p); MaPPPB = nbinrnd(MaNBB.r, MaNBB.p); MaPPPC = nbinrnd(MaNBC.r, MaNBC.p); MaPPPD = nbinrnd(MaNBD.r, MaNBD.p); MaPPPE = nbinrnd(MaNBE.r, MaNBE.p); MaPPPF = nbinrnd(MaNBF.r, MaNBF.p); MaPPPG = nbinrnd(MaNBG.r, MaNBG.p); MaPPPH = nbinrnd(MaNBH.r, MaNBH.p); MaPPPI = nbinrnd(MaNBI.r, MaNBI.p); MaPPPJ = nbinrnd(MaNBJ.r, MaNBJ.p);
        
    %Parasite density in punch biopsy
    MaPDPPA = MaPPPA/PunchVol; MaPDPPB = MaPPPB/PunchVol; MaPDPPC = MaPPPC/PunchVol; MaPDPPD = MaPPPD/PunchVol; MaPDPPE = MaPPPE/PunchVol; MaPDPPF = MaPPPF/PunchVol; MaPDPPG = MaPPPG/PunchVol; MaPDPPH = MaPPPH/PunchVol; MaPDPPI = MaPPPI/PunchVol; MaPDPPJ = MaPPPJ/PunchVol;
    
    %Expected number of parasites
    ExpParaPPA = MaPDPPA * ProbVol; ExpParaPPB = MaPDPPB * ProbVol; ExpParaPPC = MaPDPPC * ProbVol; ExpParaPPD = MaPDPPD * ProbVol; ExpParaPPE = MaPDPPE * ProbVol; ExpParaPPF = MaPDPPF * ProbVol; ExpParaPPG = MaPDPPG * ProbVol; ExpParaPPH = MaPDPPH * ProbVol; ExpParaPPI = MaPDPPI * ProbVol; ExpParaPPJ = MaPDPPJ * ProbVol; 
     
    %Number of parasites at a micro sites is a random number from the
    %exponential distribution with the mean from the macro data.
    MiPPPA = exprnd(ExpParaPPA); MiPPPB = exprnd(ExpParaPPB); MiPPPC = exprnd(ExpParaPPC); MiPPPD = exprnd(ExpParaPPD); MiPPPE = exprnd(ExpParaPPE); MiPPPF = exprnd(ExpParaPPF); MiPPPG = exprnd(ExpParaPPG); MiPPPH = exprnd(ExpParaPPH); MiPPPI = exprnd(ExpParaPPI); MiPPPJ = exprnd(ExpParaPPJ);
    
    %Number of parasites picked up by the fly on feeding on the full 0.7 uL
    %feed volume %BP Scale
    FlyParaPPA(ii) = MiPPPA; FlyParaPPB(ii) = MiPPPB; FlyParaPPC(ii) = MiPPPC; FlyParaPPD(ii) = MiPPPD; FlyParaPPE(ii) = MiPPPE; FlyParaPPF(ii) = MiPPPF; FlyParaPPG(ii) = MiPPPG; FlyParaPPH(ii) = MiPPPH; FlyParaPPI(ii) = MiPPPI; FlyParaPPJ(ii) = MiPPPJ; 
      
    %%%INFECTION STATUS
    
    %Number of parasites 
    PoTPPA(ii) = round(FlyParaPPA(ii)); PoTPPB(ii) = round(FlyParaPPB(ii)); PoTPPC(ii) = round(FlyParaPPC(ii)); PoTPPD(ii) = round(FlyParaPPD(ii)); PoTPPE(ii) = round(FlyParaPPE(ii)); PoTPPF(ii) = round(FlyParaPPF(ii)); PoTPPG(ii) = round(FlyParaPPG(ii)); PoTPPH(ii) = round(FlyParaPPH(ii)); PoTPPI(ii) = round(FlyParaPPI(ii));  PoTPPJ(ii) = round(FlyParaPPJ(ii));
    
    if PoTPPA(ii) > 1 %This is the infection threshold
        
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
    
    if PoTPPB(ii) > 1 %This is the infection threshold
        
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
    
    if PoTPPC(ii) > 1 %This is the infection threshold
        
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
    
    if PoTPPD(ii) > 1 %This is the infection threshold
        
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
    
    if PoTPPE(ii) > 1 %This is the infection threshold
        
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
    
    if PoTPPF(ii) > 1 %This is the infection threshold
        
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
    
    if PoTPPG(ii) > 1 %This is the infection threshold
        
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
    
    if PoTPPH(ii) > 1 %This is the infection threshold
        
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
    
    if PoTPPI(ii) > 1 %This is the infection threshold
        
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
    
    if PoTPPJ(ii) > 1 %This is the infection threshold
        
        pNI = 0.75;
        pInfPPJ = 1 - (pNI^PoTPPJ(ii));
        InfProbPPJ(ii) = pInfPPJ;
        
    else
        InfProbPPJ(ii) = 0;
        
    end
    
    rnPPJ = rand;
    
    if rnPPJ < InfProbPPJ(ii);
        InfYesNoPPJ(ii) = 1;
    else
        InfYesNoPPJ(ii) = 0;
    end
end

YesPPA = nnz(InfYesNoPPA); YesPPB = nnz(InfYesNoPPB); YesPPC = nnz(InfYesNoPPC); YesPPD = nnz(InfYesNoPPE); YesPPE = nnz(InfYesNoPPE); YesPPF = nnz(InfYesNoPPF); YesPPG = nnz(InfYesNoPPG); YesPPH = nnz(InfYesNoPPH); YesPPI = nnz(InfYesNoPPI); YesPPJ = nnz(InfYesNoPPJ);

for ii = 1:1000;
    
    %%%% MICRO-PATCHINESS %%%%
    
    %Number of parasites at macro site   
    %Use if macro is not patchy
    MaPNPPA = MeanA; MaPNPPB = MeanB; MaPNPPC = MeanC; MaPNPPD = MeanD; MaPNPPE = MeanE; MaPNPPF = MeanF; MaPNPPG = MeanG; MaPNPPH = MeanH; MaPNPPI = MeanI; MaPNPPJ = MeanJ;
    
    %Parasite density in punch biopsy
    MaPDNPPA = MaPNPPA/PunchVol; MaPDNPPB = MaPNPPB/PunchVol; MaPDNPPC = MaPNPPC/PunchVol; MaPDNPPD = MaPNPPD/PunchVol; MaPDNPPE = MaPNPPE/PunchVol; MaPDNPPF = MaPNPPF/PunchVol; MaPDNPPG = MaPNPPG/PunchVol; MaPDNPPH = MaPNPPH/PunchVol; MaPDNPPI = MaPNPPI/PunchVol; MaPDNPPJ = MaPNPPJ/PunchVol;
    
    %Expected number of parasites
    ExpParaNPPA = MaPDNPPA * ProbVol; ExpParaNPPB = MaPDNPPB * ProbVol; ExpParaNPPC = MaPDNPPC * ProbVol; ExpParaNPPD = MaPDNPPD * ProbVol; ExpParaNPPE = MaPDNPPE * ProbVol; ExpParaNPPF = MaPDNPPF * ProbVol; ExpParaNPPG = MaPDNPPG * ProbVol; ExpParaNPPH = MaPDNPPH * ProbVol; ExpParaNPPI = MaPDNPPI * ProbVol; ExpParaNPPJ = MaPDNPPJ * ProbVol; 
     
    %Number of parasites at a micro sites is a random number from the
    %exponential distribution with the mean from the macro data.
    MiPNPPA = exprnd(ExpParaNPPA); MiPNPPB = exprnd(ExpParaNPPB); MiPNPPC = exprnd(ExpParaNPPC); MiPNPPD = exprnd(ExpParaNPPD); MiPNPPE = exprnd(ExpParaNPPE); MiPNPPF = exprnd(ExpParaNPPF); MiPNPPG = exprnd(ExpParaNPPG); MiPNPPH = exprnd(ExpParaNPPH); MiPNPPI = exprnd(ExpParaNPPI); MiPNPPJ = exprnd(ExpParaNPPJ);
    
    %Number of parasites picked up by the fly on feeding on the full 0.7 uL
    %feed volume %BP Scale
    FlyParaNPPA(ii) = MiPNPPA; FlyParaNPPB(ii) = MiPNPPB; FlyParaNPPC(ii) = MiPNPPC; FlyParaNPPD(ii) = MiPNPPD; FlyParaNPPE(ii) = MiPNPPE; FlyParaNPPF(ii) = MiPNPPF; FlyParaNPPG(ii) = MiPNPPG; FlyParaNPPH(ii) = MiPNPPH; FlyParaNPPI(ii) = MiPNPPI; FlyParaNPPJ(ii) = MiPNPPJ; 
      
    %%%INFECTION STATUS
    
    %Number of parasites 
    PoTNPPA(ii) = round(FlyParaNPPA(ii)); PoTNPPB(ii) = round(FlyParaNPPB(ii)); PoTNPPC(ii) = round(FlyParaNPPC(ii)); PoTNPPD(ii) = round(FlyParaNPPD(ii)); PoTNPPE(ii) = round(FlyParaNPPE(ii)); PoTNPPF(ii) = round(FlyParaNPPF(ii)); PoTNPPG(ii) = round(FlyParaNPPG(ii)); PoTNPPH(ii) = round(FlyParaNPPH(ii)); PoTNPPI(ii) = round(FlyParaNPPI(ii));  PoTNPPJ(ii) = round(FlyParaNPPJ(ii));
    
    if PoTNPPA(ii) > 1 %This is the infection threshold
        
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
    
    if PoTNPPB(ii) > 1 %This is the infection threshold
        
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
    
    if PoTNPPC(ii) > 1 %This is the infection threshold
        
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
    
    if PoTNPPD(ii) > 1 %This is the infection threshold
        
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
    
    if PoTNPPE(ii) > 1 %This is the infection threshold
        
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
    
    if PoTNPPF(ii) > 1 %This is the infection threshold
        
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
    
    if PoTNPPG(ii) > 1 %This is the infection threshold
        
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
    
    if PoTNPPH(ii) > 1 %This is the infection threshold
        
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
    
    if PoTNPPI(ii) > 1 %This is the infection threshold
        
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
    
    if PoTNPPJ(ii) > 1 %This is the infection threshold
        
        pNI = 0.75;
        pInfNPPJ = 1 - (pNI^PoTNPPJ(ii));
        InfProbNPPJ(ii) = pInfNPPJ;
        
    else
        InfProbNPPJ(ii) = 0;
        
    end
    
    rnNPPJ = rand;
    
    if rnNPPJ < InfProbNPPJ(ii);
        InfYesNoNPPJ(ii) = 1;
    else
        InfYesNoNPPJ(ii) = 0;
    end
end

YesNPPA = nnz(InfYesNoNPPA); YesNPPB = nnz(InfYesNoNPPB); YesNPPC = nnz(InfYesNoNPPC); YesNPPD = nnz(InfYesNoNPPE); YesNPPE = nnz(InfYesNoNPPE); YesNPPF = nnz(InfYesNoNPPF); YesNPPG = nnz(InfYesNoNPPG); YesNPPH = nnz(InfYesNoNPPH); YesNPPI = nnz(InfYesNoNPPI); YesNPPJ = nnz(InfYesNoNPPJ);

for ii = 1:1000;    
      
    %%%% MICRO NOT PATCHY %%%%
    
    %Number of parasites at macro site
    %Macro not patchy
    MaPNPNPA = MeanA; MaPNPNPB = MeanB; MaPNPNPC = MeanC;MaPNPNPD = MeanD; MaPNPNPE = MeanE; MaPNPNPF = MeanF; MaPNPNPG = MeanG; MaPNPNPH = MeanH; MaPNPNPI = MeanI; MaPNPNPJ = MeanJ;
    
    %Parasite density in punch biopsy
    MaPDNPNPA = MaPNPNPA/PunchVol; MaPDNPNPB = MaPNPNPB/PunchVol; MaPDNPNPC = MaPNPNPC/PunchVol; MaPDNPNPD = MaPNPNPD/PunchVol; MaPDNPNPE = MaPNPNPE/PunchVol; MaPDNPNPF = MaPNPNPF/PunchVol; MaPDNPNPG = MaPNPNPG/PunchVol; MaPDNPNPH = MaPNPNPH/PunchVol; MaPDNPNPI = MaPNPNPI/PunchVol;  MaPDNPNPJ = MaPNPNPJ/PunchVol;
    
    %Expected number of parasites
    ExpParaNPNPA = MaPDNPNPA * ProbVol; ExpParaNPNPB = MaPDNPNPB * ProbVol; ExpParaNPNPC = MaPDNPNPC * ProbVol; ExpParaNPNPD = MaPDNPNPD * ProbVol; ExpParaNPNPE = MaPDNPNPE * ProbVol; ExpParaNPNPF = MaPDNPNPF * ProbVol; ExpParaNPNPG = MaPDNPNPG * ProbVol; ExpParaNPNPH = MaPDNPNPH * ProbVol; ExpParaNPNPI = MaPDNPNPI * ProbVol; ExpParaNPNPJ = MaPDNPNPJ * ProbVol;
    
    %Number of parasites picked up by the fly on feeding on the full 0.7 uL
    %feed volume
    FlyParaNPNPA(ii) = ExpParaNPNPA; FlyParaNPNPB(ii) = ExpParaNPNPB; FlyParaNPNPC(ii) = ExpParaNPNPC; FlyParaNPNPD(ii) = ExpParaNPNPD; FlyParaNPNPE(ii) = ExpParaNPNPE; FlyParaNPNPF(ii) = ExpParaNPNPF; FlyParaNPNPG(ii) = ExpParaNPNPG; FlyParaNPNPH(ii) = ExpParaNPNPH; FlyParaNPNPI(ii) = ExpParaNPNPI; FlyParaNPNPJ(ii) = ExpParaNPNPJ;     
        
    %%%INFECTION STATUS
    
    %Number of parasites 
    PoTNPNPA(ii) = round(FlyParaNPNPA(ii)); PoTNPNPB(ii) = round(FlyParaNPNPB(ii)); PoTNPNPC(ii) = round(FlyParaNPNPC(ii)); PoTNPNPD(ii) = round(FlyParaNPNPD(ii)); PoTNPNPE(ii) = round(FlyParaNPNPE(ii)); PoTNPNPF(ii) = round(FlyParaNPNPF(ii)); PoTNPNPG(ii) = round(FlyParaNPNPG(ii)); PoTNPNPH(ii) = round(FlyParaNPNPH(ii)); PoTNPNPI(ii) = round(FlyParaNPNPI(ii)); PoTNPNPJ(ii) = round(FlyParaNPNPJ(ii));
    
    if PoTNPNPA(ii) > 1 %This is the infection threshold
        
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
    
    if PoTNPNPB(ii) > 1 %This is the infection threshold
        
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
    
    if PoTNPNPC(ii) > 1 %This is the infection threshold
        
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
    
    if PoTNPNPD(ii) > 1 %This is the infection threshold
        
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
    
    if PoTNPNPE(ii) > 1 %This is the infection threshold
        
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
    
    if PoTNPNPF(ii) > 1 %This is the infection threshold
        
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
    
    if PoTNPNPG(ii) > 1 %This is the infection threshold
        
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
    
    if PoTNPNPH(ii) > 1 %This is the infection threshold
        
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
    
    if PoTNPNPI(ii) > 1 %This is the infection threshold
        
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
    
    if PoTNPNPJ(ii) > 1 %This is the infection threshold
        
        pNI = 0.75;
        pInfNPNPJ = 1 - (pNI^PoTNPNPJ(ii));
        InfProbNPNPJ(ii) = pInfNPNPJ;
        
    else
        InfProbNPNPJ(ii) = 0;
        
    end
    
    rnNPNPJ = rand;
    
    if rnNPNPJ < InfProbNPNPJ(ii);
        InfYesNoNPNPJ(ii) = 1;
    else
        InfYesNoNPNPJ(ii) = 0;
    end
end

YesNPNPA = nnz(InfYesNoNPNPA); YesNPNPB = nnz(InfYesNoNPNPB); YesNPNPC = nnz(InfYesNoNPNPC); YesNPNPD = nnz(InfYesNoNPNPE); YesNPNPE = nnz(InfYesNoNPNPE); YesNPNPF = nnz(InfYesNoNPNPF); YesNPNPG = nnz(InfYesNoNPNPG); YesNPNPH = nnz(InfYesNoNPNPH); YesNPNPI = nnz(InfYesNoNPNPI); YesNPNPJ = nnz(InfYesNoNPNPJ);

for ii = 1:1000;    
      
    %%%% MICRO NOT PATCHY %%%%
    
    %Number of parasites at macro site
    %Macro patchy
    MaPPNPA = nbinrnd(MaNBA.r, MaNBA.p); MaPPNPB = nbinrnd(MaNBB.r, MaNBB.p); MaPPNPC = nbinrnd(MaNBC.r, MaNBC.p); MaPPNPD = nbinrnd(MaNBD.r, MaNBD.p); MaPPNPE = nbinrnd(MaNBE.r, MaNBE.p); MaPPNPF = nbinrnd(MaNBF.r, MaNBF.p); MaPPNPG = nbinrnd(MaNBG.r, MaNBG.p); MaPPNPH = nbinrnd(MaNBH.r, MaNBH.p); MaPPNPI = nbinrnd(MaNBI.r, MaNBI.p); MaPPNPJ = nbinrnd(MaNBJ.r, MaNBJ.p);
    
    %Parasite density in punch biopsy
    MaPDPNPA = MaPPNPA/PunchVol; MaPDPNPB = MaPPNPB/PunchVol; MaPDPNPC = MaPPNPC/PunchVol; MaPDPNPD = MaPPNPD/PunchVol; MaPDPNPE = MaPPNPE/PunchVol; MaPDPNPF = MaPPNPF/PunchVol; MaPDPNPG = MaPPNPG/PunchVol; MaPDPNPH = MaPPNPH/PunchVol; MaPDPNPI = MaPPNPI/PunchVol;  MaPDPNPJ = MaPPNPJ/PunchVol;
    
    %Expected number of parasites
    ExpParaPNPA = MaPDPNPA * ProbVol; ExpParaPNPB = MaPDPNPB * ProbVol; ExpParaPNPC = MaPDPNPC * ProbVol; ExpParaPNPD = MaPDPNPD * ProbVol; ExpParaPNPE = MaPDPNPE * ProbVol; ExpParaPNPF = MaPDPNPF * ProbVol; ExpParaPNPG = MaPDPNPG * ProbVol; ExpParaPNPH = MaPDPNPH * ProbVol; ExpParaPNPI = MaPDPNPI * ProbVol; ExpParaPNPJ = MaPDPNPJ * ProbVol;
    
    %Number of parasites picked up by the fly on feeding on the full 0.7 uL
    %feed volume
    FlyParaPNPA(ii) = ExpParaPNPA; FlyParaPNPB(ii) = ExpParaPNPB; FlyParaPNPC(ii) = ExpParaPNPC; FlyParaPNPD(ii) = ExpParaPNPD; FlyParaPNPE(ii) = ExpParaPNPE; FlyParaPNPF(ii) = ExpParaPNPF; FlyParaPNPG(ii) = ExpParaPNPG; FlyParaPNPH(ii) = ExpParaPNPH; FlyParaPNPI(ii) = ExpParaPNPI; FlyParaPNPJ(ii) = ExpParaPNPJ;     
        
    %%%INFECTION STATUS
    
    %Number of parasites 
    PoTNPNPA(ii) = round(FlyParaNPNPA(ii)); PoTNPNPB(ii) = round(FlyParaNPNPB(ii)); PoTNPNPC(ii) = round(FlyParaNPNPC(ii)); PoTNPNPD(ii) = round(FlyParaNPNPD(ii)); PoTNPNPE(ii) = round(FlyParaNPNPE(ii)); PoTNPNPF(ii) = round(FlyParaNPNPF(ii)); PoTNPNPG(ii) = round(FlyParaNPNPG(ii)); PoTNPNPH(ii) = round(FlyParaNPNPH(ii)); PoTNPNPI(ii) = round(FlyParaNPNPI(ii)); PoTNPNPJ(ii) = round(FlyParaNPNPJ(ii));
    
    %%%INFECTION STATUS
    
    %Number of parasites 
    PoTPNPA(ii) = round(FlyParaPNPA(ii)); PoTPNPB(ii) = round(FlyParaPNPB(ii)); PoTPNPC(ii) = round(FlyParaPNPC(ii)); PoTPNPD(ii) = round(FlyParaPNPD(ii)); PoTPNPE(ii) = round(FlyParaPNPE(ii)); PoTPNPF(ii) = round(FlyParaPNPF(ii)); PoTPNPG(ii) = round(FlyParaPNPG(ii)); PoTPNPH(ii) = round(FlyParaPNPH(ii)); PoTPNPI(ii) = round(FlyParaPNPI(ii)); PoTPNPJ(ii) = round(FlyParaPNPJ(ii));
    
    if PoTPNPA(ii) > 1 %This is the infection threshold
        
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
    
    if PoTPNPB(ii) > 1 %This is the infection threshold
        
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
    
    if PoTPNPC(ii) > 1 %This is the infection threshold
        
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
    
    if PoTPNPD(ii) > 1 %This is the infection threshold
        
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
    
    if PoTPNPE(ii) > 1 %This is the infection threshold
        
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
    
    if PoTPNPF(ii) > 1 %This is the infection threshold
        
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
    
    if PoTPNPG(ii) > 1 %This is the infection threshold
        
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
    
    if PoTPNPH(ii) > 1 %This is the infection threshold
        
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
    
    if PoTPNPI(ii) > 1 %This is the infection threshold
        
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
    
    if PoTPNPJ(ii) > 1 %This is the infection threshold
        
        pNI = 0.75;
        pInfPNPJ = 1 - (pNI^PoTPNPJ(ii));
        InfProbPNPJ(ii) = pInfPNPJ;
        
    else
        InfProbPNPJ(ii) = 0;
        
    end
    
    rnPNPJ = rand;
    
    if rnPNPJ < InfProbPNPJ(ii);
        InfYesNoPNPJ(ii) = 1;
    else
        InfYesNoPNPJ(ii) = 0;
    end
end

YesPNPA = nnz(InfYesNoPNPA); YesPNPB = nnz(InfYesNoPNPB); YesPNPC = nnz(InfYesNoPNPC); YesPNPD = nnz(InfYesNoPNPE); YesPNPE = nnz(InfYesNoPNPE); YesPNPF = nnz(InfYesNoPNPF); YesPNPG = nnz(InfYesNoPNPG); YesPNPH = nnz(InfYesNoPNPH); YesPNPI = nnz(InfYesNoPNPI); YesPNPJ = nnz(InfYesNoPNPJ);

Patchiness = {'MaP MiP';'MaNP MiP';'MaNP MiNP';'MaP MiNP'};
data = [YesPPA YesPPB YesPPC YesPPD YesPPE YesPPF YesPPG YesPPH YesPPI YesPPJ; YesNPPA YesNPPB YesNPPC YesNPPD YesNPPE YesNPPF YesNPPG YesNPPH YesNPPI YesNPPJ; YesNPNPA YesNPNPB YesNPNPC YesNPNPD YesNPNPE YesNPNPF YesNPNPG YesNPNPH YesNPNPI YesNPNPJ; YesPNPA YesPNPB YesPNPC YesPNPD YesPNPE YesPNPF YesPNPG YesPNPH YesPNPI YesPNPJ]
bar(data,'grouped')
set(gca,'xticklabel',Patchiness)
ylabel('Number of Infected Flies')
title('Number of Infected Flies on the Proboscis Scale as a Function of Patchiness')