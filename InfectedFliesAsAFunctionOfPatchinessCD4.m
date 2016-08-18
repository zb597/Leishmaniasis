%Mice
MA = [6242.309	40213.92	31013.05	208242.7	60042.1	98883.61	60476.59	249944.2	10892.03	6057.031	112096.5	35391.51	86337.55	45613.8	46020.12	64171.25	56225.69	20209.33	55936.79	42591.23	88499.5	5941.761	5517.263	139639];
MB = [113652.9	23657.08	7039.903	67137.05	64529.57	131615.3	518365.6	541672.1	1800.071	14785.82	1394753	706560.7	130131.8	23352.05	208169.6	608231.4	106441.2	32235.32	44494.23	94896.09	6424.299	5645.563	27523.34	32097.01];
M14 = [18169.2	12718.1	16364.7	34639.3	102337.6	79552.5	115456.1	262826.1	81365.2	110573.1	209569.7	67009.6	66966.9	29948.9	232116.4	766783.9	63812.4	36401	77071.4	56991.7	14903.7	60864.6	8198.5	1496.9];
M2 = [44545.2	530171.5	668891.6	2469231	19219.4	571615.5	1491697	725753	42283.1	263689.3	1826778	318708.8	164573.9	91659.4	462878.5	60530.8	7163.7	47056.6	50569.9	38927.7	3931.7	6870.1	23407.3	26510.7];
M64 = [598196	1926035	488947.2	173848.7	553694.9	4061678	14662020	3776093	90795.3	3336097	3694323	9669794	293533.2	3504126	10284320	5941550	450755.2	4809105	3664616	1397573	12728.2	49876.7	67301.9	70181.8];
M8 = [3137887	6669272	1751551	2126678	6912118	5712355	15197510	17344610	7932178	4256139	80643520	11343150	8134356	6343389	7060130	6448552	3770326	2546905	6616066	6271906	185902.1	74053.5	248866.6	341170.5];
M15 = [43606.68	39610.4	5482.782	31639.34	32049.93	33995.91	37138.89	172406.8	28486.92	33673.63	50312.51	12606.38	123821	152871.1	83785.74	27238.4	476785.3	79134.71	36514.55	45777.34	13625.87	20583.27	4612.97	140612.8];
M65 = [2260.624	9063.668	3060.94	9444.282	4412.582	3038.998	2265.217	970.2748	7934.078	4313.613	7064.852	1514.515	15138.66	8959.973	9887.268	790.0161	7823.779	2604.275	6803.212	2141.352	4012.497	8171.049	906.79	4226.915];
M10 = [34873.48	178186	24155.73	12365.47	15526.27	22582.66	38169.07	6815.063	4903.425	36461.53	36999.52	27371.52	23405.3	130282.7	31532.11	41119.96	22138.27	42285.54	67661.07	53227.32	6173.713	7604.093	31383.42	45624.58];

MeanA = mean(MA); MeanB = mean(MB); MeanC = mean(M14); MeanD = mean(M2); MeanE = mean(M64); MeanF = mean(M8); MeanG = mean(M15); MeanH = mean(M65); MeanI = mean(M10); 
A = round(reshape(MA, [numel(MA),1])); B = round(reshape(MB, [numel(MB),1])); C = round(reshape(M14, [numel(M14),1])); D = round(reshape(M2, [numel(M2),1])); E = round(reshape(M64, [numel(M64),1])); F = round(reshape(M8, [numel(M8),1])); G = round(reshape(M15, [numel(M15),1])); H = round(reshape(M65, [numel(M65),1])); I = round(reshape(M10, [numel(M10),1]));

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
title('Number of Infected Flies on the BP Scale as a Function of Patchiness - CD4 Group')
legend('MouseA Exp2','MouseB Exp2','Rag1 Exp4','Rag2 Exp4','Rag6 Exp4','Rag8 Exp4','Rag1 Exp5','Rag6 Exp5','Rag10 Exp5')