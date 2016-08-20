%Mice
MC = [86988.09	2483.024	14647.37	48430.08	10765.91	9871.273	9871.273	28048.61	201290.3	13951.6	42517.28	971090.4	2863766	36039.06	42434.35	179498.3	220141	10395.69	40059.79	52496.26	44046.84	6026.122	2963.452	2269.224];
MI = [1637.44	1071.902	1039.782	1847.272	1918.505	696.3065	3563.663	2304.055	2853.535	674.1523	3040.928	322.5037	2841.976	509.4277	514.2527	1855.882	816.8228	1061.694	1126.813	1587.542	695.0186	215.0774	571.5085	7049.146];
M4 = [98782.1	23432.9	547804.1	218793.6	30498.7	42836.4	132063.2	154165	22458.6	30615.8	68735.9	56174.5	13688	30541.9	123152.8	22664.2	166397.5	24192	19657.9	42091.4	4579.2	852.9	6435.1	4415.2];
M54 = [181818.8	109668.8	412396.3	157043.2	10899.3	37303.5	52511.8	200723.1	75878.3	25746.7	16414.3	29792	96075.4	242080	28506	33780.8	46730.9	4845.2	19985.8	23996	16646.6	14755.1	4314.1	30883.5];
M7 = [867842	552477.6	388486.9	1449490	164748.1	489848.7	374590.3	1272215	65345.5	16331.7	132568.7	218877.4	119276.6	184811.2	1266134	2730995	865704.1	3535952	5601397	4343749	5967.6	10030.2	764994.3	72690];
M9 = [46414.6	219530.5	52236.3	128603.4	62248.6	19690.2	22252.3	79189.2	11349.9	14101.1	12007.4	170496.8	55922.1	158562.8	456310.2	499857.7	332860.1	99346.4	216622.8	79926.2	9206.4	12650.2	3191.1	2906.4];
M2 = [53771.23	15734.01	12067.08	35575.91	28530.69	12811.91	10118.98	13885.84	7139.758	16043.8	9844.204	6656.57	14711.64	56160.36	132494	48295.63	87786.71	68436.54	36531.82	72108.79	12688.16	14848.2	8130.676	10732.51];
M55 = [2099.204	849.3257	3339.803	2380.775	3682.561	1291.839	5476.861	1360.699	1888.724	4222.9	2481.203	1476.645	3769.24	1287.345	1493.629	1151.788	2676.89	3870.346	3575.469	6870.077	901.0186	1427.349	1226.367	1000.991];
M8 = [28611.06	3129.208	3214.273	4881.184	3551.307	17673.08	47348.29	15514.55	860.3441	3260.179	64844.03	3.775513	17518.26	19609.43	3.540197	7105.975	6733.547	9861.456	11472.41	11072.78	4183.498	14446.12	9036.298	11701.71];

MeanA = mean(MC); MeanB = mean(MI); MeanC = mean(M4); MeanD = mean(M54); MeanE = mean(M7); MeanF = mean(M9); MeanG = mean(M2); MeanH = mean(M55); MeanI = mean(M8); 
A = round(reshape(MC, [numel(MC),1])); B = round(reshape(MI, [numel(MI),1])); C = round(reshape(M4, [numel(M4),1])); D = round(reshape(M54, [numel(M54),1])); E = round(reshape(M7, [numel(M7),1])); F = round(reshape(M9, [numel(M9),1])); G = round(reshape(M2, [numel(M2),1])); H = round(reshape(M55, [numel(M55),1])); I = round(reshape(M8, [numel(M8),1]));

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
title('Number of Infected Flies on the BP Scale as a Function of Patchiness - CD8 Group')
legend('MouseC Exp2','MouseI Exp2','Rag4 Exp4','Rag5 Exp4','Rag7 Exp4','Rag9 Exp4','Rag2 Exp5','Rag5 Exp5','Rag8 Exp5')