%Punch biopsy data inserted here
Biopsies = [113652.8672	23657.07813	7039.90332	67137.04688	64529.57031	131615.2969	518365.5938	541672.125	1800.070679	14785.81934	1394753	706560.6875	130131.8047	23352.04688	208169.5938	608231.375	106441.1953	32235.32422	44494.23438	94896.08594	6424.298828	5645.563477	27523.33594	32097.01172];
MeanB = mean(Biopsies);

Biop = reshape(Biopsies, [numel(Biopsies),1]);

B = round(Biop);

%Fit Negative Binomial Distribution. Parameters are MaNB.r and MaNB.p
MaNB = fitdist(B, 'Negative Binomial');

%Volume of punch biopsy in m^3
PunchVol = 5.03 * 10^-8 ;

%Total BP volume of 1 uL (cuboidal dimensions = 1.4mm x 1.4mm x 500um)
BPVol = 1 * 10^-9;

%Proboscis Volume. Cylindrical dimensions = 75 um diameter, 350 um length).
%Cuboidal dimensions = 6.65 * 10^-5 x 6.65 * 10^-5 x 350 * 10^-6 m
ProbVol = 1.55 * 10^-12;

%Feed volume of fly
FeedVol = 0.7*10^-9;

for ii = 1:1000;
        
    %%%% MICRO NOT PATCHY
    
    %Number of parasites at macro site
    %MaP = nbinrnd(MaNB.r, MaNB.p);
    MaP = MeanB;
    
    %Parasite density in punch biopsy
    MaPD = MaP/PunchVol;
    
    %Expected number of parasites
    ExpPara = MaPD * BPVol;
    ExpParaMi = MaPD * ProbVol;
    
    %Concentration of parasites in BP
    Conc = ExpPara/BPVol; %Use for BP scale
    
    %Number of parasites picked up by the fly on feeding on the full 0.7 uL
    %feed volume
    FlyPara(ii) = Conc * FeedVol; %Use for BP scale
    FlyParaMi(ii) = ExpParaMi; %Use for prob scale
         
    %%%INFECTION STATUS1
    
    %Number of parasites
    PoT(ii) = round(FlyPara(ii));
    
    if PoT(ii) > 1 %This is the infection threshold
        
        pNI = 0.75;
        pInf = 1 - (pNI^PoT(ii));
        InfProb(ii) = pInf;
        
    else
        InfProb(ii) = 0;
        
    end
    
    rn11 = rand;
    
    if rn11 < InfProb(ii);
        InfYesNo11(ii) = 1;
    else
        InfYesNo11(ii) = 0;
    end
    
     if PoT(ii) > 10 %This is the infection threshold
        
        pNI = 0.75;
        pInf = 1 - (pNI^PoT(ii));
        InfProb(ii) = pInf;
        
    else
        InfProb(ii) = 0;
        
    end
    
    rn21 = rand;
    
    if rn21 < InfProb(ii);
        InfYesNo21(ii) = 1;
    else
        InfYesNo21(ii) = 0;
    end
    
    if PoT(ii) > 50 %This is the infection threshold
        
        pNI = 0.75;
        pInf = 1 - (pNI^PoT(ii));
        InfProb(ii) = pInf;
        
    else
        InfProb(ii) = 0;
        
    end
    
    rn31 = rand;
    
    if rn31 < InfProb(ii);
        InfYesNo31(ii) = 1;
    else
        InfYesNo31(ii) = 0;
    end
    
    if PoT(ii) > 100 %This is the infection threshold
        
        pNI = 0.75;
        pInf = 1 - (pNI^PoT(ii));
        InfProb(ii) = pInf;
        
    else
        InfProb(ii) = 0;
        
    end
    
    rn41 = rand;
    
    if rn41 < InfProb(ii);
        InfYesNo41(ii) = 1;
    else
        InfYesNo41(ii) = 0;
    end
    
    if PoT(ii) > 500 %This is the infection threshold
        
        pNI = 0.75;
        pInf = 1 - (pNI^PoT(ii));
        InfProb(ii) = pInf;
        
    else
        InfProb(ii) = 0;
        
    end
    
    rn51 = rand;
    
    if rn51 < InfProb(ii);
        InfYesNo51(ii) = 1;
    else
        InfYesNo51(ii) = 0;
    end
    
    if PoT(ii) > 1000 %This is the infection threshold
        
        pNI = 0.75;
        pInf = 1 - (pNI^PoT(ii));
        InfProb(ii) = pInf;
        
    else
        InfProb(ii) = 0;
        
    end
    
    rn61 = rand;
    
    if rn61 < InfProb(ii);
        InfYesNo61(ii) = 1;
    else
        InfYesNo61(ii) = 0;
    end
    
    %%%INFECTION STATUS2
    
    %Number of parasites
    PoTMi(ii) = round(FlyParaMi(ii));
    
    if PoTMi(ii) > 1 %This is the infection threshold
        
        pNI = 0.75;
        pInfMi = 1 - (pNI^PoTMi(ii));
        InfProbMi(ii) = pInfMi;
        
    else
        InfProbMi(ii) = 0;
        
    end
    
    rn12 = rand;
    
    if rn12 < InfProbMi(ii);
        InfYesNo12(ii) = 1;
    else
        InfYesNo12(ii) = 0;
    end
    
     if PoTMi(ii) > 10 %This is the infection threshold
        
        pNI = 0.75;
        pInfMi = 1 - (pNI^PoTMi(ii));
        InfProbMi(ii) = pInfMi;
        
    else
        InfProbMi(ii) = 0;
        
    end
    
    rn22 = rand;
    
    if rn22 < InfProbMi(ii);
        InfYesNo22(ii) = 1;
    else
        InfYesNo22(ii) = 0;
    end
    
    if PoTMi(ii) > 50 %This is the infection threshold
        
        pNI = 0.75;
        pInfMi = 1 - (pNI^PoTMi(ii));
        InfProbMi(ii) = pInfMi;
        
    else
        InfProbMi(ii) = 0;
        
    end
    
    rn32 = rand;
    
    if rn32 < InfProbMi(ii);
        InfYesNo32(ii) = 1;
    else
        InfYesNo32(ii) = 0;
    end
    
    if PoTMi(ii) > 100 %This is the infection threshold
        
        pNI = 0.75;
        pInfMi = 1 - (pNI^PoTMi(ii));
        InfProbMi(ii) = pInfMi;
        
    else
        InfProbMi(ii) = 0;
        
    end
    
    rn42 = rand;
    
    if rn42 < InfProbMi(ii);
        InfYesNo42(ii) = 1;
    else
        InfYesNo42(ii) = 0;
    end
    
    if PoTMi(ii) > 500 %This is the infection threshold
        
        pNI = 0.75;
        pInfMi = 1 - (pNI^PoTMi(ii));
        InfProbMi(ii) = pInfMi;
        
    else
        InfProbMi(ii) = 0;
        
    end
    
    rn52 = rand;
    
    if rn52 < InfProbMi(ii);
        InfYesNo52(ii) = 1;
    else
        InfYesNo52(ii) = 0;
    end
    
    if PoTMi(ii) > 1000 %This is the infection threshold
        
        pNI = 0.75;
        pInfMi = 1 - (pNI^PoTMi(ii));
        InfProbMi(ii) = pInfMi;
        
    else
        InfProbMi(ii) = 0;
        
    end
    
    rn62 = rand;
    
    if rn62 < InfProbMi(ii);
        InfYesNo62(ii) = 1;
    else
        InfYesNo62(ii) = 0;
    end
end

Yes11 = nnz(InfYesNo11);  %number of nonzero elements in matrix InfYesNo11
No11 = numel(InfYesNo11) - Yes11; %number of elements in matrix InfYesNo11
Yes21 = nnz(InfYesNo21);
No21 = numel(InfYesNo21) - Yes21;
Yes31 = nnz(InfYesNo31);
No31 = numel(InfYesNo31) - Yes31;
Yes41 = nnz(InfYesNo41);
No41 = numel(InfYesNo41) - Yes41;
Yes51 = nnz(InfYesNo51);
No51 = numel(InfYesNo51) - Yes51;
Yes61 = nnz(InfYesNo61);
No61 = numel(InfYesNo61) - Yes61;

Yes12 = nnz(InfYesNo12);  %number of nonzero elements in matrix InfYesNo12
No12 = numel(InfYesNo12) - Yes12; %number of elements in matrix InfYesNo12
Yes22 = nnz(InfYesNo22);
No22 = numel(InfYesNo22) - Yes22;
Yes32 = nnz(InfYesNo32);
No32 = numel(InfYesNo32) - Yes32;
Yes42 = nnz(InfYesNo42);
No42 = numel(InfYesNo42) - Yes42;
Yes52 = nnz(InfYesNo52);
No52 = numel(InfYesNo52) - Yes52;
Yes62 = nnz(InfYesNo62);
No62 = numel(InfYesNo62) - Yes62;

IT = {'>1';'>10';'>50';'>100';'>500';'>1000'};
data = [Yes11 Yes12; Yes21 Yes22; Yes31 Yes32; Yes41 Yes42; Yes51 Yes52; Yes61 Yes62]
bar(data,'grouped')
set(gca,'xticklabel',IT)
legend({'BP Scale';'Prob Scale'})
xlabel('Infection Threshold')
ylabel('Number of Infected Flies')
title('Number of Infected Flies as a Function of the Threshold')