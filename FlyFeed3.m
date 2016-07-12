%Punch biopsy data inserted here
Biopsies = [113652.8672	23657.07813	7039.90332	67137.04688	64529.57031	131615.2969	518365.5938	541672.125	1800.070679	14785.81934	1394753	706560.6875	130131.8047	23352.04688	208169.5938	608231.375	106441.1953	32235.32422	44494.23438	94896.08594	6424.298828	5645.563477	27523.33594	32097.01172];
MeanB = mean(Biopsies);

Biop = reshape(Biopsies, [numel(Biopsies),1]);

B = round(Biop);

%Fit Negative Binomial Distribution. Parameters are MaNB.r and MaNB.p
MaNB = fitdist(B, 'Negative Binomial');

%Volume of punch biopsy in m^3
PunchVol = 5.03 * 10^-8 ;

%Volume of a single plane of BP. Cuboidal BP split into 4micron thick
%(parasite width) planes
BPPlaneVol = 2.8 * 10^-12; %Think you can ignore this

%Number of planes in BP. 1.4mm / 4um = 350
BPPlanes = 350; %Think you can ignore this too...!

%Total BP volume of 1 uL (cuboidal dimensions = 1.4mm x 1.4mm x 500um)
BPVol = 1 * 10^-9;

%Proboscis Volume. Cylindrical dimensions = 75 um diameter, 350 um length).
%Cuboidal dimensions = 6.65 * 10^-5 x 6.65 * 10^-5 x 350 * 10^-6 m
ProbVol = 1.55 * 10^-12;

%Number of planes in Prob = 6.65*10^-5 / 4*10^-6
ProbPlanes = 16.6; %Ignore this too...

%Volume of a single proboscis plane
ProbPlaneVol = 9.31 * 10^-14; %And ignore this...!

%Feed volume of fly
FeedVol = 0.7*10^-9;

for ii = [1:1000];
    
    %%%% MICRO-PATCHINESS %%%%
    
    %Number of parasites at macro site
    MaP = nbinrnd(MaNB.r, MaNB.p);  %Use if macro is patchy
    %MaP = MeanB;   %Use if micro is not patchy
    
    %Parasite density in punch biopsy
    MaPD = MaP/PunchVol;
    
    %Expected number of parasites
    %ExpPara = MaPD * ProbVol;
    ExpPara = MaPD * BPVol;
    
     
    %Number of parasites at a micro sites is a random number from the
    %exponential distribution with the mean from the macro data.
    MiP = exprnd(ExpPara); 
    
    
    %Concentration of parasites in BP
    Conc = MiP/BPVol; %Use for BP scale
    
    %Number of parasites picked up by the fly on feeding on the full 0.7 uL
    %feed volume
    FlyPara(ii) = Conc * FeedVol; %Use for BP scale
    %FlyPara(ii) = MiP; %Use for prob scale
    
%     %%%% MICRO NOT PATCHY
%     
%     %Number of parasites at macro site
%     %MaP = nbinrnd(MaNB.r, MaNB.p);
%     MaP = MeanB;
%     
%     %Parasite density in punch biopsy
%     MaPD = MaP/PunchVol;
%     
%     %Expected number of parasites
%     %ExpPara = MaPD * BPVol;
%     ExpPara = MaPD * ProbVol;
%     
%     %Concentration of parasites in BP
%     %Conc = ExpPara/BPVol; %Use for BP scale
%     
%     %Number of parasites picked up by the fly on feeding on the full 0.7 uL
%     %feed volume
%     %FlyPara(ii) = Conc * FeedVol; %Use for BP scale
%     FlyPara(ii) = ExpPara; %Use for prob scale
    
%     
    %%%INFECTION STATUS
    
    %Number of parasites 
    PoT(ii) = round(FlyPara(ii));
    
    if PoT(ii) > 500 %This is the infection threshold
        
        pNI = 0.75;
        pInf = 1 - (pNI^PoT(ii));
        InfProb(ii) = pInf;

        
    else
        InfProb(ii) = 0;
        
    end
    
    rn = rand;
    
    if rn < InfProb(ii);
        InfYesNo(ii) = 1;
    else
        InfYesNo(ii) = 0;
    end
end

Yes = nnz(InfYesNo)  %number of nonzero elements in matrix InfYesNo
No = numel(InfYesNo) - Yes; %number of elements in matrix InfYesNo
Bar = [No, Yes];
x = [0,1];

figure
subplot(2,1,1)
histogram(PoT, 'FaceColor', 'r')
title('Number of Parasites ingested')
xlabel('Number of parasites')
ylabel('Number of flies')
subplot(2,1,2)
bar(x,Bar, 'FaceColor', 'g')
ax = gca
ax.XTickLabel = {'Not Infected', 'Infected'}
ylabel('Number of Flies')

