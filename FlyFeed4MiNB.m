%Mean number of parasites in a biopsy
Biopsies = [86988.09375
2483.02417
14647.36816
48430.07813
10765.91016
9871.273438
9871.273438
28048.60938
201290.25
13951.59961
42517.27734
971090.375
2863765.75
36039.0625
42434.35156
179498.3125
220141.0469
10395.6875
40059.78906
52496.25781
44046.84375
6026.121582
2963.452393
2269.223633];

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

%Micro Site sampling data
%BP scale
%MicroSamples = [122	219	121	2	7	9	189	60	2];
%Prob scale
MicroSamples = [7	0	4	4	2	0	33	4	9	4	29	1	0	0	0	0	0	0	0	0	0	0	0	0	0	0	1	0	0	1	0	0	4	0	2	0	5	0	1	5	2	5	0	0	0	0	0	0	0	0];

%Fit distribution to micro site data to find dispersion parameter, k, which
%in this is MiNB.r.
MicroSamp = reshape(MicroSamples, [numel(MicroSamples), 1]);
MiNB = fitdist(MicroSamp, 'Negative Binomial');

for ii = 1:1000;
    
    %%%% MICRO-PATCHINESS %%%%
    
    %Number of parasites at macro site
    %MaP = mean(Biopsies);
    MaP = nbinrnd(MaNB.r, MaNB.p);   
    
    %Parasite density in punch biopsy
    MaPD = MaP/PunchVol;
    
    %Expected number of parasites
    ExpPara = MaPD * ProbVol;
    %ExpPara = MaPD * BPVol;
    
    %k = the dispersion parameter. Equivalent to R
    k = MiNB.r*1000;
    
    %Parameter P is equivalent to mean/(mean + disp parameter)
    p = ExpPara/(ExpPara + k);
    
    %Number of parasites at a micro sites is a random number from the
    %distribution with the mean from the macro data and dispersion
    %parameter from the distribution fitted to the micro data. Work out in
    %a single plane first, then multiply by the number of planes to get
    %number of parasites in total BP
    
    MiP = nbinrnd(k, p); 
    
    %Concentration of parasites in BP
    Conc = MiP/BPVol; %Use for BP scale
    
    %Number of parasites picked up by the fly on feeding on the full 0.7 uL
    %feed volume
    %FlyPara(ii) = Conc * FeedVol; %Use for BP scale
    FlyPara(ii) = MiP; %Use for prob scale    

    %%%INFECTION STATUS
    
    %Number of parasites over time assuming that half will die before
    %infection can be assessed...
    PoT(ii) = round(FlyPara(ii));
    
    if PoT(ii) > 500
        
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

Yes = nnz(InfYesNo)
No = numel(InfYesNo) - Yes;
Bar = [No, Yes];
x = [0,1];

% x2 = 0:10;
% y2 = exppdf(x2, ExpPara);
% %y3 = nbinpdf(x2, MaNB.r, MaNB.p);
% % y4 = nbinpdf(x2, k/10, p);
% figure
% subplot(2,1,1)
% histogram(PoT, 'FaceColor', 'r')
% subplot(2,1,2)
% bar(x,Bar, 'FaceColor', 'g')
% % subplot(3,1,3)
% plot(x2,y2); hold on;
% % plot(x2, y3); hold on;
% % % plot(x2, y4)
