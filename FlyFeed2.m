%Number of parasites in the skin biopsies
Biop = [1918.327393	5368.237793	14975.54102	3805.636719	4451.460938	17625.02344	27372.53125	17554.89453	250.9739227	11183.27734	7323.999512	7664.599121	1906.455933	5972.289063	5060.821289	11847.68164	7574.185547	717.4968262	1669.182861	10734.91504	43123.44922	2413.856445	731.732666	42438.03906]; 
MB = mean(Biop);
MeanBiop = mean(Biop)/1772.5;

Data = reshape(Biop, [numel(Biop), 1]);
D = round(Data);
NB = fitdist(D, 'Negative Binomial');
M = (NB.r*(1-NB.p))/NB.p;

%Site Area. For proboscis scale, this is 2.625 x 10^-8. For BP scale this
%is 7.05 x 10^-7.
SiteArea = 7.05*10^-7;
%SiteArea = 2.625 * 10^-8;
%SiteArea = 1 * 10^-9; % 1 uL blood pool volume

%Feed size
FeedSize = 6*10^-7;
%FeedSize = 0.7 * 10^-9; % 0.7 uL feed volume

%Numbers of parasites picked up by flies
FlyPara = zeros(1,100);

%Number of infected flies
InfFlies = zeros(1,100);

MSPara = zeros(1,100);

P = zeros(1,100);

InfYesNo = zeros(1,100);

for ii = (1:100);
    
    %To find number of bites before hitting capillary
    Squares = [0,1,0,0,1,0,0,1,0,0,1,0,0,0,0,0,1,0,0,0,1,0,1,1,0,0,0,0,0,1];
    %Squares = [1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1];
    
    %attempts
    c = 1;
    
    y = datasample(Squares, 1, 'Replace', true);
    
    if y == 1; %If a capillary hit on first attempt, all well and good!
            %disp('Capillary')
            %disp(c)
            

     else y == 0; %Otherwise, it keeps sampling until it hits a capillary
            while y == 0;
                c = c + 1;
                y = datasample(Squares,1,'Replace', true);
            end
    end
    
    attempts(ii) = c;
    
    %%%MICROSCALE PATCHINESS
          
    %Number of parasites at macro-site (MaS)
%     Ma = nbinrnd(NB.r, NB.p);
%     MaNB(ii) = Ma/1772.5;
    %MaNB(ii) = (datasample(Biop,1))/1772.5;
    
    MaS = (datasample(Biop,1))/1772.5;
    %MaS = MeanBiop;
    
    %Macro site density (MSD)
    MaSD = MaS / (7*10^-6);
    %MaSD = MaS / (2.84*10^-11);
    
    %Average number of parasites at micro site. This number is the mean for
    %the exponential distribution that is formed for the microscale
    %patchiness.
    MSPara(ii) = (MaSD * SiteArea);
    
    %Pick an actual number of parasites for the micropatch from the
    %exponential distribution
    P(ii) = exprnd(MSPara(ii));
    %f = exprnd(MSPara(ii), [1,7]);
    %P(ii) = sum(f);
    %P(ii) = datasample(BPBiopsy, 1);
    %P(ii) = MSPara(ii);
    
    %Find the concentration of parasites in the blood pool formed
    Conc = P(ii)/SiteArea; 
    
    %Number of parasites picked up is the concentration x the size of the
    %feed.
    FlyPara(ii) = round(Conc * FeedSize); 
    %FlyPara(ii) = P(ii);

%     %%MICROSCALE NOT PATCHY
%     %BiopsyPara = (datasample(Biop,1))/2225;
%     BiopsyPara = MeanBiop;
%     
%     MaSD = BiopsyPara / (8.9 * 10^-6);
%     
%     FlyPara(ii) = MaSD * FeedSize;
%     %FlyPara(ii) = MaSD * SiteArea;

    %INFECTION STATUS
    
    %Probability of being infected by the parasite burden (n) is determined
    %by P(Inf) = 1 - P(NotInf)^n.
    
    NV = FlyPara(ii);
    
        if NV > 1
            
            %Probability of not being infected by 1 fly
            pNI = 0.01;
    
            %Probability of being infected by parasite burden
            PInf = 1 - ((pNI)^NV);
            
            InfFlies(ii) = PInf;
            %InfYesNo(ii) = 1;
        else
            InfFlies(ii) = 0;
            %InfYesNo(ii) = 0;
        end
        
    rn = rand;
        
    if rn < InfFlies(ii)
        InfYesNo(ii) = 1;
    else
        InfYesNo(ii) = 0;
    end
        
        
end
  
Yes = nnz(InfYesNo);
No = numel(InfYesNo) - Yes;
Bar = [No, Yes];
x = [0,1];

figure
subplot(4,1,1)
histogram(attempts, 'FaceColor', 'b')
title('Number of Bite Attempts before Capillary is located', 'FontSize', 18)
xlabel('Number of Bite Attempts', 'FontSize', 14)
ylabel('Frequency', 'FontSize', 14)

subplot(4,1,2)
histogram(FlyPara, 'FaceColor', 'r')
title('Number of Parasites picked up by Flies on feeding from Blood Pool', 'FontSize', 18)
xlabel('Number of Parasites', 'FontSize', 14)
ylabel('Frequency', 'FontSize', 14)

subplot(4,1,3)
histogram(InfFlies, 'FaceColor', 'y')
title('Probability that Flies become infected by their parasite burden', 'FontSize', 18)
xlabel('Probability of Infection', 'FontSize', 14)
ylabel('Frequency', 'FontSize', 14)

subplot(4,1,4)
bar(x, Bar, 'FaceColor', 'g')
title('Number of Infected Flies', 'FontSize', 18)
xlabel('Uninfected (0) or Infected (1)', 'FontSize', 14)
ylabel('Frequency', 'FontSize', 14)

figure
histogram(P)
title('Number of Parasites at Micro-site')
xlabel('Number of Parasites')
ylabel('Number of sites')


%Mu = dispersion/probability
%r = dispersion
%p = probability
%NB(r,p)
%I know mu from macro scale stuff, p is ?fixed?! so dispersion is mu x
%probability?!
