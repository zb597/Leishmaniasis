%Number of parasites in the skin biopsies
Biop = [1918.327393	5368.237793	14975.54102	3805.636719	4451.460938	17625.02344	27372.53125	17554.89453	250.9739227	11183.27734	7323.999512	7664.599121	1906.455933	5972.289063	5060.821289	11847.68164	7574.185547	717.4968262	1669.182861	10734.91504	43123.44922	2413.856445	731.732666	42438.03906]; 
MB = mean(Biop); %mean of elements of Biop
MeanBiop = mean(Biop)/1772.5; %1772.5 = average #parasites in a splice of skin
Data = round(reshape(Biop, [numel(Biop), 1]));
NB = fitdist(Data, 'Negative Binomial'); %find parameter values
M = (NB.r*(1-NB.p))/NB.p; %same as mean(Biop) using the parameters from neg.bin.

SiteArea = 7.05*10^-7; %Blood pool scale
FeedSize = 6*10^-7;
FlyPara = zeros(1,100); %#parasites picked up by flies
InfFlies = zeros(1,100); %#infected flies
MSPara = zeros(1,100);
P = zeros(1,100);
InfYesNo = zeros(1,100);

for ii = (1:100);

    %To find number of bites before hitting capillary
    Squares = [0,1,0,0,1,0,0,1,0,0,1,0,0,0,0,0,1,0,0,0,1,0,1,1,0,0,0,0,0,1];
        
    %attempts
    c = 1;
    
    y = datasample(Squares, 1, 'Replace', true);
    
    if y == 1; %If a capillary hit on first attempt, all well and good! %disp('Capillary') %disp(c)
            
     else y == 0; %Otherwise, it keeps sampling until it hits a capillary
            while y == 0;
                c = c + 1;
                y = datasample(Squares, 1, 'Replace', true);
            end
    end
    
    attempts(ii) = c;    
     
    MaS = (datasample(Biop,1))/1772.5;
    
    MaSD = MaS / (7*10^-6);  %Macro site density (MSD)
        
    MSPara(ii) = (MaSD * SiteArea); %Average #parasites at micro site - mean for exp. dist. that is formed for microscale patchiness
    
    %Pick an actual number of parasites for the micropatch from the exponential distribution
    P(ii) = exprnd(MSPara(ii));
            
    %Find the concentration of parasites in the blood pool formed
    Conc = P(ii)/SiteArea; 
    
    %Number of parasites picked up is the concentration x the size of the feed.
    FlyPara(ii) = round(Conc * FeedSize); 
    
    %INFECTION STATUS
    
    %Probability of being infected by the parasite burden (n) is determined by P(Inf) = 1 - P(NotInf)^n.
    NV = FlyPara(ii);
    
        if NV > 1            
            pNI = 0.01; %Probability of not being infected by 1 fly    
            PInf = 1 - ((pNI)^NV); %Probability of being infected by parasite burden            
            InfFlies(ii) = PInf;            
        else
            InfFlies(ii) = 0;
        end
        
    rn = rand;
        
    if rn < InfFlies(ii)
        InfYesNo(ii) = 1;
    else
        InfYesNo(ii) = 0;
    end        
        
end