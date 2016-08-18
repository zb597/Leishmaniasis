%CD8 mice

%% Experimental ID C Experiment 2

PreC = [0	0	0.0029	0.0247	0	1.8955	0	0.19	0	0	0.0116	0	0	173.9457	0	0	0	0	0.0143	0	0	0.0013	0.0369	0	0	0.0988	0.0002	0	0	0	0	0	0	0	0.0005	0	0	0	0	0];
PostC = [0	5.5759	0.4435	121.3444	113.6761	23.1212	0.549	0.4083	59270.81	2190.943	30.0222	3856.608	0.1187	5967.567	1402.019	0	68.2327	485.7335	5.9257	5087.115	0	0.0577	0.0361	0	124.4428	0	0	0	0	0	0	0	0.1875	0	27671.09	22774.06	0	12917.47	0.3151	2074.783];

%data
DataPreC = round(reshape(PreC, [numel(PreC), 1]));
DataPostC = round(reshape(PostC, [numel(PostC), 1]));

%log likelihood
norPreC = fitdist(DataPreC, 'Normal');
xPreC = norPreC.NLogL; %log likelihood (normal dist.)
norPostC = fitdist(DataPostC, 'Normal');
xPostC = norPostC.NLogL; %log likelihood (normal dist.)

nbPreC = fitdist(DataPreC, 'Negative Binomial');
yPreC = nbPreC.NLogL; %log likelihood (Neg. Bin. dist.)
nbPostC = fitdist(DataPostC, 'Negative Binomial');
yPostC = nbPostC.NLogL; %log likelihood (Neg. Bin. dist.)

expPreC = fitdist(DataPreC, 'Exponential');
zPreC = expPreC.NLogL; %log likelihood (exponential dist.)
expPostC = fitdist(DataPostC, 'Exponential');
zPostC = expPostC.NLogL; %log likelihood (exponential dist.)

poiPreC = fitdist(DataPreC, 'Poisson');
wPreC = poiPreC.NLogL; %log likelihood (poisson dist.)
poiPostC = fitdist(DataPostC, 'Poisson');
wPostC = poiPostC.NLogL; %log likelihood (poisson dist.)

%Akaike Information Criterion
aicxPreC = 2*2-2*log(xPreC); %AIC = 2*#parameters - 2ln(ll)
aicyPreC = 2*2-2*log(yPreC);
aiczPreC = 2*1-2*log(zPreC);
aicwPreC = 2*1-2*log(wPreC);
aicPreC = aicbic([xPreC, yPreC, zPreC, wPreC],[2, 2, 1, 1]); 
aicxPostC = 2*2-2*log(xPostC); %AIC = 2*#parameters - 2ln(ll)
aicyPostC = 2*2-2*log(yPostC);
aiczPostC = 2*1-2*log(zPostC);
aicwPostC = 2*1-2*log(wPostC);
aicPostC = aicbic([xPostC, yPostC, zPostC, wPostC],[2, 2, 1, 1]);

%Anderson-Darling test
distnPreC = makedist('normal','mu',4.4,'sigma',27.5055);
[h,p] = adtest(DataPreC,'Distribution',distnPreC);
distnPostC = makedist('normal','mu',3604.8,'sigma',10780.1);
[h,p] = adtest(DataPostC,'Distribution',distnPostC);

distnbPreC = makedist('negative binomial','r',0.00816548,'p',0.00185235);
[h,p] = adtest(DataPreC,'Distribution',distnbPreC);
distnbPostC = makedist('negative binomial','r',0.0657911,'p',1.82506e-05);
[h,p] = adtest(DataPostC,'Distribution',distnbPostC);

distePreC = makedist('exponential','mu',4.4);
[h,p] = adtest(DataPreC,'Distribution',distePreC);
distePostC = makedist('exponential','mu',3604.8);
[h,p] = adtest(DataPostC,'Distribution',distePostC);

distpPreC = makedist('poisson','lambda',4.4);
[h,p] = adtest(DataPreC,'Distribution',distpPreC);
distpPostC = makedist('poisson','lambda',3604.8);
[h,p] = adtest(DataPostC,'Distribution',distpPostC);

%fitting distributions
figure;
c1 = subplot(2,2,1);
norPreC = fitdist(DataPreC, 'Normal');
histfit(DataPreC,11,'Normal')
NormPreC = mle(DataPreC, 'distribution', 'Normal');
xlabel('Parasite Count')
ylabel('Density')
title('Normal Distribution')
axis([0,170,0,40])
suptitle('Pre Mouse C')

c2 = subplot(2,2,2);
nbPreC = fitdist(DataPreC, 'Negative Binomial');
histfit(DataPreC,11,'Negative Binomial')
NegBinPreC = mle(DataPreC, 'distribution', 'Negative Binomial');
xlabel('Parasite Count')
ylabel('Density')
title('Negative Binomial Distribution')
axis([0,170,0,40])

c3 = subplot(2,2,3);
expPreC = fitdist(DataPreC, 'Exponential');
histfit(DataPreC,11,'Exponential')
ExpoPreC = mle(DataPreC, 'distribution', 'Exponential');
xlabel('Parasite Count')
ylabel('Density')
title('Exponential Distribution')
axis([0,170,0,40])

c4 = subplot(2,2,4);
poiPreC = fitdist(DataPreC, 'Poisson');
histfit(DataPreC,11,'Poisson')
PoiPreC = mle(DataPreC, 'distribution', 'Poisson');
xlabel('Parasite Count')
ylabel('Density')
title('Poisson Distribution')
axis([0,170,0,40])

figure;
C1 = subplot(2,2,1);
norPostC = fitdist(DataPostC, 'Normal');
histfit(DataPostC,11,'Normal')
NormPostC = mle(DataPostC, 'distribution', 'Normal');
xlabel('Parasite Count')
ylabel('Density')
title('Normal Distribution')
xlim([0,inf])
suptitle('Post Mouse C')

C2 = subplot(2,2,2);
nbPostC = fitdist(DataPostC, 'Negative Binomial');
histfit(DataPostC,11,'Negative Binomial')
NegBinPostC = mle(DataPostC, 'distribution', 'Negative Binomial');
xlabel('Parasite Count')
ylabel('Density')
title('Negative Binomial Distribution')
axis([0,59400,0,40])

C3 = subplot(2,2,3);
expPostC = fitdist(DataPostC, 'Exponential');
histfit(DataPostC,11,'Exponential')
ExpoPostC = mle(DataPostC, 'distribution', 'Exponential');
xlabel('Parasite Count')
ylabel('Density')
title('Exponential Distribution')

C4 = subplot(2,2,4);
poiPostC = fitdist(DataPostC, 'Poisson');
histfit(DataPostC,11,'Poisson')
PoiPostC = mle(DataPostC, 'distribution', 'Poisson');
xlabel('Parasite Count')
ylabel('Density')
title('Poisson Distribution')

linkaxes([C1,C3,C4],'xy')

%% Experimental ID I Experiment 2

PreI = [0.4456	0	0	0	0.0082	0.2178	0.0056	0	0.8832	0.8085	1.6268	0	0.0251	0.0654	0	0	0	0	0.3765	0.0617	3.316	2.4168	1.5205	0	0.0385	0	0	0.1547	1.2092	0	0	0.0288	0	0	0	0	0.0061	22.2356	0.0042	0.1576];
PostI = [0.8487	0	0	84.0987	183.0519	0.0207	0.6195	0.3687	0.0242	0.0119	0	17.098	54.4897	0	0.2219	1.4007	0	0.1979	0	0	0.1827	0	0	1.5113	0	0.0155	0.0306	0.271	284.9258	0	0.0644	0.0081	0.1213	0.0274	0.0217	0.0707	0.0455	0	0.0098	8.3875];

%data
DataPreI = round(reshape(PreI, [numel(PreI), 1]));
DataPostI = round(reshape(PostI, [numel(PostI), 1]));

%log likelihood
norPreI = fitdist(DataPreI, 'Normal');
xPreI = norPreI.NLogL; %log likelihood (normal dist.)
norPostI = fitdist(DataPostI, 'Normal');
xPostI = norPostI.NLogL; %log likelihood (normal dist.)

nbPreI = fitdist(DataPreI, 'Negative Binomial');
yPreI = nbPreI.NLogL; %log likelihood (Neg. Bin. dist.)
nbPostI = fitdist(DataPostI, 'Negative Binomial');
yPostI = nbPostI.NLogL; %log likelihood (Neg. Bin. dist.)

expPreI = fitdist(DataPreI, 'Exponential');
zPreI = expPreI.NLogL; %log likelihood (exponential dist.)
expPostI = fitdist(DataPostI, 'Exponential');
zPostI = expPostI.NLogL; %log likelihood (exponential dist.)

poiPreI = fitdist(DataPreI, 'Poisson');
wPreI = poiPreI.NLogL; %log likelihood (poisson dist.)
poiPostI = fitdist(DataPostI, 'Poisson');
wPostI = poiPostI.NLogL; %log likelihood (poisson dist.)

%Akaike Information Criterion
aicxPreI = 2*2-2*log(xPreI); %AIC = 2*#parameters - 2ln(ll)
aicyPreI = 2*2-2*log(yPreI);
aiczPreI = 2*1-2*log(zPreI);
aicwPreI = 2*1-2*log(wPreI);
aicPreI = aicbic([xPreI, yPreI, zPreI, wPreI],[2, 2, 1, 1]); 
aicxPostI = 2*2-2*log(xPostI); %AIC = 2*#parameters - 2ln(ll)
aicyPostI = 2*2-2*log(yPostI);
aiczPostI = 2*1-2*log(zPostI);
aicwPostI = 2*1-2*log(wPostI);
aicPostI = aicbic([xPostI, yPostI, zPostI, wPostI],[2, 2, 1, 1]);

%Anderson-Darling test
distnPreI = makedist('normal','mu',0.85,'sigma',3.50494);
[h,p] = adtest(DataPreI,'Distribution',distnPreI);
distnPostI = makedist('normal','mu',15.9,'sigma',54.2864);
[h,p] = adtest(DataPostI,'Distribution',distnPostI);

distnbPreI = makedist('negative binomial','r',0.0949068,'p',0.10044);
[h,p] = adtest(DataPreI,'Distribution',distnbPreI);
distnbPostI = makedist('negative binomial','r',0.0492123,'p',0.00308557);
[h,p] = adtest(DataPostI,'Distribution',distnbPostI);

distePreB = makedist('exponential','mu',0.85);
[h,p] = adtest(DataPreB,'Distribution',distePreB);
distePostB = makedist('exponential','mu',15.9);
[h,p] = adtest(DataPostB,'Distribution',distePostB);

distpPreI = makedist('poisson','lambda',0.85);
[h,p] = adtest(DataPreI,'Distribution',distpPreI);
distpPostI = makedist('poisson','lambda',15.9);
[h,p] = adtest(DataPostI,'Distribution',distpPostI);

%fitting distributions
figure;
i1 = subplot(2,2,1);
norPreI = fitdist(DataPreI, 'Normal');
histfit(DataPreI,11,'Normal')
NormPreI = mle(DataPreI, 'distribution', 'Normal');
xlabel('Parasite Count')
ylabel('Density')
title('Normal Distribution')
xlim([0,inf])
suptitle('Pre Mouse I')

i2 = subplot(2,2,2);
nbPreI = fitdist(DataPreI, 'Negative Binomial');
histfit(DataPreI,11,'Negative Binomial')
NegBinPreI = mle(DataPreI, 'distribution', 'Negative Binomial');
xlabel('Parasite Count')
ylabel('Density')
title('Negative Binomial Distribution')
axis([0,22,0,40])

i3 = subplot(2,2,3);
expPreI = fitdist(DataPreI, 'Exponential');
histfit(DataPreI,11,'Exponential')
ExpoPreI = mle(DataPreI, 'distribution', 'Exponential');
xlabel('Parasite Count')
ylabel('Density')
title('Exponential Distribution')

i4 = subplot(2,2,4);
poiPreI = fitdist(DataPreI, 'Poisson');
histfit(DataPreI,11,'Poisson')
PoiPreI = mle(DataPreI, 'distribution', 'Poisson');
xlabel('Parasite Count')
ylabel('Density')
title('Poisson Distribution')

linkaxes([i1,i3,i4],'xy')

figure;
I1 = subplot(2,2,1);
norPostI = fitdist(DataPostI, 'Normal');
histfit(DataPostI,11,'Normal')
NormPostI = mle(DataPostI, 'distribution', 'Normal');
xlabel('Parasite Count')
ylabel('Density')
title('Normal Distribution')
xlim([0,inf])
suptitle('Post Mouse I')

I2 = subplot(2,2,2);
nbPostI = fitdist(DataPostI, 'Negative Binomial');
histfit(DataPostI,11,'Negative Binomial')
NegBinPostI = mle(DataPostI, 'distribution', 'Negative Binomial');
xlabel('Parasite Count')
ylabel('Density')
title('Negative Binomial Distribution')
axis([0,286,0,40])

I3 = subplot(2,2,3);
expPostB = fitdist(DataPostI, 'Exponential');
histfit(DataPostI,11,'Exponential')
ExpoPostI = mle(DataPostI, 'distribution', 'Exponential');
xlabel('Parasite Count')
ylabel('Density')
title('Exponential Distribution')

I4 = subplot(2,2,4);
poiPostI = fitdist(DataPostI, 'Poisson');
histfit(DataPostI,11,'Poisson')
PoiPostI = mle(DataPostI, 'distribution', 'Poisson');
xlabel('Parasite Count')
ylabel('Density')
title('Poisson Distribution')

linkaxes([I1,I3,I4],'xy')

%% Experimental ID 4 Experiment 4

Pre4 = [0.1	0.9	0	0	0	0	0	0.2	0	0	0.7	0.1	0	0	0	0.4	0	0	0	0	0	0	0	1.1	0	0	0	0	0	0	1454.3	0	0	63.4	0	0	0	0.1	0	0];
Post4 = [1.5	2.8	2.9	0.4	0.4	7528.2	3.7	471.4	2.3	1.9	0.4	0.4	206.3	2283.5	2654.9	1.7	1.6	0.1	0.4	898.2	0.4	1.1	12.4	0.6	42.4	5.7	2.5	1.9	0.3	0.5	2032.9	11.9	0.2	22.4	0.1	2.5	3.1	0.3	226.2	4.2];

%data
DataPre4 = round(reshape(Pre4, [numel(Pre4), 1]));
DataPost4 = round(reshape(Post4, [numel(Post4), 1]));

%log likelihood
norPre4 = fitdist(DataPre4, 'Normal');
xPre4 = norPre4.NLogL; %log likelihood (normal dist.)
norPost4 = fitdist(DataPost4, 'Normal');
xPost4 = norPost4.NLogL; %log likelihood (normal dist.)

nbPre4 = fitdist(DataPre4, 'Negative Binomial');
yPre4 = nbPre4.NLogL; %log likelihood (Neg. Bin. dist.)
nbPost4 = fitdist(DataPost4, 'Negative Binomial');
yPost4 = nbPost4.NLogL; %log likelihood (Neg. Bin. dist.)

expPre4 = fitdist(DataPre4, 'Exponential');
zPre4 = expPre4.NLogL; %log likelihood (exponential dist.)
expPost4 = fitdist(DataPost4, 'Exponential');
zPost4 = expPost4.NLogL; %log likelihood (exponential dist.)

poiPre4 = fitdist(DataPre4, 'Poisson');
wPre4 = poiPre4.NLogL; %log likelihood (poisson dist.)
poiPost4 = fitdist(DataPost4, 'Poisson');
wPost4 = poiPost4.NLogL; %log likelihood (poisson dist.)

%Akaike Information Criterion
aicxPre4 = 2*2-2*log(xPre4); %AIC = 2*#parameters - 2ln(ll)
aicyPre4 = 2*2-2*log(yPre4);
aiczPre4 = 2*1-2*log(zPre4);
aicwPre4 = 2*1-2*log(wPre4);
aicPre4 = aicbic([xPre4, yPre4, zPre4, wPre4],[2, 2, 1, 1]); 
aicxPost4 = 2*2-2*log(xPost4); %AIC = 2*#parameters - 2ln(ll)
aicyPost4 = 2*2-2*log(yPost4);
aiczPost4 = 2*1-2*log(zPost4);
aicwPost4 = 2*1-2*log(wPost4);
aicPost4 = aicbic([xPost4, yPost4, zPost4, wPost4],[2, 2, 1, 1]);

%Anderson-Darling test
distnPre4 = makedist('normal','mu',38,'sigma',229.845);
[h,p] = adtest(DataPre4,'Distribution',distnPre4);
distnPost4 = makedist('normal','mu',410.825,'sigma',1315.03);
[h,p] = adtest(DataPost4,'Distribution',distnPost4);

distnbPre4 = makedist('negative binomial','r',0.016876,'p',0.000443909);
[h,p] = adtest(DataPre4,'Distribution',distnbPre4);
distnbPost4 = makedist('negative binomial','r',0.123687,'p',0.00030098);
[h,p] = adtest(DataPost4,'Distribution',distnbPost4);

distePre4 = makedist('exponential','mu',38);
[h,p] = adtest(DataPre4,'Distribution',distePre4);
distePost4 = makedist('exponential','mu',410.825);
[h,p] = adtest(DataPost4,'Distribution',distePost4);

distpPre4 = makedist('poisson','lambda',38);
[h,p] = adtest(DataPre4,'Distribution',distpPre4);
distpPost4 = makedist('poisson','lambda',410.825);
[h,p] = adtest(DataPost4,'Distribution',distpPost4);

%fitting distributions
figure;
r41 = subplot(2,2,1);
norPre4 = fitdist(DataPre4, 'Normal');
histfit(DataPre4,11,'Normal')
NormPre4 = mle(DataPre4, 'distribution', 'Normal');
xlabel('Parasite Count')
ylabel('Density')
title('Normal Distribution')
xlim([0,inf])
suptitle('Pre Mouse RAG 4 Experiment 4')

r42 = subplot(2,2,2);
nbPre4 = fitdist(DataPre4, 'Negative Binomial');
histfit(DataPre4,11,'Negative Binomial')
NegBinPre4 = mle(DataPre4, 'distribution', 'Negative Binomial');
xlabel('Parasite Count')
ylabel('Density')
title('Negative Binomial Distribution')
axis([0,1540,0,40])

r43 = subplot(2,2,3);
expPre4 = fitdist(DataPre4, 'Exponential');
histfit(DataPre4,11,'Exponential')
ExpoPre4 = mle(DataPre4, 'distribution', 'Exponential');
xlabel('Parasite Count')
ylabel('Density')
title('Exponential Distribution')

r44 = subplot(2,2,4);
poiPre4 = fitdist(DataPre4, 'Poisson');
histfit(DataPre4,11,'Poisson')
PoiPre4 = mle(DataPre4, 'distribution', 'Poisson');
xlabel('Parasite Count')
ylabel('Density')
title('Poisson Distribution')

linkaxes([r41,r43,r44],'xy')

figure;
R41 = subplot(2,2,1);
norPost4 = fitdist(DataPost4, 'Normal');
histfit(DataPost4,11,'Normal')
NormPost4 = mle(DataPost4, 'distribution', 'Normal');
xlabel('Parasite Count')
ylabel('Density')
title('Normal Distribution')
xlim([0,inf])
suptitle('Post Mouse RAG 4 Experiment 4')

R42 = subplot(2,2,2);
nbPost = fitdist(DataPost4, 'Negative Binomial');
histfit(DataPost4,11,'Negative Binomial')
NegBinPost4 = mle(DataPost4, 'distribution', 'Negative Binomial');
xlabel('Parasite Count')
ylabel('Density')
title('Negative Binomial Distribution')
axis([0,7590,0,40])

R43 = subplot(2,2,3);
expPost4 = fitdist(DataPost4, 'Exponential');
histfit(DataPost4,11,'Exponential')
ExpoPost4 = mle(DataPost4, 'distribution', 'Exponential');
xlabel('Parasite Count')
ylabel('Density')
title('Exponential Distribution')

R44 = subplot(2,2,4);
poiPost4 = fitdist(DataPost4, 'Poisson');
histfit(DataPost4,11,'Poisson')
PoiPost4 = mle(DataPost4, 'distribution', 'Poisson');
xlabel('Parasite Count')
ylabel('Density')
title('Poisson Distribution')

linkaxes([R41,R43,R44],'xy')

%% Experimental ID 5 Experiment 4

Pre54 = [0.3	0.1	0	24.4	2.2	0.3	0	0	0.1	2.4	9	2.9	0	508.2	0	1232.7	0.2	0	0	0	0.1	0	0.1	0	0	42.2	0.3	0	0	9.4	1168.3	0	0	1.2	0	0.2	0.6	0.6	0.1	0.1];
Post54 = [863.9	4821.4	242.7	0.1	0	1.4	0.1	1.6	164.7	1.2	0	0.5	0	8668.3	0	0.1	43.3	273.4	0.3	2336.4	1475.9	0	0.4	0.2	19.8	2676.5	944.5	8491.8	294.2	0.1	0.5	3.7	1.3	208.4	1.3	233.1	0	4.1	5.4	0.4];

%data
DataPre54 = round(reshape(Pre54, [numel(Pre54), 1]));
DataPost54 = round(reshape(Post54, [numel(Post54), 1]));

%log likelihood
norPre54 = fitdist(DataPre54, 'Normal');
xPre54 = norPre54.NLogL; %log likelihood (normal dist.)
norPost54 = fitdist(DataPost54, 'Normal');
xPost54 = norPost54.NLogL; %log likelihood (normal dist.)

nbPre54 = fitdist(DataPre54, 'Negative Binomial');
yPre54 = nbPre54.NLogL; %log likelihood (Neg. Bin. dist.)
nbPost54 = fitdist(DataPost54, 'Negative Binomial');
yPost54 = nbPost54.NLogL; %log likelihood (Neg. Bin. dist.)

expPre54 = fitdist(DataPre54, 'Exponential');
zPre54 = expPre54.NLogL; %log likelihood (exponential dist.)
expPost54 = fitdist(DataPost54, 'Exponential');
zPost54 = expPost54.NLogL; %log likelihood (exponential dist.)

poiPre54 = fitdist(DataPre54, 'Poisson');
wPre54 = poiPre54.NLogL; %log likelihood (poisson dist.)
poiPost54 = fitdist(DataPost54, 'Poisson');
wPost54 = poiPost54.NLogL; %log likelihood (poisson dist.)

%Akaike Information Criterion
aicxPre54 = 2*2-2*log(xPre54); %AIC = 2*#parameters - 2ln(ll)
aicyPre54 = 2*2-2*log(yPre54);
aiczPre54 = 2*1-2*log(zPre54);
aicwPre54 = 2*1-2*log(wPre54);
aicPre54 = aicbic([xPre54, yPre54, zPre54, wPre54],[2, 2, 1, 1]); 
aicxPost54 = 2*2-2*log(xPost54); %AIC = 2*#parameters - 2ln(ll)
aicyPost54 = 2*2-2*log(yPost54);
aiczPost54 = 2*1-2*log(zPost54);
aicwPost54 = 2*1-2*log(wPost54);
aicPost54 = aicbic([xPost54, yPost54, zPost54, wPost54],[2, 2, 1, 1]);

%Anderson-Darling test
distnPre54 = makedist('normal','mu',75.075,'sigma',273.611);
[h,p] = adtest(DataPre54,'Distribution',distnPre54);
distnPost54 = makedist('normal','mu',794.475,'sigma',2038.29);
[h,p] = adtest(DataPost54,'Distribution',distnPost54);

distnbPre54 = makedist('negative binomial','r',0.0513579,'p',0.00068362);
[h,p] = adtest(DataPre54,'Distribution',distnbPre54);
distnbPost54 = makedist('negative binomial','r',0.107996,'p',0.000135915);
[h,p] = adtest(DataPost54,'Distribution',distnbPost54);

distePre54 = makedist('exponential','mu',75.075);
[h,p] = adtest(DataPre54,'Distribution',distePre54);
distePost54 = makedist('exponential','mu',794.475);
[h,p] = adtest(DataPost54,'Distribution',distePost54);

distpPre54 = makedist('poisson','lambda',75.075);
[h,p] = adtest(DataPre54,'Distribution',distpPre54);
distpPost54 = makedist('poisson','lambda',794.475);
[h,p] = adtest(DataPost54,'Distribution',distpPost54);

%fitting distributions
figure;
r541 = subplot(2,2,1);
norPre54 = fitdist(DataPre54, 'Normal');
histfit(DataPre54,11,'Normal')
NormPre54 = mle(DataPre54, 'distribution', 'Normal');
xlabel('Parasite Count')
ylabel('Density')
title('Normal Distribution')
xlim([0,inf])
suptitle('Pre Mouse RAG 5 Experiment 4')

r542 = subplot(2,2,2);
nbPre54 = fitdist(DataPre54, 'Negative Binomial');
histfit(DataPre54,11,'Negative Binomial')
NegBinPre54 = mle(DataPre54, 'distribution', 'Negative Binomial');
xlabel('Parasite Count')
ylabel('Density')
title('Negative Binomial Distribution')
axis([0,1320,0,40])

r543 = subplot(2,2,3);
expPre54 = fitdist(DataPre54, 'Exponential');
histfit(DataPre54,11,'Exponential')
ExpoPre54 = mle(DataPre54, 'distribution', 'Exponential');
xlabel('Parasite Count')
ylabel('Density')
title('Exponential Distribution')

r544 = subplot(2,2,4);
poiPre54 = fitdist(DataPre54, 'Poisson');
histfit(DataPre54,11,'Poisson')
PoiPre54 = mle(DataPre54, 'distribution', 'Poisson');
xlabel('Parasite Count')
ylabel('Density')
title('Poisson Distribution')

linkaxes([r541,r543,r544],'xy')

figure;
R541 = subplot(2,2,1);
norPost54 = fitdist(DataPost54, 'Normal');
histfit(DataPost54,11,'Normal')
NormPost54 = mle(DataPost54, 'distribution', 'Normal');
xlabel('Parasite Count')
ylabel('Density')
title('Normal Distribution')
xlim([0,inf])
suptitle('Post Mouse RAG 5 Experiment 4')

R542 = subplot(2,2,2);
nbPost54 = fitdist(DataPost54, 'Negative Binomial');
histfit(DataPost54,11,'Negative Binomial')
NegBinPost54 = mle(DataPost54, 'distribution', 'Negative Binomial');
xlabel('Parasite Count')
ylabel('Density')
title('Negative Binomial Distribution')
axis([0,8690,0,40])

R543 = subplot(2,2,3);
expPost54 = fitdist(DataPost54, 'Exponential');
histfit(DataPost54,11,'Exponential')
ExpoPost54 = mle(DataPost54, 'distribution', 'Exponential');
xlabel('Parasite Count')
ylabel('Density')
title('Exponential Distribution')

R544 = subplot(2,2,4);
poiPost54 = fitdist(DataPost54, 'Poisson');
histfit(DataPost54,11,'Poisson')
PoiPost54 = mle(DataPost54, 'distribution', 'Poisson');
xlabel('Parasite Count')
ylabel('Density')
title('Poisson Distribution')

linkaxes([R541,R543,R544],'xy')

%% Experimental ID 7 Experiment 4

Pre7 = [60894.7	29312.6	149.9	329.2	0.8	828.8	1.1	7.2	23167.8	46533.2	716.8	10084.5	12668.4	2778	23218	20898.7	1.3	2162.7	5618	31991.6	2385	7459.2	24130.7	12464.5	9145.1	78655	25509	13988.6	23600.3	3259.9	59350.7	20476.9	3.3	3629.2	198.5	2430.9	3.8	114.8	27108.9	0.3];
Post7 = [12.5	49844.3	925.7	51171.1	8598	1.4	2846.3	131405.3	51407.1	33692.4	18234.4	327061.6	27193.8	44500.2	390562.3	44369.2	68823.3	15961.9	52160.3	7226.3	3	40.4	10832.6	9240.8	3427.2	27649.3	24212.2	154125.4	20413.9	29544.5	4001.6	24486.1	10991.4	48774.8	28130.8	416.6	55300.7	37953.5	51720.2	174.2];

%data
DataPre7 = round(reshape(Pre7, [numel(Pre7), 1]));
DataPost7= round(reshape(Post7, [numel(Post7), 1]));

%log likelihood
norPre7 = fitdist(DataPre7, 'Normal');
xPre7 = norPre7.NLogL; %log likelihood (normal dist.)
norPost7 = fitdist(DataPost7, 'Normal');
xPost7 = norPost7.NLogL; %log likelihood (normal dist.)

nbPre7 = fitdist(DataPre7, 'Negative Binomial');
yPre7 = nbPre7.NLogL; %log likelihood (Neg. Bin. dist.)
nbPost7 = fitdist(DataPost7, 'Negative Binomial');
yPost7 = nbPost7.NLogL; %log likelihood (Neg. Bin. dist.)

expPre7 = fitdist(DataPre7, 'Exponential');
zPre7 = expPre7.NLogL; %log likelihood (exponential dist.)
expPost7 = fitdist(DataPost7, 'Exponential');
zPost7 = expPost7.NLogL; %log likelihood (exponential dist.)

poiPre7 = fitdist(DataPre7, 'Poisson');
wPre7 = poiPre7.NLogL; %log likelihood (poisson dist.)
poiPost7 = fitdist(DataPost7, 'Poisson');
wPost7 = poiPost7.NLogL; %log likelihood (poisson dist.)

%Akaike Information Criterion
aicxPre7 = 2*2-2*log(xPre7); %AIC = 2*#parameters - 2ln(ll)
aicyPre7 = 2*2-2*log(yPre7);
aiczPre7 = 2*1-2*log(zPre7);
aicwPre7 = 2*1-2*log(wPre7);
aicPre7 = aicbic([xPre7, yPre7, zPre7, wPre7],[2, 2, 1, 1]); 
aicxPost7 = 2*2-2*log(xPost7); %AIC = 2*#parameters - 2ln(ll)
aicyPost7 = 2*2-2*log(yPost7);
aiczPost7 = 2*1-2*log(zPost7);
aicwPost7 = 2*1-2*log(wPost7);
aicPost7 = aicbic([xPost7, yPost7, zPost7, wPost7],[2, 2, 1, 1]);

%Anderson-Darling test
distnPre7 = makedist('normal','mu',14632,'sigma',19075.8);
[h,p] = adtest(DataPre7,'Distribution',distnPre7);
distnPost7 = makedist('normal','mu',46685.9,'sigma',79908.8);
[h,p] = adtest(DataPost7,'Distribution',distnPost7);

distnbPre7 = makedist('negative binomial','r',0.299062,'p',0.0000204384);
[h,p] = adtest(DataPre7,'Distribution',distnbPre7);
distnbPost7 = makedist('negative binomial','r',0.402308,'p',0.00000861726);
[h,p] = adtest(DataPost7,'Distribution',distnbPost7);

distePre7 = makedist('exponential','mu',14632);
[h,p] = adtest(DataPre7,'Distribution',distePre7);
distePost7 = makedist('exponential','mu',46685.9);
[h,p] = adtest(DataPost7,'Distribution',distePost7);

distpPre7 = makedist('poisson','lambda',14632);
[h,p] = adtest(DataPre7,'Distribution',distpPre7);
distpPost7 = makedist('poisson','lambda',46685.9   );
[h,p] = adtest(DataPost7,'Distribution',distpPost7);

%fitting distributions
figure;
r71 = subplot(2,2,1);
norPre7 = fitdist(DataPre7, 'Normal');
histfit(DataPre7,11,'Normal')
NormPre7 = mle(DataPre7, 'distribution', 'Normal');
xlabel('Parasite Count')
ylabel('Density')
title('Normal Distribution')
xlim([0,inf])
suptitle('Pre Mouse RAG 7 Experiment 4')

r72 = subplot(2,2,2);
nbPre7 = fitdist(DataPre7, 'Negative Binomial');
histfit(DataPre7,11,'Negative Binomial')
NegBinPre7 = mle(DataPre7, 'distribution', 'Negative Binomial');
xlabel('Parasite Count')
ylabel('Density')
title('Negative Binomial Distribution')
axis([0,79200,0,20])

r73 = subplot(2,2,3);
expPre7 = fitdist(DataPre7, 'Exponential');
histfit(DataPre7,11,'Exponential')
ExpoPre7 = mle(DataPre7, 'distribution', 'Exponential');
xlabel('Parasite Count')
ylabel('Density')
title('Exponential Distribution')

r74 = subplot(2,2,4);
poiPre7 = fitdist(DataPre7, 'Poisson');
histfit(DataPre7,11,'Poisson')
PoiPre7 = mle(DataPre7, 'distribution', 'Poisson');
xlabel('Parasite Count')
ylabel('Density')
title('Poisson Distribution')

linkaxes([r71,r73,r74],'xy')

figure;
R71 = subplot(2,2,1);
norPost7 = fitdist(DataPost7, 'Normal');
histfit(DataPost7,11,'Normal')
NormPost7 = mle(DataPost7, 'distribution', 'Normal');
xlabel('Parasite Count')
ylabel('Density')
title('Normal Distribution')
xlim([0,inf])
suptitle('Post Mouse RAG 7 Experiment 4')

R72 = subplot(2,2,2);
nbPost7 = fitdist(DataPost7, 'Negative Binomial');
histfit(DataPost7,11,'Negative Binomial')
NegBinPost7 = mle(DataPost7, 'distribution', 'Negative Binomial');
xlabel('Parasite Count')
ylabel('Density')
title('Negative Binomial Distribution')
axis([0,396000,0,30])

R73 = subplot(2,2,3);
expPost7 = fitdist(DataPost7, 'Exponential');
histfit(DataPost7,11,'Exponential')
ExpoPost7 = mle(DataPost7, 'distribution', 'Exponential');
xlabel('Parasite Count')
ylabel('Density')
title('Exponential Distribution')

R74 = subplot(2,2,4);
poiPost7 = fitdist(DataPost7, 'Poisson');
histfit(DataPost7,11,'Poisson')
PoiPost7 = mle(DataPost7, 'distribution', 'Poisson');
xlabel('Parasite Count')
ylabel('Density')
title('Poisson Distribution')

linkaxes([R71,R73,R74],'xy')

%% Experimental ID 9 Experiment 4

Pre9 = [0.4	0.1	0.1	0	0.4	0.1	136.7	0.1	0.6	0.1	0	1.5	0.6	0	2	0	1.4	0.1	0	1.7	0	0.6	1.3	0	0.3	0.1	0.2	0.9	0	0	0.3	0.1	1.2	1.1	0	2079.3	0.7	0.7	0.2	0.1];
Post9 = [191.1	4.3	271.7	3	624	0	782.5	0.1	0.1	0.5	4	0.2	0	0	0	0	0	0.9	0.7	4.9	0.5	0.4	93	0.7	0.7	47.9	155.6	0.4	1.1	1.5	0	0	0.1	0.6	0	0.1	0.1	0.5	1.7	4720.3];

%data
DataPre9 = round(reshape(Pre9, [numel(Pre9), 1]));
DataPost9 = round(reshape(Post9, [numel(Post9), 1]));

%log likelihood
norPre9 = fitdist(DataPre9, 'Normal');
xPre9 = norPre9.NLogL; %log likelihood (normal dist.)
norPost9 = fitdist(DataPost9, 'Normal');
xPost9 = norPost9.NLogL; %log likelihood (normal dist.)

nbPre9 = fitdist(DataPre9, 'Negative Binomial');
yPre9 = nbPre9.NLogL; %log likelihood (Neg. Bin. dist.)
nbPost9 = fitdist(DataPost9, 'Negative Binomial');
yPost9 = nbPost9.NLogL; %log likelihood (Neg. Bin. dist.)

expPre9 = fitdist(DataPre9, 'Exponential');
zPre9 = expPre9.NLogL; %log likelihood (exponential dist.)
expPost9 = fitdist(DataPost9, 'Exponential');
zPost9 = expPost9.NLogL; %log likelihood (exponential dist.)

poiPre9 = fitdist(DataPre9, 'Poisson');
wPre9 = poiPre9.NLogL; %log likelihood (poisson dist.)
poiPost9 = fitdist(DataPost9, 'Poisson');
wPost9 = poiPost9.NLogL; %log likelihood (poisson dist.)

%Akaike Information Criterion
aicxPre9 = 2*2-2*log(xPre9); %AIC = 2*#parameters - 2ln(ll)
aicyPre9 = 2*2-2*log(yPre9);
aiczPre9 = 2*1-2*log(zPre9);
aicwPre9 = 2*1-2*log(wPre9);
aicPre9 = aicbic([xPre9, yPre9, zPre9, wPre9],[2, 2, 1, 1]); 
aicxPost9 = 2*2-2*log(xPost9); %AIC = 2*#parameters - 2ln(ll)
aicyPost9 = 2*2-2*log(yPost9);
aiczPost9 = 2*1-2*log(zPost9);
aicwPost9 = 2*1-2*log(wPost9);
aicPost9 = aicbic([xPost9, yPost9, zPost9, wPost9],[2, 2, 1, 1]);

%Anderson-Darling test
distnPre9 = makedist('normal','mu',55.8,'sigma',328.809);
[h,p] = adtest(DataPre9,'Distribution',distnPre9);
distnPost9 = makedist('normal','mu',172.9,'sigma',754.999);
[h,p] = adtest(DataPost9,'Distribution',distnPost9);

distnbPre9 = makedist('negative binomial','r',0.0580477,'p',0.0010392);
[h,p] = adtest(DataPre9,'Distribution',distnbPre9);
distnbPost9 = makedist('negative binomial','r',0.0948693,'p',0.000548394);
[h,p] = adtest(DataPost9,'Distribution',distnbPost9);

distePre9 = makedist('exponential','mu',55.8);
[h,p] = adtest(DataPre9,'Distribution',distePre9);
distePost9 = makedist('exponential','mu',172.9);
[h,p] = adtest(DataPost9,'Distribution',distePost9);

distpPre9 = makedist('poisson','lambda',55.8);
[h,p] = adtest(DataPre9,'Distribution',distpPre9);
distpPost9 = makedist('poisson','lambda',172.9);
[h,p] = adtest(DataPost9,'Distribution',distpPost9);

%fitting distributions
figure;
r91 = subplot(2,2,1);
norPre9 = fitdist(DataPre9, 'Normal');
histfit(DataPre9,11,'Normal')
NormPre9 = mle(DataPre9, 'distribution', 'Normal');
xlabel('Parasite Count')
ylabel('Density')
title('Normal Distribution')
xlim([0,inf])
suptitle('Pre Mouse RAG 9 Experiment 4')

r92 = subplot(2,2,2);
nbPre9 = fitdist(DataPre9, 'Negative Binomial');
histfit(DataPre9,11,'Negative Binomial')
NegBinPre9 = mle(DataPre9, 'distribution', 'Negative Binomial');
xlabel('Parasite Count')
ylabel('Density')
title('Negative Binomial Distribution')
axis([0,2090,0,40])

r93 = subplot(2,2,3);
expPre9 = fitdist(DataPre9, 'Exponential');
histfit(DataPre9,11,'Exponential')
ExpoPre9 = mle(DataPre9, 'distribution', 'Exponential');
xlabel('Parasite Count')
ylabel('Density')
title('Exponential Distribution')

r94 = subplot(2,2,4);
poiPre9 = fitdist(DataPre9, 'Poisson');
histfit(DataPre9,11,'Poisson')
PoiPre9 = mle(DataPre9, 'distribution', 'Poisson');
xlabel('Parasite Count')
ylabel('Density')
title('Poisson Distribution')

linkaxes([r91,r93,r94],'xy')

figure;
R91 = subplot(2,2,1);
norPost9 = fitdist(DataPost9, 'Normal');
histfit(DataPost9,11,'Normal')
NormPost9 = mle(DataPost9, 'distribution', 'Normal');
xlabel('Parasite Count')
ylabel('Density')
title('Normal Distribution')
xlim([0,inf])
suptitle('Post Mouse RAG 9 Experiment 4')

R92 = subplot(2,2,2);
nbPost9 = fitdist(DataPost9, 'Negative Binomial');
histfit(DataPost9,11,'Negative Binomial')
NegBinPost9 = mle(DataPost9, 'distribution', 'Negative Binomial');
xlabel('Parasite Count')
ylabel('Density')
title('Negative Binomial Distribution')
axis([0,4730,0,40])

R93 = subplot(2,2,3);
expPost9 = fitdist(DataPost9, 'Exponential');
histfit(DataPost9,11,'Exponential')
ExpoPost9 = mle(DataPost9, 'distribution', 'Exponential');
xlabel('Parasite Count')
ylabel('Density')
title('Exponential Distribution')

R94 = subplot(2,2,4);
poiPost9 = fitdist(DataPost9, 'Poisson');
histfit(DataPost9,11,'Poisson')
PoiPost9 = mle(DataPost9, 'distribution', 'Poisson');
xlabel('Parasite Count')
ylabel('Density')
title('Poisson Distribution')

linkaxes([R91,R93,R94],'xy')

%% Experimental ID 2 Experiment 5

Pre2 = [20.5	0	0	0	0	0	0	0	0	0	0	0	7	0	0	0	91.4	0	2020.9	0	0	0	0	0	0	0	0	0	421	0	0	1237.4	0	0	0	0	35.9	0	0	1182.8];
Post2 = [0.4	13.4	0	28.8	901.9	129.7	98.1	9.1	0.3	0.1	1.8	0	2.8	3667	329.5	3.2	0.2	0.4	0.3	0.6	0.4	2.3	228.6	3	4.1	0	0	22.3	0.2	8.8	5.7	0.2	0.2	0	3.6	1.9	0.2	130.9	1.1	0.3];

%data
DataPre2 = round(reshape(Pre2, [numel(Pre2), 1]));
DataPost2 = round(reshape(Post2, [numel(Post2), 1]));

%log likelihood
norPre2 = fitdist(DataPre2, 'Normal');
xPre2 = norPre2.NLogL; %log likelihood (normal dist.)
norPost2 = fitdist(DataPost2, 'Normal');
xPost2 = norPost2.NLogL; %log likelihood (normal dist.)

nbPre2 = fitdist(DataPre2, 'Negative Binomial');
yPre2 = nbPre2.NLogL; %log likelihood (Neg. Bin. dist.)
nbPost2 = fitdist(DataPost2, 'Negative Binomial');
yPost2 = nbPost2.NLogL; %log likelihood (Neg. Bin. dist.)

expPre2 = fitdist(DataPre2, 'Exponential');
zPre2 = expPre2.NLogL; %log likelihood (exponential dist.)
expPost2 = fitdist(DataPost2, 'Exponential');
zPost2 = expPost2.NLogL; %log likelihood (exponential dist.)

poiPre2 = fitdist(DataPre2, 'Poisson');
wPre2 = poiPre2.NLogL; %log likelihood (poisson dist.)
poiPost2 = fitdist(DataPost2, 'Poisson');
wPost2 = poiPost2.NLogL; %log likelihood (poisson dist.)

%Akaike Information Criterion
aicxPre2 = 2*2-2*log(xPre2); %AIC = 2*#parameters - 2ln(ll)
aicyPre2 = 2*2-2*log(yPre2);
aiczPre2 = 2*1-2*log(zPre2);
aicwPre2 = 2*1-2*log(wPre2);
aicPre2 = aicbic([xPre2, yPre2, zPre2, wPre2],[2, 2, 1, 1]); 
aicxPost2 = 2*2-2*log(xPost2); %AIC = 2*#parameters - 2ln(ll)
aicyPost2 = 2*2-2*log(yPost2);
aiczPost2 = 2*1-2*log(zPost2);
aicwPost2 = 2*1-2*log(wPost2);
aicPost2 = aicbic([xPost2, yPost2, zPost2, wPost2],[2, 2, 1, 1]);

%Anderson-Darling test
distnPre2 = makedist('normal','mu',125.425,'sigma',410.508);
[h,p] = adtest(DataPre2,'Distribution',distnPre2);
distnPost2 = makedist('normal','mu',140,'sigma',592.328);
[h,p] = adtest(DataPost2,'Distribution',distnPost2);

distnbPre2 = makedist('negative binomial','r',0.0273885,'p',0.000218318);
[h,p] = adtest(DataPre2,'Distribution',distnbPre2);
distnbPost2 = makedist('negative binomial','r',0.106422,'p',0.000759578);
[h,p] = adtest(DataPost2,'Distribution',distnbPost2);

distePre2 = makedist('exponential','mu',125.425);
[h,p] = adtest(DataPre2,'Distribution',distePre2);
distePost2 = makedist('exponential','mu',140);
[h,p] = adtest(DataPost2,'Distribution',distePost2);

distpPre2 = makedist('poisson','lambda',125.425);
[h,p] = adtest(DataPre2,'Distribution',distpPre2);
distpPost2 = makedist('poisson','lambda',140);
[h,p] = adtest(DataPost2,'Distribution',distpPost2);

%fitting distributions
figure;
r21 = subplot(2,2,1);
norPre2 = fitdist(DataPre2, 'Normal');
histfit(DataPre2,11,'Normal')
NormPre2 = mle(DataPre2, 'distribution', 'Normal');
xlabel('Parasite Count')
ylabel('Density')
title('Normal Distribution')
axis([0,2090,0,40])
suptitle('Pre Mouse RAG 2 Experiment 5')

r22 = subplot(2,2,2);
nbPre2 = fitdist(DataPre2, 'Negative Binomial');
histfit(DataPre2,11,'Negative Binomial')
NegBinPre2 = mle(DataPre2, 'distribution', 'Negative Binomial');
xlabel('Parasite Count')
ylabel('Density')
title('Negative Binomial Distribution')
axis([0,2090,0,40])

r23 = subplot(2,2,3);
expPre2 = fitdist(DataPre2, 'Exponential');
histfit(DataPre2,11,'Exponential')
ExpoPre2 = mle(DataPre2, 'distribution', 'Exponential');
xlabel('Parasite Count')
ylabel('Density')
title('Exponential Distribution')
axis([0,2090,0,40])

r24 = subplot(2,2,4);
poiPre2 = fitdist(DataPre2, 'Poisson');
histfit(DataPre2,11,'Poisson')
PoiPre2 = mle(DataPre2, 'distribution', 'Poisson');
xlabel('Parasite Count')
ylabel('Density')
title('Poisson Distribution')
axis([0,2090,0,40])

linkaxes([r21,r23,r24],'xy')

% figure;
R21 = subplot(2,2,1);
norPost2 = fitdist(DataPost2, 'Normal');
histfit(DataPost2,11,'Normal')
NormPost2 = mle(DataPost2, 'distribution', 'Normal');
xlabel('Parasite Count')
ylabel('Density')
title('Normal Distribution')
xlim([0,inf])
suptitle('Post Mouse RAG 2 Experiment 5')

R22 = subplot(2,2,2);
nbPost2 = fitdist(DataPost2, 'Negative Binomial');
histfit(DataPost2,11,'Negative Binomial')
NegBinPost2 = mle(DataPost2, 'distribution', 'Negative Binomial');
xlabel('Parasite Count')
ylabel('Density')
title('Negative Binomial Distribution')
axis([0,3740,0,40])

R23 = subplot(2,2,3);
expPost2 = fitdist(DataPost2, 'Exponential');
histfit(DataPost2,11,'Exponential')
ExpoPost2 = mle(DataPost2, 'distribution', 'Exponential');
xlabel('Parasite Count')
ylabel('Density')
title('Exponential Distribution')

R24 = subplot(2,2,4);
poiPost2 = fitdist(DataPost2, 'Poisson');
histfit(DataPost2,11,'Poisson')
PoiPost2 = mle(DataPost2, 'distribution', 'Poisson');
xlabel('Parasite Count')
ylabel('Density')
title('Poisson Distribution')

linkaxes([R21,R23,R24],'xy')

%% Experimental ID 5 Experiment 5

Pre55 = [10.7	32.2	24.3	8.8	18.3	15.7	9.8	89.6	21.1	7.7	35.7	12	12.2	14.1	6.2	16.7	0.7	13.7	6.6	1.4	1	9.3	1.1	3	1.3	21	11.2	0.3	8.4	0.8	7.6	1.1	3.9	1	3.4	2.3	4	0.5	0.4	25.5];
Post55 = [0.3	0.3	1.2	0	0.1	0	8.9	1.5	0.2	1.3	0.1	0.2	0	4.6	0.3	0.1	0.2	0.7	0.2	0.1	0.4	3.1	1.1	1.9	1.5	0.5	0.2	0	0.1	1.5	0.3	0.3	717.9	0.7	0.1	0.1	1.7	0.4	0.1	1];

%data
DataPre55 = round(reshape(Pre55, [numel(Pre55), 1]));
DataPost55 = round(reshape(Post55, [numel(Post55), 1]));

%log likelihood
norPre55 = fitdist(DataPre55, 'Normal');
xPre55 = norPre55.NLogL; %log likelihood (normal dist.)
norPost55 = fitdist(DataPost55, 'Normal');
xPost55 = norPost55.NLogL; %log likelihood (normal dist.)

nbPre55 = fitdist(DataPre55, 'Negative Binomial');
yPre55 = nbPre55.NLogL; %log likelihood (Neg. Bin. dist.)
nbPost55 = fitdist(DataPost55, 'Negative Binomial');
yPost55 = nbPost55.NLogL; %log likelihood (Neg. Bin. dist.)

expPre55 = fitdist(DataPre55, 'Exponential');
zPre55 = expPre55.NLogL; %log likelihood (exponential dist.)
expPost55 = fitdist(DataPost55, 'Exponential');
zPost55 = expPost55.NLogL; %log likelihood (exponential dist.)

poiPre55 = fitdist(DataPre55, 'Poisson');
wPre55 = poiPre55.NLogL; %log likelihood (poisson dist.)
poiPost55 = fitdist(DataPost55, 'Poisson');
wPost55 = poiPost55.NLogL; %log likelihood (poisson dist.)

%Akaike Information Criterion
aicxPre55 = 2*2-2*log(xPre55); %AIC = 2*#parameters - 2ln(ll)
aicyPre55 = 2*2-2*log(yPre55);
aiczPre55 = 2*1-2*log(zPre55);
aicwPre55 = 2*1-2*log(wPre55);
aicPre55 = aicbic([xPre55, yPre55, zPre55, wPre55],[2, 2, 1, 1]); 
aicxPost55 = 2*2-2*log(xPost55); %AIC = 2*#parameters - 2ln(ll)
aicyPost55 = 2*2-2*log(yPost55);
aiczPost55 = 2*1-2*log(zPost55);
aicwPost55 = 2*1-2*log(wPost55);
aicPost55 = aicbic([xPost55, yPost55, zPost55, wPost55],[2, 2, 1, 1]);

%Anderson-Darling test
distnPre55 = makedist('normal','mu',11.625,'sigma',15.6413);
[h,p] = adtest(DataPre55,'Distribution',distnPre55);
distnPost55 = makedist('normal','mu',18.8,'sigma',113.401);
[h,p] = adtest(DataPost55,'Distribution',distnPost55);

distnbPre55 = makedist('negative binomial','r',0.818208,'p',0.0657554);
[h,p] = adtest(DataPre55,'Distribution',distnbPre55);
distnbPost55 = makedist('negative binomial','r',0.079551,'p',0.00421361);
[h,p] = adtest(DataPost55,'Distribution',distnbPost55);

distePre55 = makedist('exponential','mu',11.625);
[h,p] = adtest(DataPre55,'Distribution',distePre55);
distePost55 = makedist('exponential','mu',18.8);
[h,p] = adtest(DataPost55,'Distribution',distePost55);

distpPre55 = makedist('poisson','lambda',11.625);
[h,p] = adtest(DataPre55,'Distribution',distpPre55);
distpPost55 = makedist('poisson','lambda',18.8);
[h,p] = adtest(DataPost55,'Distribution',distpPost55);

%fitting distributions
figure;
r551 = subplot(2,2,1);
norPre55 = fitdist(DataPre55, 'Normal');
histfit(DataPre55,11,'Normal')
NormPre55 = mle(DataPre55, 'distribution', 'Normal');
xlabel('Parasite Count')
ylabel('Density')
title('Normal Distribution')
xlim([0,inf])
suptitle('Pre Mouse RAG 5 Experiment 5')

r552 = subplot(2,2,2);
nbPre55 = fitdist(DataPre55, 'Negative Binomial');
histfit(DataPre55,11,'Negative Binomial')
NegBinPre55 = mle(DataPre55, 'distribution', 'Negative Binomial');
xlabel('Parasite Count')
ylabel('Density')
title('Negative Binomial Distribution')
axis([0,90,0,30])

r553 = subplot(2,2,3);
expPre55 = fitdist(DataPre55, 'Exponential');
histfit(DataPre55,11,'Exponential')
ExpoPre55 = mle(DataPre55, 'distribution', 'Exponential');
xlabel('Parasite Count')
ylabel('Density')
title('Exponential Distribution')

r554 = subplot(2,2,4);
poiPre55 = fitdist(DataPre55, 'Poisson');
histfit(DataPre55,11,'Poisson')
PoiPre55 = mle(DataPre55, 'distribution', 'Poisson');
xlabel('Parasite Count')
ylabel('Density')
title('Poisson Distribution')

linkaxes([r551,r553,r554],'xy')

figure;
R551 = subplot(2,2,1);
norPost55 = fitdist(DataPost55, 'Normal');
histfit(DataPost55,11,'Normal')
NormPost55 = mle(DataPost55, 'distribution', 'Normal');
xlabel('Parasite Count')
ylabel('Density')
title('Normal Distribution')
xlim([0,inf])
suptitle('Post Mouse RAG 5 Experiment 5')

R552 = subplot(2,2,2);
nbPost55 = fitdist(DataPost55, 'Negative Binomial');
histfit(DataPost55,11,'Negative Binomial')
NegBinPost55 = mle(DataPost55, 'distribution', 'Negative Binomial');
xlabel('Parasite Count')
ylabel('Density')
title('Negative Binomial Distribution')
axis([0,726,0,40])

R553 = subplot(2,2,3);
expPost55 = fitdist(DataPost55, 'Exponential');
histfit(DataPost55,11,'Exponential')
ExpoPost55 = mle(DataPost55, 'distribution', 'Exponential');
xlabel('Parasite Count')
ylabel('Density')
title('Exponential Distribution')

R554 = subplot(2,2,4);
poiPost55 = fitdist(DataPost55, 'Poisson');
histfit(DataPost55,11,'Poisson')
PoiPost55 = mle(DataPost55, 'distribution', 'Poisson');
xlabel('Parasite Count')
ylabel('Density')
title('Poisson Distribution')

linkaxes([R551,R553,R554],'xy')

%% Experimental ID 8 Experiment 5

Pre8 = [0	0	6	0	0	0	0	0	0.4	0	1	0	0	0.3	0	17.8	19	16.9	51.4	17.6	18.7	45.7	49.4	34.5	6.9	24.5	76.6	10.2	5.3	14.2	7.2	5.3	13.8	13	12.6	26	17.1	11.3	7.1	581.5];
Post8 = [0	0.9	0	0	0	0.4	0	0.3	2.9	0.6	0	0	0	0	0	0.1	0	0	0.3	0.1	0.1	1.4	0	0.1	0.1	0.9	0.4	0	0	0	0	0.1	0	0.2	0.8	0.1	0	8.4	0.2	0.2];

%data
DataPre8 = round(reshape(Pre8, [numel(Pre8), 1]));
DataPost8 = round(reshape(Post8, [numel(Post8), 1]));

%log likelihood
norPre8 = fitdist(DataPre8, 'Normal');
xPre8 = norPre8.NLogL; %log likelihood (normal dist.)
norPost8 = fitdist(DataPost8, 'Normal');
xPost8 = norPost8.NLogL; %log likelihood (normal dist.)

nbPre8 = fitdist(DataPre8, 'Negative Binomial');
yPre8 = nbPre8.NLogL; %log likelihood (Neg. Bin. dist.)
nbPost8 = fitdist(DataPost8, 'Negative Binomial');
yPost8 = nbPost8.NLogL; %log likelihood (Neg. Bin. dist.)

expPre8 = fitdist(DataPre8, 'Exponential');
zPre8 = expPre8.NLogL; %log likelihood (exponential dist.)
expPost8 = fitdist(DataPost8, 'Exponential');
zPost8 = expPost8.NLogL; %log likelihood (exponential dist.)

poiPre8 = fitdist(DataPre8, 'Poisson');
wPre8 = poiPre8.NLogL; %log likelihood (poisson dist.)
poiPost8 = fitdist(DataPost8, 'Poisson');
wPost8 = poiPost8.NLogL; %log likelihood (poisson dist.)

%Akaike Information Criterion
aicxPre8 = 2*2-2*log(xPre8); %AIC = 2*#parameters - 2ln(ll)
aicyPre8 = 2*2-2*log(yPre8);
aiczPre8 = 2*1-2*log(zPre8);
aicwPre8 = 2*1-2*log(wPre8);
aicPre8 = aicbic([xPre8, yPre8, zPre8, wPre8],[2, 2, 1, 1]); 
aicxPost8 = 2*2-2*log(xPost8); %AIC = 2*#parameters - 2ln(ll)
aicyPost8 = 2*2-2*log(yPost8);
aiczPost8 = 2*1-2*log(zPost8);
aicwPost8 = 2*1-2*log(wPost8);
aicPost8 = aicbic([xPost8, yPost8, zPost8, wPost8],[2, 2, 1, 1]);

%Anderson-Darling test
distnPre8 = makedist('normal','mu',27.8,'sigma',91.4998);
[h,p] = adtest(DataPre8,'Distribution',distnPre8);
distnPost8 = makedist('normal','mu',0.4,'sigma',1.35495);
[h,p] = adtest(DataPost8,'Distribution',distnPost8);

distnbPre8 = makedist('negative binomial','r',0.250973,'p',0.00894704);
[h,p] = adtest(DataPre8,'Distribution',distnbPre8);
distnbPost8 = makedist('negative binomial','r',0.138334,'p',0.256967);
[h,p] = adtest(DataPost8,'Distribution',distnbPost8);

distePre8 = makedist('exponential','mu',27.8);
[h,p] = adtest(DataPre8,'Distribution',distePre8);
distePost8 = makedist('exponential','mu',0.4);
[h,p] = adtest(DataPost8,'Distribution',distePost8);

distpPre8 = makedist('poisson','lambda',27.8);
[h,p] = adtest(DataPre8,'Distribution',distpPre8);
distpPost8 = makedist('poisson','lambda',0.4);
[h,p] = adtest(DataPost8,'Distribution',distpPost8);

%fitting distributions
figure;
r81 = subplot(2,2,1);
norPre8 = fitdist(DataPre8, 'Normal');
histfit(DataPre8,11,'Normal')
NormPre8 = mle(DataPre8, 'distribution', 'Normal');
xlabel('Parasite Count')
ylabel('Density')
title('Normal Distribution')
xlim([0,inf])
suptitle('Pre Mouse RAG 8 Experiment 5')

r82 = subplot(2,2,2);
nbPre8 = fitdist(DataPre8, 'Negative Binomial');
histfit(DataPre8,11,'Negative Binomial')
NegBinPre8 = mle(DataPre8, 'distribution', 'Negative Binomial');
xlabel('Parasite Count')
ylabel('Density')
title('Negative Binomial Distribution')
axis([0,575,0,40])

r83 = subplot(2,2,3);
expPre8 = fitdist(DataPre8, 'Exponential');
histfit(DataPre8,11,'Exponential')
ExpoPre8 = mle(DataPre8, 'distribution', 'Exponential');
xlabel('Parasite Count')
ylabel('Density')
title('Exponential Distribution')

r84 = subplot(2,2,4);
poiPre8 = fitdist(DataPre8, 'Poisson');
histfit(DataPre8,11,'Poisson')
PoiPre8 = mle(DataPre8, 'distribution', 'Poisson');
xlabel('Parasite Count')
ylabel('Density')
title('Poisson Distribution')

linkaxes([r81,r83,r84],'xy')

figure;
R81 = subplot(2,2,1);
norPost8 = fitdist(DataPost8, 'Normal');
histfit(DataPost8,11,'Normal')
NormPost8 = mle(DataPost8, 'distribution', 'Normal');
xlabel('Parasite Count')
ylabel('Density')
title('Normal Distribution')
xlim([0,inf])
suptitle('Post Mouse RAG 8 Experiment 5')

R82 = subplot(2,2,2);
nbPost8 = fitdist(DataPost8, 'Negative Binomial');
histfit(DataPost8,11,'Negative Binomial')
NegBinPost8 = mle(DataPost8, 'distribution', 'Negative Binomial');
xlabel('Parasite Count')
ylabel('Density')
title('Negative Binomial Distribution')
axis([0,8,0,40])

R83 = subplot(2,2,3);
expPost8 = fitdist(DataPost8, 'Exponential');
histfit(DataPost8,11,'Exponential')
ExpoPost8 = mle(DataPost8, 'distribution', 'Exponential');
xlabel('Parasite Count')
ylabel('Density')
title('Exponential Distribution')

R84 = subplot(2,2,4);
poiPost8 = fitdist(DataPost8, 'Poisson');
histfit(DataPost8,11,'Poisson')
PoiPost8 = mle(DataPost8, 'distribution', 'Poisson');
xlabel('Parasite Count')
ylabel('Density')
title('Poisson Distribution')

linkaxes([R81,R83,R84],'xy')
