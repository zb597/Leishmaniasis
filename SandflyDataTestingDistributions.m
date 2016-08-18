%% mouse A -  sandfly data
PreA = [27.1983	0.1372	0.0722	0.311	0.0185	0.1038	0	0.0115	2.5457	0	0	0.0676	0.0383	0	0.0293	0	0.0557	0	0.0216	0.4142	0.0604	7.8228	0	0	0.0611	0.0848	0.0205	0.0065	0.1266	0.0255	0.1392	0.0452	0.0358	0	352.9986	0.0434	395.6194	0.1813	0.4578	0.0797];
PostA = [2330.658	0.211676	0.012102	1.382218	0.270469	496.018	0.129177	9216.162	14183.15	594.3729	0.972096	108.5778	24313.63	1.979764	6015.981	0.23771	3.772378	3.264987	7.608706	39334.38	39.58626	4877.153	18.16791	0.088745	210.875	0.661749	1774.815	136537	0.139089	52032.7	1.851245	0.431662	0.135417	192.2832	551.4379	0.583188	267.3093	328.6461	0.177631	13412.1];

%data
DataPreA = round(reshape(PreA, [numel(PreA), 1]));
DataPostA = round(reshape(PostA, [numel(PostA), 1]));

%log likelihood
norPreA = fitdist(DataPreA, 'Normal');
xPreA = norPreA.NLogL; %log likelihood (normal dist.)
norPostA = fitdist(DataPostA, 'Normal');
xPostA = norPostA.NLogL; %log likelihood (normal dist.)

nbPreA = fitdist(DataPreA, 'Negative Binomial');
yPreA = nbPreA.NLogL; %log likelihood (Neg. Bin. dist.)
nbPostA = fitdist(DataPostA, 'Negative Binomial');
yPostA = nbPostA.NLogL; %log likelihood (Neg. Bin. dist.)

expPreA = fitdist(DataPreA, 'Exponential');
zPreA = expPreA.NLogL; %log likelihood (exponential dist.)
expPostA = fitdist(DataPostA, 'Exponential');
zPostA = expPostA.NLogL; %log likelihood (exponential dist.)

poiPreA = fitdist(DataPreA, 'Poisson');
wPreA = poiPreA.NLogL; %log likelihood (poisson dist.)
poiPostA = fitdist(DataPostA, 'Poisson');
wPostA = poiPostA.NLogL; %log likelihood (poisson dist.)

%Akaike Information Criterion
aicxPreA = 2*2-2*log(xPreA); %AIC = 2*#parameters - 2ln(ll)
aicyPreA = 2*2-2*log(yPreA);
aiczPreA = 2*1-2*log(zPreA);
aicwPreA = 2*1-2*log(wPreA);
aicPreA = aicbic([xPreA, yPreA, zPreA, wPreA],[2, 2, 1, 1]); 
aicxPostA = 2*2-2*log(xPostA); %AIC = 2*#parameters - 2ln(ll)
aicyPostA = 2*2-2*log(yPostA);
aiczPostA = 2*1-2*log(zPostA);
aicwPostA = 2*1-2*log(wPostA);
aicPostA = aicbic([xPostA, yPostA, zPostA, wPostA],[2, 2, 1, 1]);

%Anderson-Darling test
distnPreA = makedist('normal','mu',19.675,'sigma',82.7016);
[h,p] = adtest(DataPreA,'Distribution',distnPreA);
distnPostA = makedist('normal','mu',7671.45,'sigma',23571.9);
[h,p] = adtest(DataPostA,'Distribution',distnPostA);

distnbPreA = makedist('negative binomial','r',0.0195457,'p',0.00099244);
[h,p] = adtest(DataPreA,'Distribution',distnbPreA);
distnbPostA = makedist('negative binomial','r',0.109374,'p',1.42571e-05);
[h,p] = adtest(DataPostA,'Distribution',distnbPostA);

distePreA = makedist('exponential','mu',19.675);
[h,p] = adtest(DataPreA,'Distribution',distePreA);
distePostA = makedist('exponential','mu',7671.45);
[h,p] = adtest(DataPostA,'Distribution',distePostA);

distpPreA = makedist('poisson','lambda',19.675);
[h,p] = adtest(DataPreA,'Distribution',distpPreA);
distpPostA = makedist('poisson','lambda',7671.45);
[h,p] = adtest(DataPostA,'Distribution',distpPostA);

%fitting distributions
figure;
a1 = subplot(2,2,1);
norPreA = fitdist(DataPreA, 'Normal');
histfit(DataPreA,11,'Normal')
NormPreA = mle(DataPreA, 'distribution', 'Normal');
xlabel('Parasite Count')
ylabel('Density')
title('Normal Distribution')
axis([0,400,0,40])
suptitle('Pre Mouse A')

a2 = subplot(2,2,2);
nbPreA = fitdist(DataPreA, 'Negative Binomial');
histfit(DataPreA,11,'Negative Binomial')
NegBinPreA = mle(DataPreA, 'distribution', 'Negative Binomial');
xlabel('Parasite Count')
ylabel('Density')
title('Negative Binomial Distribution')
axis([0,400,0,40])

a3 = subplot(2,2,3);
expPreA = fitdist(DataPreA, 'Exponential');
histfit(DataPreA,11,'Exponential')
ExpoPreA = mle(DataPreA, 'distribution', 'Exponential');
xlabel('Parasite Count')
ylabel('Density')
title('Exponential Distribution')
axis([0,400,0,40])

a4 = subplot(2,2,4);
poiPreA = fitdist(DataPreA, 'Poisson');
histfit(DataPreA,11,'Poisson')
PoiPreA = mle(DataPreA, 'distribution', 'Poisson');
xlabel('Parasite Count')
ylabel('Density')
title('Poisson Distribution')
axis([0,400,0,40])

figure;
A1 = subplot(2,2,1);
norPostA = fitdist(DataPostA, 'Normal');
histfit(DataPostA,11,'Normal')
NormPostA = mle(DataPostA, 'distribution', 'Normal');
xlabel('Parasite Count')
ylabel('Density')
title('Normal Distribution')
xlim([0,inf])
suptitle('Post Mouse A')

A2 = subplot(2,2,2);
nbPostA = fitdist(DataPostA, 'Negative Binomial');
histfit(DataPostA,11,'Negative Binomial')
NegBinPostA = mle(DataPostA, 'distribution', 'Negative Binomial');
xlabel('Parasite Count')
ylabel('Density')
title('Negative Binomial Distribution')
axis([0,136500,0,40])

A3 = subplot(2,2,3);
expPostA = fitdist(DataPostA, 'Exponential');
histfit(DataPostA,11,'Exponential')
ExpoPostA = mle(DataPostA, 'distribution', 'Exponential');
xlabel('Parasite Count')
ylabel('Density')
title('Exponential Distribution')

A4 = subplot(2,2,4);
poiPostA = fitdist(DataPostA, 'Poisson');
histfit(DataPostA,11,'Poisson')
PoiPostA = mle(DataPostA, 'distribution', 'Poisson');
xlabel('Parasite Count')
ylabel('Density')
title('Poisson Distribution')

linkaxes([A1,A3,A4],'xy')

%% mouse B
PreB = [0	0	0.1416	0	0	0	0	8.8597	0	0	0	0	0	0	9.8638	95.5238	0.0069	265.0301	0	78.2494	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0];
PostB = [0.0078	95.0941	0.0356	108.6664	0.0833	0	0.0333	0.0572	0	0.1846	69.317	0.0142	475.459	0.114	0.2222	0.1291	0.0523	136.3398	0	0.0821	518.3716	0.0084	0.2384	0.1386	0.3097	4.348	186.5302	0	65207.98	556.7891	116.8306	0.0101	1.2092	229.806	0.0399	1554.047	0	0	1.3197	86.1103];

%data
DataPreB = round(reshape(PreB, [numel(PreB), 1]));
DataPostB = round(reshape(PostB, [numel(PostB), 1]));

%log likelihood
norPreB = fitdist(DataPreB, 'Normal');
xPreB = norPreB.NLogL; %log likelihood (normal dist.)
norPostB = fitdist(DataPostB, 'Normal');
xPostB = norPostB.NLogL; %log likelihood (normal dist.)

nbPreB = fitdist(DataPreB, 'Negative Binomial');
yPreB = nbPreB.NLogL; %log likelihood (Neg. Bin. dist.)
nbPostB = fitdist(DataPostB, 'Negative Binomial');
yPostB = nbPostB.NLogL; %log likelihood (Neg. Bin. dist.)

expPreB = fitdist(DataPreB, 'Exponential');
zPreB = expPreB.NLogL; %log likelihood (exponential dist.)
expPostB = fitdist(DataPostB, 'Exponential');
zPostB = expPostB.NLogL; %log likelihood (exponential dist.)

poiPreB = fitdist(DataPreB, 'Poisson');
wPreB = poiPreB.NLogL; %log likelihood (poisson dist.)
poiPostB = fitdist(DataPostB, 'Poisson');
wPostB = poiPostB.NLogL; %log likelihood (poisson dist.)

%Akaike Information Criterion
aicxPreB = 2*2-2*log(xPreB); %AIC = 2*#parameters - 2ln(ll)
aicyPreB = 2*2-2*log(yPreB);
aiczPreB = 2*1-2*log(zPreB);
aicwPreB = 2*1-2*log(wPreB);
aicPreB = aicbic([xPreB, yPreB, zPreB, wPreB],[2, 2, 1, 1]); 
aicxPostB = 2*2-2*log(xPostB); %AIC = 2*#parameters - 2ln(ll)
aicyPostB = 2*2-2*log(yPostB);
aiczPostB = 2*1-2*log(zPostB);
aicwPostB = 2*1-2*log(wPostB);
aicPostB = aicbic([xPostB, yPostB, zPostB, wPostB],[2, 2, 1, 1]);

%Anderson-Darling test
distnPreB = makedist('normal','mu',11.45,'sigma',45.4216);
[h,p] = adtest(DataPreB,'Distribution',distnPreB);
distnPostB = makedist('normal','mu',1733.68,'sigma',10297.2);
[h,p] = adtest(DataPostB,'Distribution',distnPostB);

distnbPreB = makedist('negative binomial','r',0.021837,'p',0.00190353);
[h,p] = adtest(DataPreB,'Distribution',distnbPreB);
distnbPostB = makedist('negative binomial','r',0.0475553,'p',2.74296e-05);
[h,p] = adtest(DataPostB,'Distribution',distnbPostB);

distePreB = makedist('exponential','mu',11.45);
[h,p] = adtest(DataPreB,'Distribution',distePreB);
distePostB = makedist('exponential','mu',1733.68);
[h,p] = adtest(DataPostB,'Distribution',distePostB);

distpPreB = makedist('poisson','lambda',11.45);
[h,p] = adtest(DataPreB,'Distribution',distpPreB);
distpPostB = makedist('poisson','lambda',1733.68);
[h,p] = adtest(DataPostB,'Distribution',distpPostB);

% fitting distributions
figure;
b1 = subplot(2,2,1);
norPreB = fitdist(DataPreB, 'Normal');
histfit(DataPreB,11,'Normal')
NormPreB = mle(DataPreB, 'distribution', 'Normal');
xlabel('Parasite Count')
ylabel('Density')
title('Normal Distribution')
axis([0,400,0,40])
suptitle('Pre Mouse B')

b2 = subplot(2,2,2);
nbPreB = fitdist(DataPreB, 'Negative Binomial');
histfit(DataPreB,11,'Negative Binomial')
NegBinPreB = mle(DataPreB, 'distribution', 'Negative Binomial');
xlabel('Parasite Count')
ylabel('Density')
title('Negative Binomial Distribution')
axis([0,400,0,40])

b3 = subplot(2,2,3);
expPreB = fitdist(DataPreB, 'Exponential');
histfit(DataPreB,11,'Exponential')
ExpoPreB = mle(DataPreB, 'distribution', 'Exponential');
xlabel('Parasite Count')
ylabel('Density')
title('Exponential Distribution')
axis([0,400,0,40])

b4 = subplot(2,2,4);
poiPreB = fitdist(DataPreB, 'Poisson');
histfit(DataPreB,11,'Poisson')
PoiPreB = mle(DataPreB, 'distribution', 'Poisson');
xlabel('Parasite Count')
ylabel('Density')
title('Poisson Distribution')
axis([0,400,0,40])

figure;
B1 = subplot(2,2,1);
norPostB = fitdist(DataPostB, 'Normal');
histfit(DataPostB,11,'Normal')
NormPostB = mle(DataPostB, 'distribution', 'Normal');
xlabel('Parasite Count')
ylabel('Density')
title('Normal Distribution')
xlim([0,inf])
suptitle('Post Mouse B')

B2 = subplot(2,2,2);
nbPostB = fitdist(DataPostB, 'Negative Binomial');
histfit(DataPostB,11,'Negative Binomial')
NegBinPostB = mle(DataPostB, 'distribution', 'Negative Binomial');
xlabel('Parasite Count')
ylabel('Density')
title('Negative Binomial Distribution')
axis([0,63000,0,40])

B3 = subplot(2,2,3);
expPostB = fitdist(DataPostB, 'Exponential');
histfit(DataPostB,11,'Exponential')
ExpoPostB = mle(DataPostB, 'distribution', 'Exponential');
xlabel('Parasite Count')
ylabel('Density')
title('Exponential Distribution')

B4 = subplot(2,2,4);
poiPostB = fitdist(DataPostB, 'Poisson');
histfit(DataPostB,11,'Poisson')
PoiPostB = mle(DataPostB, 'distribution', 'Poisson');
xlabel('Parasite Count')
ylabel('Density')
title('Poisson Distribution')

linkaxes([B1,B3,B4],'xy')

%% mouse C
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
axis([0,56700,0,40])

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

%% mouse D
PreD = [0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0.2714	0];
PostD = [0	0	0	0.2034	0.0523	0.5188	0.0928	0	0	0.4741	0	0	0	16.005	0	0.1429	0.0364	0.0324	0.5092	0.1279	0.0303	1.8155	0.225	0.1006	0	0.0295	0.0435	0	0.3187	0.3095	0.1015	0.0575	0	0	0	0	0	0.1158	0.1267	0];

%data
DataPreD = round(reshape(PreD, [numel(PreD), 1]));
DataPostD = round(reshape(PostD, [numel(PostD), 1]));

%log likelihood
% norPreD = fitdist(DataPreD, 'Normal'); %%%Data = zero vector
% xPreD = norPreD.NLogL; %log likelihood (normal dist.)
norPostD = fitdist(DataPostD, 'Normal');
xPostD = norPostD.NLogL; %log likelihood (normal dist.)

% nbPreD = fitdist(DataPreD, 'Negative Binomial'); %%%Data = zero vector
% yPreD = nbPreD.NLogL; %log likelihood (Neg. Bin. dist.)
nbPostD = fitdist(DataPostD, 'Negative Binomial');
yPostD = nbPostD.NLogL; %log likelihood (Neg. Bin. dist.)

% expPreD = fitdist(DataPreD, 'Exponential'); %%%Data = zero vector
% zPreD = expPreD.NLogL; %log likelihood (exponential dist.)
expPostD = fitdist(DataPostD, 'Exponential');
zPostD = expPostD.NLogL; %log likelihood (exponential dist.)

% poiPreD = fitdist(DataPreD, 'Poisson'); %%%Data = zero vector
% wPreD = poiPreD.NLogL; %log likelihood (poisson dist.)
poiPostD = fitdist(DataPostD, 'Poisson');
wPostD = poiPostD.NLogL; %log likelihood (poisson dist.)

%Akaike Information Criterion
% aicxPreD = 2*2-2*log(xPreD); %AIC = 2*#parameters - 2ln(ll)
% aicyPreD = 2*2-2*log(yPreD);
% aiczPreD = 2*1-2*log(zPreD);
% aicwPreD = 2*1-2*log(wPreD);
% aicPreD = aicbic([xPreD, yPreD, zPreD, wPreD],[2, 2, 1, 1]); 
aicxPostD = 2*2-2*log(xPostD); %AIC = 2*#parameters - 2ln(ll)
aicyPostD = 2*2-2*log(yPostD);
aiczPostD = 2*1-2*log(zPostD);
aicwPostD = 2*1-2*log(wPostD);
aicPostD = aicbic([xPostD, yPostD, zPostD, wPostD],[2, 2, 1, 1]);

%Anderson-Darling test
% distnPreD = makedist('normal','mu',0,'sigma',0); %%%Data = zero vector
% [h,p] = adtest(DataPreD,'Distribution',distnPreD);
distnPostD = makedist('normal','mu',0.5,'sigma',2.54196);
[h,p] = adtest(DataPostD,'Distribution',distnPostD);

% distnbPreD = makedist('negative binomial','r',Inf,'p',1); %%%Data = zero vector
% [h,p] = adtest(DataPreD,'Distribution',distnbPreD);
distnbPostD = makedist('negative binomial','r',0.0400706,'p',0.0741951);
[h,p] = adtest(DataPostD,'Distribution',distnbPostD);

% distePreD = makedist('exponential','mu',0); %%%Data = zero vector
% [h,p] = adtest(DataPreD,'Distribution',distePreD);
distePostD = makedist('exponential','mu',0.5);
[h,p] = adtest(DataPostD,'Distribution',distePostD);

% distpPreD = makedist('poisson','lambda',0); %%%Data = zero vector
% [h,p] = adtest(DataPreD,'Distribution',distpPreD);
distpPostD = makedist('poisson','lambda',0.5);
[h,p] = adtest(DataPostD,'Distribution',distpPostD);

%fitting distributions
% figure;
% d1 = subplot(2,2,1)
% norPreD = fitdist(DataPreD, 'Normal');
% histfit(DataPreD,11,'Normal')
% NormPreD = mle(DataPreD, 'distribution', 'Normal');
% xlabel('Parasite Count')
% ylabel('Density')
% title('Normal Distribution')
% axis([0,400,0,40])
% suptitle('Pre Mouse D')
% 
% d2 = subplot(2,2,2)
% nbPreD = fitdist(DataPreD, 'Negative Binomial');
% histfit(DataPreD,11,'Negative Binomial')
% NegBinPreD = mle(DataPreD, 'distribution', 'Negative Binomial');
% xlabel('Parasite Count')
% ylabel('Density')
% title('Negative Binomial Distribution')
% axis([0,400,0,40])
% 
% d3 = subplot(2,2,3)
% expPreD = fitdist(DataPreD, 'Exponential');
% histfit(DataPreD,11,'Exponential')
% ExpoPreD = mle(DataPreD, 'distribution', 'Exponential');
% xlabel('Parasite Count')
% ylabel('Density')
% title('Exponential Distribution')
% axis([0,400,0,40])
% 
% d4 = subplot(2,2,4)
% poiPreD = fitdist(DataPreD, 'Poisson');
% histfit(DataPreD,11,'Poisson')
% PoiPreD = mle(DataPreD, 'distribution', 'Poisson');
% xlabel('Parasite Count')
% ylabel('Density')
% title('Poisson Distribution')
% axis([0,400,0,40])

figure;
D1 = subplot(2,2,1);
norPostD = fitdist(DataPostD, 'Normal');
histfit(DataPostD,11,'Normal')
NormPostD = mle(DataPostD, 'distribution', 'Normal');
xlabel('Parasite Count')
ylabel('Density')
title('Normal Distribution')
xlim([0,inf])
suptitle('Post Mouse D')

D2 = subplot(2,2,2);
nbPostD = fitdist(DataPostD, 'Negative Binomial');
histfit(DataPostD,11,'Negative Binomial')
NegBinPostD = mle(DataPostD, 'distribution', 'Negative Binomial');
xlabel('Parasite Count')
ylabel('Density')
title('Negative Binomial Distribution')
axis([0,16,0,40])

D3 = subplot(2,2,3);
expPostD = fitdist(DataPostD, 'Exponential');
histfit(DataPostD,11,'Exponential')
ExpoPostD = mle(DataPostD, 'distribution', 'Exponential');
xlabel('Parasite Count')
ylabel('Density')
title('Exponential Distribution')

D4 = subplot(2,2,4);
poiPostD = fitdist(DataPostD, 'Poisson');
histfit(DataPostD,11,'Poisson')
PoiPostD = mle(DataPostD, 'distribution', 'Poisson');
xlabel('Parasite Count')
ylabel('Density')
title('Poisson Distribution')

linkaxes([D1,D3,D4],'xy')

%% mouse E
PreE = [0	0	0	0	0	0	0	0	0.0623	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0];
PostE = [0	0.0998	0	0	0	0	0.0292	0	0.065	0	0	10.0138	0	0.1518	0	0	0.9843	0	0	0.0125	0	0	0	0.0103	0.2649	0.0077	4.9609	0	0	0	0	0	0	0	0	0	0	0	0	0];

%data
DataPreE = round(reshape(PreE, [numel(PreE), 1]));
DataPostE = round(reshape(PostE, [numel(PostE), 1]));

%log likelihood
% norPreE = fitdist(DataPreE, 'Normal'); %%%Data = zero vector
% xPreE = norPreE.NLogL; %log likelihood (normal dist.)
norPostE = fitdist(DataPostE, 'Normal');
xPostE = norPostE.NLogL; %log likelihood (normal dist.)

% nbPreE = fitdist(DataPreE, 'Negative Binomial'); %%%Data = zero vector
% yPreE = nbPreE.NLogL; %log likelihood (Neg. Bin. dist.)
nbPostE = fitdist(DataPostE, 'Negative Binomial');
yPostE = nbPostE.NLogL; %log likelihood (Neg. Bin. dist.)

% expPreE = fitdist(DataPreE, 'Exponential'); %%%Data = zero vector
% zPreE = expPreE.NLogL; %log likelihood (exponential dist.)
expPostE = fitdist(DataPostE, 'Exponential');
zPostE = expPostE.NLogL; %log likelihood (exponential dist.)

% poiPreE = fitdist(DataPreE, 'Poisson'); %%%Data = zero vector
% wPreE = poiPreE.NLogL; %log likelihood (poisson dist.)
poiPostE = fitdist(DataPostE, 'Poisson');
wPostE = poiPostE.NLogL; %log likelihood (poisson dist.)

%Akaike Information Criterion
% aicxPreE = 2*2-2*log(xPreE); %AIC = 2*#parameters - 2ln(ll)
% aicyPreE = 2*2-2*log(yPreE);
% aiczPreE = 2*1-2*log(zPreE);
% aicwPreE = 2*1-2*log(wPreE);
% aicPreE = aicbic([xPreE, yPreE, zPreE, wPreE],[2, 2, 1, 1]); 
aicxPostE = 2*2-2*log(xPostE); %AIC = 2*#parameters - 2ln(ll)
aicyPostE = 2*2-2*log(yPostE);
aiczPostE = 2*1-2*log(zPostE);
aicwPostE = 2*1-2*log(wPostE);
aicPostE = aicbic([xPostE, yPostE, zPostE, wPostE],[2, 2, 1, 1]);

%Anderson-Darling test
% distnPreE = makedist('normal','mu',0,'sigma',0); %%%Data = zero vector
% [h,p] = adtest(DataPreE,'Distribution',distnPreE);
distnPostE = makedist('normal','mu',0.4,'sigma',1.75119);
[h,p] = adtest(DataPostE,'Distribution',distnPostE);

% distnbPreE = makedist('negative binomial','r',Inf,'p',1); %%%Data = zero vector
% [h,p] = adtest(DataPreE,'Distribution',distnbPreE);
distnbPostE = makedist('negative binomial','r',0.0292258,'p',0.0680896);
[h,p] = adtest(DataPostE,'Distribution',distnbPostE);

% distePreE = makedist('exponential','mu',0); %%%Data = zero vector
% [h,p] = adtest(DataPreE,'Distribution',distePreE);
distePostE = makedist('exponential','mu',0.4);
[h,p] = adtest(DataPostE,'Distribution',distePostE);

% distpPreE = makedist('poisson','lambda',0); %%%Data = zero vector
% [h,p] = adtest(DataPreE,'Distribution',distpPreE);
distpPostE = makedist('poisson','lambda',0.4);
[h,p] = adtest(DataPostE,'Distribution',distpPostE);

% %fitting distributions
% figure;
% e1 = subplot(2,2,1)
% norPreE = fitdist(DataPreE, 'Normal');
% histfit(DataPreE,11,'Normal')
% NormPreE = mle(DataPreE, 'distribution', 'Normal');
% xlabel('Parasite Count')
% ylabel('Density')
% title('Normal Distribution')
% axis([0,400,0,40])
% suptitle('Pre Mouse E')
% 
% e2 = subplot(2,2,2)
% nbPreE = fitdist(DataPreE, 'Negative Binomial');
% histfit(DataPreE,11,'Negative Binomial')
% NegBinPreE = mle(DataPreE, 'distribution', 'Negative Binomial');
% xlabel('Parasite Count')
% ylabel('Density')
% title('Negative Binomial Distribution')
% axis([0,400,0,40])
% 
% e3 = subplot(2,2,3)
% expPreE = fitdist(DataPreE, 'Exponential');
% histfit(DataPreE,11,'Exponential')
% ExpoPreE = mle(DataPreE, 'distribution', 'Exponential');
% xlabel('Parasite Count')
% ylabel('Density')
% title('Exponential Distribution')
% axis([0,400,0,40])
% 
% e4 = subplot(2,2,4)
% poiPreE = fitdist(DataPreE, 'Poisson');
% histfit(DataPreE,11,'Poisson')
% PoiPreE = mle(DataPreE, 'distribution', 'Poisson');
% xlabel('Parasite Count')
% ylabel('Density')
% title('Poisson Distribution')
% axis([0,400,0,40])

figure;
E1 = subplot(2,2,1);
norPostE = fitdist(DataPostE, 'Normal');
histfit(DataPostE,11,'Normal')
NormPostE = mle(DataPostE, 'distribution', 'Normal');
xlabel('Parasite Count')
ylabel('Density')
title('Normal Distribution')
xlim([0,inf])
suptitle('Post Mouse E')

E2 = subplot(2,2,2);
nbPostE = fitdist(DataPostE, 'Negative Binomial');
histfit(DataPostE,11,'Negative Binomial')
NegBinPostE = mle(DataPostE, 'distribution', 'Negative Binomial');
xlabel('Parasite Count')
ylabel('Density')
title('Negative Binomial Distribution')
axis([0,10,0,40])

E3 = subplot(2,2,3);
expPostE = fitdist(DataPostE, 'Exponential');
histfit(DataPostE,11,'Exponential')
ExpoPostE = mle(DataPostE, 'distribution', 'Exponential');
xlabel('Parasite Count')
ylabel('Density')
title('Exponential Distribution')

E4 = subplot(2,2,4);
poiPostE = fitdist(DataPostE, 'Poisson');
histfit(DataPostE,11,'Poisson')
PoiPostE = mle(DataPostE, 'distribution', 'Poisson');
xlabel('Parasite Count')
ylabel('Density')
title('Poisson Distribution')

linkaxes([E1,E3,E4],'xy')

%% mouse F
PreF = [0	0	0	0	0	0	0.0254	0	0	0	0	0	0	0	0	0	0	0.0089	0	0	0	0	0	5.2044	0.006	0	0	0	0	0	0.0094	0	0	0	0	0	0	0	0	0];
PostF = [0	0	0	0	0	0	0	0	0	0	0	0	0.1917	0	0	0	0	0	0.0265	0	0	0.3034	0	0	0	0	0	0	0	0	0	7.113	0	0	0	0	0	0	0	6.349];

%data
DataPreF = round(reshape(PreF, [numel(PreF), 1]));
DataPostF = round(reshape(PostF, [numel(PostF), 1]));

%log likelihood
norPreF = fitdist(DataPreF, 'Normal');
xPreF = norPreF.NLogL; %log likelihood (normal dist.)
norPostF = fitdist(DataPostF, 'Normal');
xPostF = norPostF.NLogL; %log likelihood (normal dist.)

nbPreF = fitdist(DataPreF, 'Negative Binomial');
yPreF = nbPreF.NLogL; %log likelihood (Neg. Bin. dist.)
nbPostF = fitdist(DataPostF, 'Negative Binomial');
yPostF = nbPostF.NLogL; %log likelihood (Neg. Bin. dist.)

expPreF = fitdist(DataPreF, 'Exponential');
zPreF = expPreF.NLogL; %log likelihood (exponential dist.)
expPostF = fitdist(DataPostF, 'Exponential');
zPostF = expPostF.NLogL; %log likelihood (exponential dist.)

poiPreF = fitdist(DataPreF, 'Poisson');
wPreF = poiPreF.NLogL; %log likelihood (poisson dist.)
poiPostF = fitdist(DataPostF, 'Poisson');
wPostF = poiPostF.NLogL; %log likelihood (poisson dist.)

%Akaike Information Criterion
aicxPreF = 2*2-2*log(xPreF); %AIC = 2*#parameters - 2ln(ll)
aicyPreF = 2*2-2*log(yPreF);
aiczPreF = 2*1-2*log(zPreF);
aicwPreF = 2*1-2*log(wPreF);
aicPreF = aicbic([xPreF, yPreF, zPreF, wPreF],[2, 2, 1, 1]); 
aicxPostF = 2*2-2*log(xPostF); %AIC = 2*#parameters - 2ln(ll)
aicyPostF = 2*2-2*log(yPostF);
aiczPostF = 2*1-2*log(zPostF);
aicwPostF = 2*1-2*log(wPostF);
aicPostF = aicbic([xPostF, yPostF, zPostF, wPostF],[2, 2, 1, 1]);

%Anderson-Darling test
distnPreF = makedist('normal','mu',0.125,'sigma',0.790569);
[h,p] = adtest(DataPreF,'Distribution',distnPreF);
distnPostF = makedist('normal','mu',0.325,'sigma',1.43915);
[h,p] = adtest(DataPostF,'Distribution',distnPostF);

distnbPreF = makedist('negative binomial','r',0.00968933,'p',0.0719384);
[h,p] = adtest(DataPreF,'Distribution',distnbPreF);
distnbPostF = makedist('negative binomial','r',0.0175009,'p',0.0510974);
[h,p] = adtest(DataPostF,'Distribution',distnbPostF);

distePreF = makedist('exponential','mu',0.125);
[h,p] = adtest(DataPreF,'Distribution',distePreF);
distePostF = makedist('exponential','mu',0.325);
[h,p] = adtest(DataPostF,'Distribution',distePostF);

distpPreF = makedist('poisson','lambda',0.125);
[h,p] = adtest(DataPreF,'Distribution',distpPreF);
distpPostF = makedist('poisson','lambda',0.325);
[h,p] = adtest(DataPostF,'Distribution',distpPostF);

%fitting distributions
figure;
f1 = subplot(2,2,1);
norPreF = fitdist(DataPreF, 'Normal');
histfit(DataPreF,11,'Normal')
NormPreF = mle(DataPreF, 'distribution', 'Normal');
xlabel('Parasite Count')
ylabel('Density')
title('Normal Distribution')
xlim([0,inf])
suptitle('Pre Mouse F')

f2 = subplot(2,2,2);
nbPreF = fitdist(DataPreF, 'Negative Binomial');
histfit(DataPreF,11,'Negative Binomial')
NegBinPreF = mle(DataPreF, 'distribution', 'Negative Binomial');
xlabel('Parasite Count')
ylabel('Density')
title('Negative Binomial Distribution')
axis([0,5,0,40])

f3 = subplot(2,2,3);
expPreF = fitdist(DataPreF, 'Exponential');
histfit(DataPreF,11,'Exponential')
ExpoPreF = mle(DataPreF, 'distribution', 'Exponential');
xlabel('Parasite Count')
ylabel('Density')
title('Exponential Distribution')

f4 = subplot(2,2,4);
poiPreF = fitdist(DataPreF, 'Poisson');
histfit(DataPreF,11,'Poisson')
PoiPreF = mle(DataPreF, 'distribution', 'Poisson');
xlabel('Parasite Count')
ylabel('Density')
title('Poisson Distribution')

linkaxes([f1,f3,f4],'xy')

figure;
F1 = subplot(2,2,1);
norPostF = fitdist(DataPostF, 'Normal');
histfit(DataPostF,11,'Normal')
NormPostF = mle(DataPostF, 'distribution', 'Normal');
xlabel('Parasite Count')
ylabel('Density')
title('Normal Distribution')
xlim([0,inf])
suptitle('Post Mouse F')

F2 = subplot(2,2,2);
nbPostF = fitdist(DataPostF, 'Negative Binomial');
histfit(DataPostF,11,'Negative Binomial')
NegBinPostF = mle(DataPostF, 'distribution', 'Negative Binomial');
xlabel('Parasite Count')
ylabel('Density')
title('Negative Binomial Distribution')
axis([0,7,0,40])

F3 = subplot(2,2,3);
expPostF = fitdist(DataPostF, 'Exponential');
histfit(DataPostF,11,'Exponential')
ExpoPostF = mle(DataPostF, 'distribution', 'Exponential');
xlabel('Parasite Count')
ylabel('Density')
title('Exponential Distribution')

F4 = subplot(2,2,4);
poiPostF = fitdist(DataPostF, 'Poisson');
histfit(DataPostF,11,'Poisson')
PoiPostF = mle(DataPostF, 'distribution', 'Poisson');
xlabel('Parasite Count')
ylabel('Density')
title('Poisson Distribution')

linkaxes([F1,F3,F4],'xy')

%% mouse G
PreG = [0.0273	0	0	0	0	0	0	0	0	0	0.1284	0	0	0	0	0	0	0	0	0	0	0	0	1.2784	0	0	0.1589	0	0	0	0	0	0	0	0	0	0	0	0	0];
PostG = [0.8379	15.1819	0	0.1302	0.3051	0.2715	0.0485	0.039	0.0824	0	0.0719	0.0441	0	0.1073	0.0954	0	0.1162	0.2096	0.0288	0	0	0.2031	0	0.138	0.0333	0.1276	0.2739	0.0227	0.0139	0.1696	0.5549	0	0	0.012	0.0358	0	0	0.034	25.0815	0.4392];

%data
DataPreG = round(reshape(PreG, [numel(PreG), 1]));
DataPostG = round(reshape(PostG, [numel(PostG), 1]));

%log likelihood
norPreG = fitdist(DataPreG, 'Normal');
xPreG = norPreG.NLogL; %log likelihood (normal dist.)
norPostG = fitdist(DataPostG, 'Normal');
xPostG = norPostG.NLogL; %log likelihood (normal dist.)

% nbPreG = fitdist(DataPreG, 'Negative Binomial'); %%%Issue!!
% yPreG = nbPreG.NLogL; %log likelihood (Neg. Bin. dist.)
nbPostG = fitdist(DataPostG, 'Negative Binomial');
yPostG = nbPostG.NLogL; %log likelihood (Neg. Bin. dist.)

expPreG = fitdist(DataPreG, 'Exponential');
zPreG = expPreG.NLogL; %log likelihood (exponential dist.)
expPostG = fitdist(DataPostG, 'Exponential');
zPostG = expPostG.NLogL; %log likelihood (exponential dist.)

poiPreG = fitdist(DataPreG, 'Poisson');
wPreG = poiPreG.NLogL; %log likelihood (poisson dist.)
poiPostG = fitdist(DataPostG, 'Poisson');
wPostG = poiPostG.NLogL; %log likelihood (poisson dist.)

%Akaike Information Criterion
aicxPreG = 2*2-2*log(xPreG); %AIC = 2*#parameters - 2ln(ll)
% aicyPreG = 2*2-2*log(yPreG);
aiczPreG = 2*1-2*log(zPreG);
aicwPreG = 2*1-2*log(wPreG);
aicPreG = aicbic([xPreG, zPreG, wPreG],[2, 1, 1]); 
aicxPostG = 2*2-2*log(xPostG); %AIC = 2*#parameters - 2ln(ll)
aicyPostG = 2*2-2*log(yPostG);
aiczPostG = 2*1-2*log(zPostG);
aicwPostG = 2*1-2*log(wPostG);
aicPostG = aicbic([xPostG, yPostG, zPostG, wPostG],[2, 2, 1, 1]);

%Anderson-Darling test
distnPreG = makedist('normal','mu',0.025,'sigma',0.158114);
[h,p] = adtest(DataPreG,'Distribution',distnPreG);
distnPostG = makedist('normal','mu',1.05,'sigma',4.55142);
[h,p] = adtest(DataPostG,'Distribution',distnPostG);

% distnbPreG = makedist('negative binomial','r',Inf,'p',1); %%%Issue!!
% [h,p] = adtest(DataPreG,'Distribution',distnbPreG);
distnbPostG = makedist('negative binomial','r',0.0290615,'p',0.0269322);
[h,p] = adtest(DataPostG,'Distribution',distnbPostG);

distePreG = makedist('exponential','mu',0.025);
[h,p] = adtest(DataPreG,'Distribution',distePreG);
distePostG = makedist('exponential','mu',1.05);
[h,p] = adtest(DataPostG,'Distribution',distePostG);

distpPreG = makedist('poisson','lambda',0.025);
[h,p] = adtest(DataPreG,'Distribution',distpPreG);
distpPostG = makedist('poisson','lambda',1.05);
[h,p] = adtest(DataPostG,'Distribution',distpPostG);

%fitting distributions
figure;
g1 = subplot(3,1,1);
norPreG = fitdist(DataPreG, 'Normal');
histfit(DataPreG,11,'Normal')
NormPreG = mle(DataPreG, 'distribution', 'Normal');
xlabel('Parasite Count')
ylabel('Density')
title('Normal Distribution')
xlim([0,Inf])
suptitle('Pre Mouse G')

g2 = subplot(3,1,2);
expPreG = fitdist(DataPreG, 'Exponential');
histfit(DataPreG,11,'Exponential')
ExpoPreG = mle(DataPreG, 'distribution', 'Exponential');
xlabel('Parasite Count')
ylabel('Density')
title('Exponential Distribution')

g3 = subplot(3,1,3);
poiPreG = fitdist(DataPreG, 'Poisson');
histfit(DataPreG,11,'Poisson')
PoiPreG = mle(DataPreG, 'distribution', 'Poisson');
xlabel('Parasite Count')
ylabel('Density')
title('Poisson Distribution')

linkaxes([g1,g2,g3],'xy')

figure;
G1 = subplot(2,2,1);
norPostG = fitdist(DataPostG, 'Normal');
histfit(DataPostG,11,'Normal')
NormPostG = mle(DataPostG, 'distribution', 'Normal');
xlabel('Parasite Count')
ylabel('Density')
title('Normal Distribution')
xlim([0,inf])
suptitle('Post Mouse G')

G2 = subplot(2,2,2);
nbPostG = fitdist(DataPostG, 'Negative Binomial');
histfit(DataPostG,11,'Negative Binomial')
NegBinPostG = mle(DataPostG, 'distribution', 'Negative Binomial');
xlabel('Parasite Count')
ylabel('Density')
axis([0,25,0,40])
title('Negative Binomial Distribution')

G3 = subplot(2,2,3);
expPostG = fitdist(DataPostG, 'Exponential');
histfit(DataPostG,11,'Exponential')
ExpoPostG = mle(DataPostG, 'distribution', 'Exponential');
xlabel('Parasite Count')
ylabel('Density')
title('Exponential Distribution')

G4 = subplot(2,2,4);
poiPostG = fitdist(DataPostG, 'Poisson');
histfit(DataPostG,11,'Poisson')
PoiPostG = mle(DataPostG, 'distribution', 'Poisson');
xlabel('Parasite Count')
ylabel('Density')
title('Poisson Distribution')

linkaxes([G1,G3,G4],'xy')

%% mouse H
PreH = [0.084	0	0	0.3188	0	0.1276	0.0043	0.0352	0	0	0	0.0917	0.133	0.011	0.003	0	0	0.0064	0	0.1357	0.0379	0.0849	0.4297	0.0392	0.0363	0	0	0	0.0123	0.0093	0	0.0135	0	0.1318	0.0034	0.0097	0.0068	0.0304	0	0];
PostH = [0	0.6721	0.2423	0.0186	0	0.0569	0.73	3.7859	0.0178	0.0337	0.1867	0	0.0242	0.0635	0.0323	0.3772	0.0291	0.0389	0	0	0.1153	0	0	0.0232	0	0	0.5628	0.2002	0.0116	0	0.0723	0.0635	0	0	0.0966	0.0132	0	0	0.1435	0.059];

%data
DataPreH = round(reshape(PreH, [numel(PreH), 1]));
DataPostH = round(reshape(PostH, [numel(PostH), 1]));

%log likelihood
% norPreH = fitdist(DataPreH, 'Normal'); %%%Data = zero vector
% xPreH = norPreH.NLogL; %log likelihood (normal dist.)
norPostH = fitdist(DataPostH, 'Normal');
xPostH = norPostH.NLogL; %log likelihood (normal dist.)

% nbPreH = fitdist(DataPreH, 'Negative Binomial'); %%%Data = zero vector
% yPreH = nbPreH.NLogL; %log likelihood (Neg. Bin. dist.)
nbPostH = fitdist(DataPostH, 'Negative Binomial');
yPostH = nbPostH.NLogL; %log likelihood (Neg. Bin. dist.)

% expPreH = fitdist(DataPreH, 'Exponential'); %%%Data = zero vector
% zPreH = expPreH.NLogL; %log likelihood (exponential dist.)
expPostH = fitdist(DataPostH, 'Exponential');
zPostH = expPostH.NLogL; %log likelihood (exponential dist.)

% poiPreH = fitdist(DataPreH, 'Poisson'); %%%Data = zero vector
% wPreH = poiPreH.NLogL; %log likelihood (poisson dist.)
poiPostH = fitdist(DataPostH, 'Poisson');
wPostH = poiPostH.NLogL; %log likelihood (poisson dist.)

%Akaike Information Criterion
% aicxPreH = 2*2-2*log(xPreH); %AIC = 2*#parameters - 2ln(ll)
% aicyPreH = 2*2-2*log(yPreH);
% aiczPreH = 2*1-2*log(zPreH);
% aicwPreH = 2*1-2*log(wPreH);
% aicPreH = aicbic([xPreH, yPreH, zPreH, wPreH],[2, 2, 1, 1]); 
aicxPostH = 2*2-2*log(xPostH); %AIC = 2*#parameters - 2ln(ll)
aicyPostH = 2*2-2*log(yPostH);
aiczPostH = 2*1-2*log(zPostH);
aicwPostH = 2*1-2*log(wPostH);
aicPostH = aicbic([xPostH, yPostH, zPostH, wPostH],[2, 2, 1, 1]);

%Anderson-Darling test
% distnPreH = makedist('normal','mu',0,'sigma',0); %%%Data = zero vector
% [h,p] = adtest(DataPreH,'Distribution',distnPreH);
distnPostH = makedist('normal','mu',0.175,'sigma',0.675107);
[h,p] = adtest(DataPostH,'Distribution',distnPostH);

% distnbPreH = makedist('negative binomial','r',Inf,'p',1); %%%Data = zero vector
% [h,p] = adtest(DataPreH,'Distribution',distnbPreH);
distnbPostH= makedist('negative binomial','r',0.109813,'p',0.385562);
[h,p] = adtest(DataPostH,'Distribution',distnbPostH);

% distePreH = makedist('exponential','mu',0); %%%Data = zero vector
% [h,p] = adtest(DataPreH,'Distribution',distePreH);
distePostH = makedist('exponential','mu',0.175);
[h,p] = adtest(DataPostB,'Distribution',distePostB);

% distpPreH = makedist('poisson','lambda',0); %%%Data = zero vector
% [h,p] = adtest(DataPreH,'Distribution',distpPreH);
distpPostH = makedist('poisson','lambda',0.175);
[h,p] = adtest(DataPostH,'Distribution',distpPostH);

%fitting distributions
% figure;
% h1 = subplot(2,2,1)
% norPreH = fitdist(DataPreH, 'Normal');
% histfit(DataPreH,11,'Normal')
% NormPreH = mle(DataPreH, 'distribution', 'Normal');
% xlabel('Parasite Count')
% ylabel('Density')
% title('Normal Distribution')
% axis([0,400,0,40])
% suptitle('Pre Mouse H')
% 
% h2 = subplot(2,2,2)
% nbPreH = fitdist(DataPreH, 'Negative Binomial');
% histfit(DataPreH,11,'Negative Binomial')
% NegBinPreH = mle(DataPreH, 'distribution', 'Negative Binomial');
% xlabel('Parasite Count')
% ylabel('Density')
% title('Negative Binomial Distribution')
% axis([0,400,0,40])
% 
% h3 = subplot(2,2,3)
% expPreH = fitdist(DataPreH, 'Exponential');
% histfit(DataPreH,11,'Exponential')
% ExpoPreH = mle(DataPreH, 'distribution', 'Exponential');
% xlabel('Parasite Count')
% ylabel('Density')
% title('Exponential Distribution')
% axis([0,400,0,40])
% 
% h4 = subplot(2,2,4)
% poiPreH = fitdist(DataPreH, 'Poisson');
% histfit(DataPreH,11,'Poisson')
% PoiPreH = mle(DataPreH, 'distribution', 'Poisson');
% xlabel('Parasite Count')
% ylabel('Density')
% title('Poisson Distribution')
% axis([0,400,0,40])

figure;
H1 = subplot(2,2,1);
norPostH = fitdist(DataPostH, 'Normal');
histfit(DataPostH,11,'Normal')
NormPostH = mle(DataPostH, 'distribution', 'Normal');
xlabel('Parasite Count')
ylabel('Density')
title('Normal Distribution')
xlim([0,inf])
suptitle('Post Mouse H')

H2 = subplot(2,2,2);
nbPostH = fitdist(DataPostH, 'Negative Binomial');
histfit(DataPostH,11,'Negative Binomial')
NegBinPostH = mle(DataPostH, 'distribution', 'Negative Binomial');
xlabel('Parasite Count')
ylabel('Density')
title('Negative Binomial Distribution')
axis([0,4,0,40])

H3 = subplot(2,2,3);
expPostH = fitdist(DataPostH, 'Exponential');
histfit(DataPostH,11,'Exponential')
ExpoPostH = mle(DataPostH, 'distribution', 'Exponential');
xlabel('Parasite Count')
ylabel('Density')
title('Exponential Distribution')

H4 = subplot(2,2,4);
poiPostH = fitdist(DataPostH, 'Poisson');
histfit(DataPostH,11,'Poisson')
PoiPostH = mle(DataPostH, 'distribution', 'Poisson');
xlabel('Parasite Count')
ylabel('Density')
title('Poisson Distribution')

linkaxes([H1,H3,H4],'xy')

%% mouse I
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
axis([0,21,0,40])

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
axis([0,273,0,40])

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

linkaxes([I1,I2,I3,I4],'xy')

%% mouse J
PreJ = [0	0	0	0	0	0	0.0899	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	3.4093	0	0	0	0	0	0	0];
PostJ = [2.0842	0	0	0	0	0	41.4374	0	0	0	0	0	0	0.6836	0	0	0	0.541	0	0	0	0	2.6215	0	0	0	0	0	0	0.0198	0	29.5184	0	0	0	0	0	0	0	0];

%data
DataPreJ = round(reshape(PreJ, [numel(PreJ), 1]));
DataPostJ = round(reshape(PostJ, [numel(PostJ), 1]));

%log likelihood
norPreJ = fitdist(DataPreJ, 'Normal');
xPreJ = norPreJ.NLogL; %log likelihood (normal dist.)
norPostJ = fitdist(DataPostJ, 'Normal');
xPostJ = norPostJ.NLogL; %log likelihood (normal dist.)

nbPreJ = fitdist(DataPreJ, 'Negative Binomial');
yPreJ = nbPreJ.NLogL; %log likelihood (Neg. Bin. dist.)
nbPostJ = fitdist(DataPostJ, 'Negative Binomial');
yPostJ = nbPostJ.NLogL; %log likelihood (Neg. Bin. dist.)

expPreJ = fitdist(DataPreJ, 'Exponential');
zPreJ = expPreJ.NLogL; %log likelihood (exponential dist.)
expPostJ = fitdist(DataPostJ, 'Exponential');
zPostJ = expPostJ.NLogL; %log likelihood (exponential dist.)

poiPreJ = fitdist(DataPreJ, 'Poisson');
wPreJ = poiPreJ.NLogL; %log likelihood (poisson dist.)
poiPostJ = fitdist(DataPostJ, 'Poisson');
wPostJ = poiPostJ.NLogL; %log likelihood (poisson dist.)

%Akaike Information Criterion
aicxPreJ = 2*2-2*log(xPreJ); %AIC = 2*#parameters - 2ln(ll)
aicyPreJ = 2*2-2*log(yPreJ);
aiczPreJ = 2*1-2*log(zPreJ);
aicwPreJ = 2*1-2*log(wPreJ);
aicPreJ = aicbic([xPreJ, yPreJ, zPreJ, wPreJ],[2, 2, 1, 1]); 
aicxPostJ = 2*2-2*log(xPostJ); %AIC = 2*#parameters - 2ln(ll)
aicyPostJ = 2*2-2*log(yPostJ);
aiczPostJ = 2*1-2*log(zPostJ);
aicwPostJ = 2*1-2*log(wPostJ);
aicPostJ = aicbic([xPostJ, yPostJ, zPostJ, wPostJ],[2, 2, 1, 1]);

%Anderson-Darling test
distnPreJ = makedist('normal','mu',0.075,'sigma',0.474342);
[h,p] = adtest(DataPreJ,'Distribution',distnPreJ);
distnPostJ = makedist('normal','mu',1.95,'sigma',7.91607);
[h,p] = adtest(DataPostJ,'Distribution',distnPostJ);

distnbPreJ = makedist('negative binomial','r',0.0136176,'p',0.153667);
[h,p] = adtest(DataPreJ,'Distribution',distnbPreJ);
distnbPostJ = makedist('negative binomial','r',0.041603,'p',0.0208892);
[h,p] = adtest(DataPostJ,'Distribution',distnbPostJ);

distePreJ = makedist('exponential','mu',0.075);
[h,p] = adtest(DataPreJ,'Distribution',distePreJ);
distePostJ = makedist('exponential','mu',1.95);
[h,p] = adtest(DataPostJ,'Distribution',distePostJ);

distpPreJ = makedist('poisson','lambda',0.075);
[h,p] = adtest(DataPreJ,'Distribution',distpPreJ);
distpPostJ = makedist('poisson','lambda',1.95);
[h,p] = adtest(DataPostJ,'Distribution',distpPostJ);

%fitting distributions
figure;
j1 = subplot(2,2,1);
norPreJ = fitdist(DataPreJ, 'Normal');
histfit(DataPreJ,11,'Normal')
NormPreJ = mle(DataPreJ, 'distribution', 'Normal');
xlabel('Parasite Count')
ylabel('Density')
title('Normal Distribution')
xlim([0,inf])
suptitle('Pre Mouse J')

j2 = subplot(2,2,2);
nbPreJ = fitdist(DataPreJ, 'Negative Binomial');
histfit(DataPreJ,11,'Negative Binomial')
NegBinPreJ = mle(DataPreJ, 'distribution', 'Negative Binomial');
xlabel('Parasite Count')
ylabel('Density')
title('Negative Binomial Distribution')
axis([0,3,0,40])

j3 = subplot(2,2,3);
expPreJ = fitdist(DataPreJ, 'Exponential');
histfit(DataPreJ,11,'Exponential')
ExpoPreJ = mle(DataPreJ, 'distribution', 'Exponential');
xlabel('Parasite Count')
ylabel('Density')
title('Exponential Distribution')

j4 = subplot(2,2,4);
poiPreJ = fitdist(DataPreJ, 'Poisson');
histfit(DataPreJ,11,'Poisson')
PoiPreJ = mle(DataPreJ, 'distribution', 'Poisson');
xlabel('Parasite Count')
ylabel('Density')
title('Poisson Distribution')

linkaxes([j1,j3,j4],'xy')

figure;
J1 = subplot(2,2,1);
norPostJ = fitdist(DataPostJ, 'Normal');
histfit(DataPostJ,11,'Normal')
NormPostJ = mle(DataPostJ, 'distribution', 'Normal');
xlabel('Parasite Count')
ylabel('Density')
title('Normal Distribution')
xlim([0,inf])
suptitle('Post Mouse J')

J2 = subplot(2,2,2);
nbPostJ = fitdist(DataPostJ, 'Negative Binomial');
histfit(DataPostJ,11,'Negative Binomial')
NegBinPostJ = mle(DataPostJ, 'distribution', 'Negative Binomial');
xlabel('Parasite Count')
ylabel('Density')
title('Negative Binomial Distribution')
axis([0,40,0,40])

J3 = subplot(2,2,3);
expPostJ = fitdist(DataPostJ, 'Exponential');
histfit(DataPostJ,11,'Exponential')
ExpoPostJ = mle(DataPostJ, 'distribution', 'Exponential');
xlabel('Parasite Count')
ylabel('Density')
title('Exponential Distribution')

J4 = subplot(2,2,4);
poiPostJ = fitdist(DataPostJ, 'Poisson');
histfit(DataPostJ,11,'Poisson')
PoiPostJ = mle(DataPostJ, 'distribution', 'Poisson');
xlabel('Parasite Count')
ylabel('Density')
title('Poisson Distribution')

linkaxes([J1,J3,J4],'xy')
