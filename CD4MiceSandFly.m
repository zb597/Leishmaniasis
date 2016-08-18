%CD4 mice

%% Experimental ID A Experiment 2

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
histfit(DataPreA,11,'Normal')
NormPreA = mle(DataPreA, 'distribution', 'Normal');
xlabel('Parasite Count')
ylabel('Density')
title('Normal Distribution')
axis([0,400,0,40])
suptitle('Pre Mouse A')

a2 = subplot(2,2,2);
histfit(DataPreA,11,'Negative Binomial')
NegBinPreA = mle(DataPreA, 'distribution', 'Negative Binomial');
xlabel('Parasite Count')
ylabel('Density')
title('Negative Binomial Distribution')
axis([0,400,0,40])

a3 = subplot(2,2,3);
histfit(DataPreA,11,'Exponential')
ExpoPreA = mle(DataPreA, 'distribution', 'Exponential');
xlabel('Parasite Count')
ylabel('Density')
title('Exponential Distribution')
axis([0,400,0,40])

a4 = subplot(2,2,4);
histfit(DataPreA,11,'Poisson')
PoiPreA = mle(DataPreA, 'distribution', 'Poisson');
xlabel('Parasite Count')
ylabel('Density')
title('Poisson Distribution')
axis([0,400,0,40])

figure;
A1 = subplot(2,2,1);
histfit(DataPostA,11,'Normal')
NormPostA = mle(DataPostA, 'distribution', 'Normal');
xlabel('Parasite Count')
ylabel('Density')
title('Normal Distribution')
xlim([0,inf])
suptitle('Post Mouse A')

A2 = subplot(2,2,2);
histfit(DataPostA,11,'Negative Binomial')
NegBinPostA = mle(DataPostA, 'distribution', 'Negative Binomial');
xlabel('Parasite Count')
ylabel('Density')
title('Negative Binomial Distribution')
axis([0,136500,0,40])

A3 = subplot(2,2,3);
histfit(DataPostA,11,'Exponential')
ExpoPostA = mle(DataPostA, 'distribution', 'Exponential');
xlabel('Parasite Count')
ylabel('Density')
title('Exponential Distribution')

A4 = subplot(2,2,4);
histfit(DataPostA,11,'Poisson')
PoiPostA = mle(DataPostA, 'distribution', 'Poisson');
xlabel('Parasite Count')
ylabel('Density')
title('Poisson Distribution')

linkaxes([A1,A3,A4],'xy')

%% Experimental ID B Experiment 2

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
histfit(DataPreB,11,'Normal')
NormPreB = mle(DataPreB, 'distribution', 'Normal');
xlabel('Parasite Count')
ylabel('Density')
title('Normal Distribution')
axis([0,400,0,40])
suptitle('Pre Mouse B')

b2 = subplot(2,2,2);
histfit(DataPreB,11,'Negative Binomial')
NegBinPreB = mle(DataPreB, 'distribution', 'Negative Binomial');
xlabel('Parasite Count')
ylabel('Density')
title('Negative Binomial Distribution')
axis([0,400,0,40])

b3 = subplot(2,2,3);
histfit(DataPreB,11,'Exponential')
ExpoPreB = mle(DataPreB, 'distribution', 'Exponential');
xlabel('Parasite Count')
ylabel('Density')
title('Exponential Distribution')
axis([0,400,0,40])

b4 = subplot(2,2,4);
histfit(DataPreB,11,'Poisson')
PoiPreB = mle(DataPreB, 'distribution', 'Poisson');
xlabel('Parasite Count')
ylabel('Density')
title('Poisson Distribution')
axis([0,400,0,40])

figure;
B1 = subplot(2,2,1);
histfit(DataPostB,11,'Normal')
NormPostB = mle(DataPostB, 'distribution', 'Normal');
xlabel('Parasite Count')
ylabel('Density')
title('Normal Distribution')
xlim([0,inf])
suptitle('Post Mouse B')

B2 = subplot(2,2,2);
histfit(DataPostB,11,'Negative Binomial')
NegBinPostB = mle(DataPostB, 'distribution', 'Negative Binomial');
xlabel('Parasite Count')
ylabel('Density')
title('Negative Binomial Distribution')
axis([0,63000,0,40])

B3 = subplot(2,2,3);
histfit(DataPostB,11,'Exponential')
ExpoPostB = mle(DataPostB, 'distribution', 'Exponential');
xlabel('Parasite Count')
ylabel('Density')
title('Exponential Distribution')

B4 = subplot(2,2,4);
histfit(DataPostB,11,'Poisson')
PoiPostB = mle(DataPostB, 'distribution', 'Poisson');
xlabel('Parasite Count')
ylabel('Density')
title('Poisson Distribution')

linkaxes([B1,B3,B4],'xy')

%% Experimental ID 1 Experiment 4

Pre14 = [1.4	1.7	0.3	232.7	380	0.1	1096.3	1112.8	5.1	0.1	1319.7	19107.4	0.8	1.2	2617.6	0.2	1.3	9.1	0.8	0.1	0	684.4	308.5	1.7	0	3151.7	0.1	0.2	127.6	0	0.1	0.7	0	11747.5	2.1	3993.1	2.6	0.5	0.8	0.4];
Post14 = [20.6	0.8	190.3	269.9	0	3.2	0.3	0	0.4	83.9	29.9	0.6	31.2	2.4	0	2.7	448.6	5.3	120.7	205.1	0	180.7	0.7	0	0	93.9	0	0	0	84.1	63.9	0	0	0.6	100.5	16.6	0	10.8	1.8	613.3];

%data
DataPre14 = round(reshape(Pre14, [numel(Pre14), 1]));
DataPost14 = round(reshape(Post14, [numel(Post14), 1]));

%log likelihood
norPre14 = fitdist(DataPre14, 'Normal');
xPre14 = norPre14.NLogL; %log likelihood (normal dist.)
norPost14 = fitdist(DataPost14, 'Normal');
xPost14 = norPost14.NLogL; %log likelihood (normal dist.)

nbPre14 = fitdist(DataPre14, 'Negative Binomial');
yPre14 = nbPre14.NLogL; %log likelihood (Neg. Bin. dist.)
nbPost14 = fitdist(DataPost14, 'Negative Binomial');
yPost14 = nbPost14.NLogL; %log likelihood (Neg. Bin. dist.)

expPre14 = fitdist(DataPre14, 'Exponential');
zPre14 = expPre14.NLogL; %log likelihood (exponential dist.)
expPost14 = fitdist(DataPost14, 'Exponential');
zPost14 = expPost14.NLogL; %log likelihood (exponential dist.)

poiPre14 = fitdist(DataPre14, 'Poisson');
wPre14 = poiPre14.NLogL; %log likelihood (poisson dist.)
poiPost14 = fitdist(DataPost14, 'Poisson');
wPost14 = poiPost14.NLogL; %log likelihood (poisson dist.)

%Akaike Information Criterion
aicxPre14 = 2*2-2*log(xPre14); %AIC = 2*#parameters - 2ln(ll)
aicyPre14 = 2*2-2*log(yPre14);
aiczPre14 = 2*1-2*log(zPre14);
aicwPre14 = 2*1-2*log(wPre14);
aicPre14 = aicbic([xPre14, yPre14, zPre14, wPre14],[2, 2, 1, 1]); 
aicxPost14 = 2*2-2*log(xPost14); %AIC = 2*#parameters - 2ln(ll)
aicyPost14 = 2*2-2*log(yPost14);
aiczPost14 = 2*1-2*log(zPost14);
aicwPost14 = 2*1-2*log(wPost14);
aicPost14 = aicbic([xPost14, yPost14, zPost14, wPost14],[2, 2, 1, 1]);

%Anderson-Darling test
distnPre14 = makedist('normal','mu',1147.8,'sigma',3537.67);
[h,p] = adtest(DataPre14,'Distribution',distnPre14);
distnPost14 = makedist('normal','mu',64.625,'sigma',128.976);
[h,p] = adtest(DataPost14,'Distribution',distnPost14);

distnbPre14 = makedist('negative binomial','r',0.101468,'p',0.0000883946);
[h,p] = adtest(DataPre14,'Distribution',distnbPre14);
distnbPost14 = makedist('negative binomial','r',0.17425,'p',0.00268908);
[h,p] = adtest(DataPost14,'Distribution',distnbPost14);

distePre14 = makedist('exponential','mu',1147.8);
[h,p] = adtest(DataPre14,'Distribution',distePre14);
distePost14 = makedist('exponential','mu',64.625);
[h,p] = adtest(DataPost14,'Distribution',distePost14);

distpPre14 = makedist('poisson','lambda',1147.8);
[h,p] = adtest(DataPre14,'Distribution',distpPre14);
distpPost14 = makedist('poisson','lambda',64.625);
[h,p] = adtest(DataPost14,'Distribution',distpPost14);

%fitting distributions
figure;
r141 = subplot(2,2,1);
histfit(DataPre14,11,'Normal')
NormPre14 = mle(DataPre14, 'distribution', 'Normal');
xlabel('Parasite Count')
ylabel('Density')
title('Normal Distribution')
xlim([0,inf])
suptitle('Pre Mouse RAG 1 Experiment 4')

r142 = subplot(2,2,2);
histfit(DataPre14,11,'Negative Binomial')
NegBinPre14 = mle(DataPre14, 'distribution', 'Negative Binomial');
xlabel('Parasite Count')
ylabel('Density')
title('Negative Binomial Distribution')
axis([0,20700,0,40])

r143 = subplot(2,2,3);
histfit(DataPre14,11,'Exponential')
ExpoPre14 = mle(DataPre14, 'distribution', 'Exponential');
xlabel('Parasite Count')
ylabel('Density')
title('Exponential Distribution')

r144 = subplot(2,2,4);
histfit(DataPre14,11,'Poisson')
PoiPre14 = mle(DataPre14, 'distribution', 'Poisson');
xlabel('Parasite Count')
ylabel('Density')
title('Poisson Distribution')

linkaxes([r141,r143,r144],'xy')

figure;
R141 = subplot(2,2,1);
histfit(DataPost14,11,'Normal')
NormPost14 = mle(DataPost14, 'distribution', 'Normal');
xlabel('Parasite Count')
ylabel('Density')
title('Normal Distribution')
xlim([0,inf])
suptitle('Post Mouse RAG 1 Experiment 4')

R142 = subplot(2,2,2);
histfit(DataPost14,11,'Negative Binomial')
NegBinPost14 = mle(DataPost14, 'distribution', 'Negative Binomial');
xlabel('Parasite Count')
ylabel('Density')
title('Negative Binomial Distribution')
axis([0,644,0,30])

R143 = subplot(2,2,3);
histfit(DataPost14,11,'Exponential')
ExpoPost14 = mle(DataPost14, 'distribution', 'Exponential');
xlabel('Parasite Count')
ylabel('Density')
title('Exponential Distribution')

R144 = subplot(2,2,4);
histfit(DataPost14,11,'Poisson')
PoiPost14 = mle(DataPost14, 'distribution', 'Poisson');
xlabel('Parasite Count')
ylabel('Density')
title('Poisson Distribution')

linkaxes([R141,R143,R144],'xy')

%% Experimental ID 2 Experiment 4

Pre2 = [0	0.4	1.2	0.1	0	0	0	1	37.3	0.1	400.1	0	0	0	0	0	4395.1	0	0.1	0	0	0	0	0	0	6.2	0.1	0	0.3	3.5	0	0	0	0.1	0	0.6	0	0.1	0	3330.3];
Post2 = [3.5	1.3	0.4	99.3	0.1	2.9	3.3	2191.4	14.4	0.7	0.3	4.1	0.1	0.2	330.7	2.6	1.8	0.3	0.2	0.1	0.7	0.2	0.2	0	3.7	1.7	9	7.2	1	1.6	0.2	1	3055.3	522.4	0.2	0.4	1.6	25559.2	0.4	0.4];

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
distnPre2 = makedist('normal','mu',204.375,'sigma',860.763);
[h,p] = adtest(DataPre2,'Distribution',distnPre2);
distnPost2 = makedist('normal','mu',795.525,'sigma',4058.78);
[h,p] = adtest(DataPost2,'Distribution',distnPost2);

distnbPre2 = makedist('negative binomial','r',0.0277886,'p',0.00013595);
[h,p] = adtest(DataPre2,'Distribution',distnbPre2);
distnbPost2 = makedist('negative binomial','r',0.0794412,'p',0.0000998501);
[h,p] = adtest(DataPost2,'Distribution',distnbPost2);

distePre2 = makedist('exponential','mu',204.375);
[h,p] = adtest(DataPre2,'Distribution',distePre2);
distePost2 = makedist('exponential','mu',795.525);
[h,p] = adtest(DataPost2,'Distribution',distePost2);

distpPre2 = makedist('poisson','lambda',204.375);
[h,p] = adtest(DataPre2,'Distribution',distpPre2);
distpPost2 = makedist('poisson','lambda',795.525);
[h,p] = adtest(DataPost2,'Distribution',distpPost2);

%fitting distributions
figure;
r21 = subplot(2,2,1);
histfit(DataPre2,11,'Normal')
NormPre2 = mle(DataPre2, 'distribution', 'Normal');
xlabel('Parasite Count')
ylabel('Density')
title('Normal Distribution')
xlim([0,inf])
suptitle('Pre Mouse RAG 2 Experiment 4')

r22 = subplot(2,2,2);
histfit(DataPre2,11,'Negative Binomial')
NegBinPre2 = mle(DataPre2, 'distribution', 'Negative Binomial');
xlabel('Parasite Count')
ylabel('Density')
title('Negative Binomial Distribution')
axis([0,4600,0,40])

r23 = subplot(2,2,3);
histfit(DataPre2,11,'Exponential')
ExpoPre2 = mle(DataPre2, 'distribution', 'Exponential');
xlabel('Parasite Count')
ylabel('Density')
title('Exponential Distribution')

r24 = subplot(2,2,4);
histfit(DataPre2,11,'Poisson')
PoiPre2 = mle(DataPre2, 'distribution', 'Poisson');
xlabel('Parasite Count')
ylabel('Density')
title('Poisson Distribution')

linkaxes([r21,r23,r24],'xy')

figure;
R21 = subplot(2,2,1);
histfit(DataPost2,11,'Normal')
NormPost2 = mle(DataPost2, 'distribution', 'Normal');
xlabel('Parasite Count')
ylabel('Density')
title('Normal Distribution')
xlim([0,inf])
suptitle('Post Mouse RAG 2 Experiment 4')

R22 = subplot(2,2,2);
histfit(DataPost2,11,'Negative Binomial')
NegBinPost2 = mle(DataPost2, 'distribution', 'Negative Binomial');
xlabel('Parasite Count')
ylabel('Density')
title('Negative Binomial Distribution')
axis([0,27600,0,40])

R23 = subplot(2,2,3);
histfit(DataPost2,11,'Exponential')
ExpoPost2 = mle(DataPost2, 'distribution', 'Exponential');
xlabel('Parasite Count')
ylabel('Density')
title('Exponential Distribution')

R24 = subplot(2,2,4);
histfit(DataPost2,11,'Poisson')
PoiPost2 = mle(DataPost2, 'distribution', 'Poisson');
xlabel('Parasite Count')
ylabel('Density')
title('Poisson Distribution')

linkaxes([R21,R23,R24],'xy')

%% Experimental ID 6 Experiment 4

Pre64 = [7477.4	14329	6075.6	1.7	187.1	7487.5	20163.4	0.4	26311.9	9.9	1.4	17999.6	1285.4	107.9	15895.1	15150.7	1.7	2449	18093.9	17692.7	12182.9	19800.4	6272.9	22233.9	8125.6	9961.6	22065.3	22622.6	3286.7	27884.8	63332.6	13730.3	9791.7	678	10504.9	2622.6	3201.1	1.6	8167.8	840.4];
Post64 = [3546.9	62531.9	10708.5	20164.9	15823.2	30946.8	43496	673.7	2359.5	16694	2.1	2061.7	0.7	2907.7	31219.4	38228.3	60933.9	3.9	3.5	3798.3	4	42926.6	42	79377.3	9.3	36725.8	33753.8	55109.4	18997.8	26817.2	37516.3	7567.2	55.9	33335.7	15379.3	34329	25426.8	35.1	77003.3	18.2];

%data
DataPre64 = round(reshape(Pre64, [numel(Pre64), 1]));
DataPost64 = round(reshape(Post64, [numel(Post64), 1]));

%log likelihood
norPre64 = fitdist(DataPre64, 'Normal');
xPre64 = norPre64.NLogL; %log likelihood (normal dist.)
norPost64 = fitdist(DataPost64, 'Normal');
xPost64 = norPost64.NLogL; %log likelihood (normal dist.)

nbPre64 = fitdist(DataPre64, 'Negative Binomial');
yPre64 = nbPre64.NLogL; %log likelihood (Neg. Bin. dist.)
nbPost64 = fitdist(DataPost64, 'Negative Binomial');
yPost64 = nbPost64.NLogL; %log likelihood (Neg. Bin. dist.)

expPre64 = fitdist(DataPre64, 'Exponential');
zPre64 = expPre64.NLogL; %log likelihood (exponential dist.)
expPost64 = fitdist(DataPost64, 'Exponential');
zPost64 = expPost64.NLogL; %log likelihood (exponential dist.)

poiPre64 = fitdist(DataPre64, 'Poisson');
wPre64 = poiPre64.NLogL; %log likelihood (poisson dist.)
poiPost64 = fitdist(DataPost64, 'Poisson');
wPost64 = poiPost64.NLogL; %log likelihood (poisson dist.)

%Akaike Information Criterion
aicxPre64 = 2*2-2*log(xPre64); %AIC = 2*#parameters - 2ln(ll)
aicyPre64 = 2*2-2*log(yPre64);
aiczPre64 = 2*1-2*log(zPre64);
aicwPre64 = 2*1-2*log(wPre64);
aicPre64 = aicbic([xPre64, yPre64, zPre64, wPre64],[2, 2, 1, 1]); 
aicxPost64 = 2*2-2*log(xPost64); %AIC = 2*#parameters - 2ln(ll)
aicyPost64 = 2*2-2*log(yPost64);
aiczPost64 = 2*1-2*log(zPost64);
aicwPost64 = 2*1-2*log(wPost64);
aicPost64 = aicbic([xPost64, yPost64, zPost64, wPost64],[2, 2, 1, 1]);

%Anderson-Darling test
distnPre64 = makedist('normal','mu',10950.8,'sigma',12012.7);
[h,p] = adtest(DataPre64,'Distribution',distnPre64);
distnPost64 = makedist('normal','mu',21763.4,'sigma',22856.7);
[h,p] = adtest(DataPost64,'Distribution',distnPost64);

distnbPre64 = makedist('negative binomial','r',0.384226,'p',0.0000350854);
[h,p] = adtest(DataPre64,'Distribution',distnbPre64);
distnbPost64 = makedist('negative binomial','r',0.322317,'p',0.0000148098);
[h,p] = adtest(DataPost64,'Distribution',distnbPost64);

distePre64 = makedist('exponential','mu',10950.8);
[h,p] = adtest(DataPre64,'Distribution',distePre64);
distePost64 = makedist('exponential','mu',21763.4);
[h,p] = adtest(DataPost64,'Distribution',distePost64);

distpPre64 = makedist('poisson','lambda',10950.8);
[h,p] = adtest(DataPre64,'Distribution',distpPre64);
distpPost64 = makedist('poisson','lambda',21763.4);
[h,p] = adtest(DataPost64,'Distribution',distpPost64);

%fitting distributions
figure;
r641 = subplot(2,2,1);
histfit(DataPre64,11,'Normal')
NormPre64 = mle(DataPre64, 'distribution', 'Normal');
xlabel('Parasite Count')
ylabel('Density')
title('Normal Distribution')
axis([0,66700,0,15])
suptitle('Pre Mouse RAG 6 Experiment 4')

r642 = subplot(2,2,2);
histfit(DataPre64,11,'Negative Binomial')
NegBinPre64 = mle(DataPre64, 'distribution', 'Negative Binomial');
xlabel('Parasite Count')
ylabel('Density')
title('Negative Binomial Distribution')
axis([0,66700,0,15])

r643 = subplot(2,2,3);
histfit(DataPre64,11,'Exponential')
ExpoPre64 = mle(DataPre64, 'distribution', 'Exponential');
xlabel('Parasite Count')
ylabel('Density')
title('Exponential Distribution')
axis([0,66700,0,15])

r644 = subplot(2,2,4);
histfit(DataPre64,11,'Poisson')
PoiPre64 = mle(DataPre64, 'distribution', 'Poisson');
xlabel('Parasite Count')
ylabel('Density')
title('Poisson Distribution')
axis([0,66700,0,15])

figure;
R641 = subplot(2,2,1);
histfit(DataPost64,11,'Normal')
NormPost64 = mle(DataPost64, 'distribution', 'Normal');
xlabel('Parasite Count')
ylabel('Density')
title('Normal Distribution')
axis([0,76650,0,20])
suptitle('Post Mouse RAG 6 Experiment 4')

R642 = subplot(2,2,2);
histfit(DataPost64,11,'Negative Binomial')
NegBinPost64 = mle(DataPost64, 'distribution', 'Negative Binomial');
xlabel('Parasite Count')
ylabel('Density')
title('Negative Binomial Distribution')
axis([0,76650,0,20])

R643 = subplot(2,2,3);
histfit(DataPost64,11,'Exponential')
ExpoPost64 = mle(DataPost64, 'distribution', 'Exponential');
xlabel('Parasite Count')
ylabel('Density')
title('Exponential Distribution')
axis([0,76650,0,20])

R644 = subplot(2,2,4);
histfit(DataPost64,11,'Poisson')
PoiPost64 = mle(DataPost64, 'distribution', 'Poisson');
xlabel('Parasite Count')
ylabel('Density')
title('Poisson Distribution')
axis([0,76650,0,20])

%% Experimental ID 8 Experiment 4

Pre8 = [685.9	192.5	13783.4	4162.7	8695.3	19029.3	48.2	10060.3	3123.3	6614.5	188.6	78.1	228.3	0	14.6	0.9	9008.1	2183.3	8739.7	2184.3	341	852.5	2.5	26812	16293.6	3410.1	16014.4	13022.2	1950.5	20.6	1508.5	1780.1	46395.4	7452.9	35	235.8	245.9	7176.8	7467.9	1558.5];
Post8 = [88906	19149.5	18599.4	9562.3	14940.8	36613.4	20706.3	9008.5	10162.4	34724.1	0.5	12324.1	31943.8	8432.6	22005.8	34875.9	5.2	3496.7	1021.9	30555.2	3.4	9530.1	175.9	10234.8	15064.4	2524.2	4107.9	5326.5	37661.2	10965.5	485.7	24218.8	5724.7	27079.4	3414.6	10154.8	84.5	17578.7	16548.6	38973.2];

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
distnPre8 = makedist('normal','mu',6040,'sigma',9140.11);
[h,p] = adtest(DataPre8,'Distribution',distnPre8);
distnPost8 = makedist('normal','mu',16172.4,'sigma',16872.3);
[h,p] = adtest(DataPost8,'Distribution',distnPost8);

distnbPre8 = makedist('negative binomial','r',0.353602,'p',0.0000585399);
[h,p] = adtest(DataPre8,'Distribution',distnbPre8);
distnbPost8 = makedist('negative binomial','r',0.542371,'p',0.0000335358);
[h,p] = adtest(DataPost8,'Distribution',distnbPost8);

distePre8 = makedist('exponential','mu',6040);
[h,p] = adtest(DataPre8,'Distribution',distePre8);
distePost8 = makedist('exponential','mu',16172.4);
[h,p] = adtest(DataPost8,'Distribution',distePost8);

distpPre8 = makedist('poisson','lambda',6040);
[h,p] = adtest(DataPre8,'Distribution',distpPre8);
distpPost8 = makedist('poisson','lambda',16172.4);
[h,p] = adtest(DataPost8,'Distribution',distpPost8);

%fitting distributions
figure;
r81 = subplot(2,2,1);
histfit(DataPre8,11,'Normal')
NormPre8 = mle(DataPre8, 'distribution', 'Normal');
xlabel('Parasite Count')
ylabel('Density')
title('Normal Distribution')
xlim([0,inf])
suptitle('Pre Mouse RAG 8 Experiment 4')

r82 = subplot(2,2,2);
histfit(DataPre8,11,'Negative Binomial')
NegBinPre8 = mle(DataPre8, 'distribution', 'Negative Binomial');
xlabel('Parasite Count')
ylabel('Density')
title('Negative Binomial Distribution')
axis([0,49450,0,30])

r83 = subplot(2,2,3);
histfit(DataPre8,11,'Exponential')
ExpoPre8 = mle(DataPre8, 'distribution', 'Exponential');
xlabel('Parasite Count')
ylabel('Density')
title('Exponential Distribution')

r84 = subplot(2,2,4);
histfit(DataPre8,11,'Poisson')
PoiPre8 = mle(DataPre8, 'distribution', 'Poisson');
xlabel('Parasite Count')
ylabel('Density')
title('Poisson Distribution')

linkaxes([r81,r83,r84],'xy')

figure;
R81 = subplot(2,2,1);
histfit(DataPost8,11,'Normal')
NormPost8 = mle(DataPost8, 'distribution', 'Normal');
xlabel('Parasite Count')
ylabel('Density')
title('Normal Distribution')
xlim([0,inf])
suptitle('Post Mouse RAG 8 Experiment 4')

R82 = subplot(2,2,2);
histfit(DataPost8,11,'Negative Binomial')
NegBinPost8 = mle(DataPost8, 'distribution', 'Negative Binomial');
xlabel('Parasite Count')
ylabel('Density')
title('Negative Binomial Distribution')
axis([0,93150,0,15])

R83 = subplot(2,2,3);
histfit(DataPost8,11,'Exponential')
ExpoPost8 = mle(DataPost8, 'distribution', 'Exponential');
xlabel('Parasite Count')
ylabel('Density')
title('Exponential Distribution')
axis([0,93150,0,15])

R84 = subplot(2,2,4);
histfit(DataPost8,11,'Poisson')
PoiPost8 = mle(DataPost8, 'distribution', 'Poisson');
xlabel('Parasite Count')
ylabel('Density')
title('Poisson Distribution')

linkaxes([R81,R84],'xy')

%% Experimental ID 1 Experiment 5

Pre15 = [676.4	4.5	1.8	3.8	615.4	131.6	1.5	30.6	1.8	1	1.5	5.7	113.9	6.3	171.3	9	239	0	5.3	0.1	0.9	515.6	4.8	8.6	63.6	136	0.6	36.1	8.3	1.3	0.2	4.9	6.1	30.2	82.4	27.8	25.1	0.2	50	127.1];
Post15 = [1.4	17	106.7	1935.1	1440.5	4.2	0.9	0.5	0.9	1.5	2.8	0.2	4506	0.1	699	3	0.3	0.1	29.8	0	14	12.5	0.1	0.6	0.4	851.1	77.1	28.9	2280.2	44.5	545.4	0.3	0.1	6.1	2.5	190.8	0	0	0	0.2];

%data
DataPre15 = round(reshape(Pre15, [numel(Pre15), 1]));
DataPost15 = round(reshape(Post15, [numel(Post15), 1]));

%log likelihood
norPre15 = fitdist(DataPre15, 'Normal');
xPre15 = norPre15.NLogL; %log likelihood (normal dist.)
norPost15 = fitdist(DataPost15, 'Normal');
xPost15 = norPost15.NLogL; %log likelihood (normal dist.)

nbPre15 = fitdist(DataPre15, 'Negative Binomial');
yPre15 = nbPre15.NLogL; %log likelihood (Neg. Bin. dist.)
nbPost15 = fitdist(DataPost15, 'Negative Binomial');
yPost15 = nbPost15.NLogL; %log likelihood (Neg. Bin. dist.)

expPre15 = fitdist(DataPre15, 'Exponential');
zPre15 = expPre15.NLogL; %log likelihood (exponential dist.)
expPost15 = fitdist(DataPost15, 'Exponential');
zPost15 = expPost15.NLogL; %log likelihood (exponential dist.)

poiPre15 = fitdist(DataPre15, 'Poisson');
wPre15 = poiPre15.NLogL; %log likelihood (poisson dist.)
poiPost15 = fitdist(DataPost15, 'Poisson');
wPost15 = poiPost15.NLogL; %log likelihood (poisson dist.)

%Akaike Information Criterion
aicxPre15 = 2*2-2*log(xPre15); %AIC = 2*#parameters - 2ln(ll)
aicyPre15 = 2*2-2*log(yPre15);
aiczPre15 = 2*1-2*log(zPre15);
aicwPre15 = 2*1-2*log(wPre15);
aicPre15 = aicbic([xPre15, yPre15, zPre15, wPre15],[2, 2, 1, 1]); 
aicxPost15 = 2*2-2*log(xPost15); %AIC = 2*#parameters - 2ln(ll)
aicyPost15 = 2*2-2*log(yPost15);
aiczPost15 = 2*1-2*log(zPost15);
aicwPost15 = 2*1-2*log(wPost15);
aicPost15 = aicbic([xPost15, yPost15, zPost15, wPost15],[2, 2, 1, 1]);

%Anderson-Darling test
distnPre15 = makedist('normal','mu',78.8,'sigma',161.947);
[h,p] = adtest(DataPre15,'Distribution',distnPre15);
distnPost15 = makedist('normal','mu',320.15,'sigma',859.211);
[h,p] = adtest(DataPost15,'Distribution',distnPost15);

distnbPre15 = makedist('negative binomial','r',0.309354,'p',0.00391046);
[h,p] = adtest(DataPre15,'Distribution',distnbPre15);
distnbPost15 = makedist('negative binomial','r',0.125262,'p',0.000391106);
[h,p] = adtest(DataPost15,'Distribution',distnbPost15);

distePre15 = makedist('exponential','mu',6040);
[h,p] = adtest(DataPre15,'Distribution',distePre15);
distePost15 = makedist('exponential','mu',320.15);
[h,p] = adtest(DataPost15,'Distribution',distePost15);

distpPre15 = makedist('poisson','lambda',78.8);
[h,p] = adtest(DataPre15,'Distribution',distpPre15);
distpPost15 = makedist('poisson','lambda',320.15);
[h,p] = adtest(DataPost15,'Distribution',distpPost15);

%fitting distributions
figure;
r151 = subplot(2,2,1);
histfit(DataPre15,11,'Normal')
NormPre15 = mle(DataPre15, 'distribution', 'Normal');
xlabel('Parasite Count')
ylabel('Density')
title('Normal Distribution')
xlim([0,inf])
suptitle('Pre Mouse RAG 1 Experiment 5')

r152 = subplot(2,2,2);
histfit(DataPre15,11,'Negative Binomial')
NegBinPre15 = mle(DataPre15, 'distribution', 'Negative Binomial');
xlabel('Parasite Count')
ylabel('Density')
title('Negative Binomial Distribution')
axis([0,713,0,30])

r153 = subplot(2,2,3);
histfit(DataPre15,11,'Exponential')
ExpoPre15 = mle(DataPre15, 'distribution', 'Exponential');
xlabel('Parasite Count')
ylabel('Density')
title('Exponential Distribution')

r154 = subplot(2,2,4);
histfit(DataPre15,11,'Poisson')
PoiPre15 = mle(DataPre15, 'distribution', 'Poisson');
xlabel('Parasite Count')
ylabel('Density')
title('Poisson Distribution')

linkaxes([r151,r153,r154],'xy')

figure;
R151 = subplot(2,2,1);
histfit(DataPost15,11,'Normal')
NormPost15 = mle(DataPost15, 'distribution', 'Normal');
xlabel('Parasite Count')
ylabel('Density')
title('Normal Distribution')
xlim([0,inf])
suptitle('Post Mouse RAG 1 Experiment 5')

R152 = subplot(2,2,2);
histfit(DataPost15,11,'Negative Binomial')
NegBinPost15 = mle(DataPost15, 'distribution', 'Negative Binomial');
xlabel('Parasite Count')
ylabel('Density')
title('Negative Binomial Distribution')
axis([0,4715,0,40])

R153 = subplot(2,2,3);
histfit(DataPost15,11,'Exponential')
ExpoPost15 = mle(DataPost15, 'distribution', 'Exponential');
xlabel('Parasite Count')
ylabel('Density')
title('Exponential Distribution')

R154 = subplot(2,2,4);
histfit(DataPost15,11,'Poisson')
PoiPost15 = mle(DataPost15, 'distribution', 'Poisson');
xlabel('Parasite Count')
ylabel('Density')
title('Poisson Distribution')

linkaxes([R151,R153,R154],'xy')

%% Experimental ID 6 Experiment 5

Pre65 = [222.2	16.7	10.2	26.5	16.4	14.7	440.4	23.2	6.1	3.9	1.7	0.3	0.7	0.1	0.6	0.8	5.1	14.9	0	0.3	0.4	5.1	1.1	0.3	0.2	0.6	0.1	4.6	0.1	0	0	0.1	0	1	0	0	14.5	0.2	0.1	1.1];
Post65 = [616	381.5	0	0.3	1.1	1	0	26.8	0.1	0.2	0	1.1	0.1	0	0.4	0	0	294	0	0	1.9	0	0	0.6	0	57.3	0	0.3	0	0	0	0.3	3145.1	31.5	0	0	0	0	0.9	0];

%data
DataPre65 = round(reshape(Pre65, [numel(Pre65), 1]));
DataPost65 = round(reshape(Post65, [numel(Post65), 1]));

%log likelihood
norPre65 = fitdist(DataPre65, 'Normal');
xPre65 = norPre65.NLogL; %log likelihood (normal dist.)
norPost65 = fitdist(DataPost65, 'Normal');
xPost65 = norPost65.NLogL; %log likelihood (normal dist.)

nbPre65 = fitdist(DataPre65, 'Negative Binomial');
yPre65 = nbPre65.NLogL; %log likelihood (Neg. Bin. dist.)
nbPost65 = fitdist(DataPost65, 'Negative Binomial');
yPost65 = nbPost65.NLogL; %log likelihood (Neg. Bin. dist.)

expPre65 = fitdist(DataPre65, 'Exponential');
zPre65 = expPre65.NLogL; %log likelihood (exponential dist.)
expPost65 = fitdist(DataPost65, 'Exponential');
zPost65 = expPost65.NLogL; %log likelihood (exponential dist.)

poiPre65 = fitdist(DataPre65, 'Poisson');
wPre65 = poiPre65.NLogL; %log likelihood (poisson dist.)
poiPost65 = fitdist(DataPost65, 'Poisson');
wPost65 = poiPost65.NLogL; %log likelihood (poisson dist.)

%Akaike Information Criterion
aicxPre65 = 2*2-2*log(xPre65); %AIC = 2*#parameters - 2ln(ll)
aicyPre65 = 2*2-2*log(yPre65);
aiczPre65 = 2*1-2*log(zPre65);
aicwPre65 = 2*1-2*log(wPre65);
aicPre65 = aicbic([xPre65, yPre65, zPre65, wPre65],[2, 2, 1, 1]); 
aicxPost65 = 2*2-2*log(xPost65); %AIC = 2*#parameters - 2ln(ll)
aicyPost65 = 2*2-2*log(yPost65);
aiczPost65 = 2*1-2*log(zPost65);
aicwPost65 = 2*1-2*log(wPost65);
aicPost65 = aicbic([xPost65, yPost65, zPost65, wPost65],[2, 2, 1, 1]);

%Anderson-Darling test
distnPre65 = makedist('normal','mu',20.85,'sigma',76.4955);
[h,p] = adtest(DataPre65,'Distribution',distnPre65);
distnPost65 = makedist('normal','mu',114,'sigma',506.066);
[h,p] = adtest(DataPost65,'Distribution',distnPost65);

distnbPre65 = makedist('negative binomial','r',0.154355,'p',0.00734872);
[h,p] = adtest(DataPre65,'Distribution',distnbPre65);
distnbPost65 = makedist('negative binomial','r',0.0482573,'p',0.000423131);
[h,p] = adtest(DataPost65,'Distribution',distnbPost65);

distePre65 = makedist('exponential','mu',20.85);
[h,p] = adtest(DataPre65,'Distribution',distePre65);
distePost65 = makedist('exponential','mu',114);
[h,p] = adtest(DataPost65,'Distribution',distePost65);

distpPre65 = makedist('poisson','lambda',20.85);
[h,p] = adtest(DataPre65,'Distribution',distpPre65);
distpPost65 = makedist('poisson','lambda',114);
[h,p] = adtest(DataPost65,'Distribution',distpPost65);

%fitting distributions
figure;
r651 = subplot(2,2,1);
histfit(DataPre65,11,'Normal')
NormPre65 = mle(DataPre65, 'distribution', 'Normal');
xlabel('Parasite Count')
ylabel('Density')
title('Normal Distribution')
xlim([0,inf])
suptitle('Pre Mouse RAG 6 Experiment 5')

r652 = subplot(2,2,2);
histfit(DataPre65,11,'Negative Binomial')
NegBinPre65 = mle(DataPre65, 'distribution', 'Negative Binomial');
xlabel('Parasite Count')
ylabel('Density')
title('Negative Binomial Distribution')
axis([0,440,0,40])

r653 = subplot(2,2,3);
histfit(DataPre65,11,'Exponential')
ExpoPre65 = mle(DataPre65, 'distribution', 'Exponential');
xlabel('Parasite Count')
ylabel('Density')
title('Exponential Distribution')

r654 = subplot(2,2,4);
histfit(DataPre65,11,'Poisson')
PoiPre65 = mle(DataPre65, 'distribution', 'Poisson');
xlabel('Parasite Count')
ylabel('Density')
title('Poisson Distribution')

linkaxes([r651,r653,r654],'xy')

figure;
R651 = subplot(2,2,1);
histfit(DataPost65,11,'Normal')
NormPost65 = mle(DataPost65, 'distribution', 'Normal');
xlabel('Parasite Count')
ylabel('Density')
title('Normal Distribution')
xlim([0,inf])
suptitle('Post Mouse RAG 6 Experiment 5')

R652 = subplot(2,2,2);
histfit(DataPost65,11,'Negative Binomial')
NegBinPost65 = mle(DataPost65, 'distribution', 'Negative Binomial');
xlabel('Parasite Count')
ylabel('Density')
title('Negative Binomial Distribution')
axis([0,3190,0,40])

R653 = subplot(2,2,3);
histfit(DataPost65,11,'Exponential')
ExpoPost65 = mle(DataPost65, 'distribution', 'Exponential');
xlabel('Parasite Count')
ylabel('Density')
title('Exponential Distribution')

R654 = subplot(2,2,4);
histfit(DataPost65,11,'Poisson')
PoiPost65 = mle(DataPost65, 'distribution', 'Poisson');
xlabel('Parasite Count')
ylabel('Density')
title('Poisson Distribution')

linkaxes([R651,R653,R654],'xy')

%% Experimental ID 10 Experiment 5

Pre10 = [18.3	11.1	134.4	8.9	6.8	16.6	23.2	12.4	5.9	74.7	15.5	19	13.8	9.7	18.2	5.7	13	3.6	26	11.2	6.4	6	8.1	7.8	32.2	48	29.2	32.5	90	6.5	3.8	4.9	4.3	20.3	6.7	7.3	12.4	10.9	4.5	9.2];
Post10 = [0.1	0.3	5.9	2.9	0.1	2.4	10.8	0.1	19.3	0.8	1.4	0.1	0	0	89.2	0.1	1.6	0.1	0.1	0.1	0.1	0	23.7	0	0.1	0.3	0.4	0.1	0	0	0	0.1	5161.1	0.2	0.2	0	0	0.1	1124.5	0.3];

%data
DataPre10 = round(reshape(Pre10, [numel(Pre10), 1]));
DataPost10 = round(reshape(Post10, [numel(Post10), 1]));

%log likelihood
norPre10 = fitdist(DataPre10, 'Normal');
xPre10 = norPre10.NLogL; %log likelihood (normal dist.)
norPost10 = fitdist(DataPost10, 'Normal');
xPost10 = norPost10.NLogL; %log likelihood (normal dist.)

nbPre10 = fitdist(DataPre10, 'Negative Binomial');
yPre10 = nbPre10.NLogL; %log likelihood (Neg. Bin. dist.)
nbPost10 = fitdist(DataPost10, 'Negative Binomial');
yPost10 = nbPost10.NLogL; %log likelihood (Neg. Bin. dist.)

expPre10 = fitdist(DataPre10, 'Exponential');
zPre10 = expPre10.NLogL; %log likelihood (exponential dist.)
expPost10 = fitdist(DataPost10, 'Exponential');
zPost10 = expPost10.NLogL; %log likelihood (exponential dist.)

poiPre10 = fitdist(DataPre10, 'Poisson');
wPre10 = poiPre10.NLogL; %log likelihood (poisson dist.)
poiPost10 = fitdist(DataPost10, 'Poisson');
wPost10 = poiPost10.NLogL; %log likelihood (poisson dist.)

%Akaike Information Criterion
aicxPre10 = 2*2-2*log(xPre10); %AIC = 2*#parameters - 2ln(ll)
aicyPre10 = 2*2-2*log(yPre10);
aiczPre10 = 2*1-2*log(zPre10);
aicwPre10 = 2*1-2*log(wPre10);
aicPre10 = aicbic([xPre10, yPre10, zPre10, wPre10],[2, 2, 1, 1]); 
aicxPost10 = 2*2-2*log(xPost10); %AIC = 2*#parameters - 2ln(ll)
aicyPost10 = 2*2-2*log(yPost10);
aiczPost10 = 2*1-2*log(zPost10);
aicwPost10 = 2*1-2*log(wPost10);
aicPost10 = aicbic([xPost10, yPost10, zPost10, wPost10],[2, 2, 1, 1]);

%Anderson-Darling test
distnPre10 = makedist('normal','mu',20.0256,'sigma',26.1277);
[h,p] = adtest(DataPre10,'Distribution',distnPre10);
distnPost10 = makedist('normal','mu',161.1,'sigma',830.083);
[h,p] = adtest(DataPost10,'Distribution',distnPost10);

distnbPre10 = makedist('negative binomial','r',1.28653,'p',0.0603658);
[h,p] = adtest(DataPre10,'Distribution',distnbPre10);
distnbPost10 = makedist('negative binomial','r',0.0408038,'p',0.000253219);
[h,p] = adtest(DataPost10,'Distribution',distnbPost10);

distePre10 = makedist('exponential','mu',20.0256);
[h,p] = adtest(DataPre10,'Distribution',distePre10);
distePost10 = makedist('exponential','mu',161.1);
[h,p] = adtest(DataPost10,'Distribution',distePost10);

distpPre10 = makedist('poisson','lambda',20.0256);
[h,p] = adtest(DataPre10,'Distribution',distpPre10);
distpPost10 = makedist('poisson','lambda',161.1);
[h,p] = adtest(DataPost10,'Distribution',distpPost10);

%fitting distributions
figure;
r101 = subplot(2,2,1);
histfit(DataPre10,11,'Normal')
NormPre10 = mle(DataPre10, 'distribution', 'Normal');
xlabel('Parasite Count')
ylabel('Density')
title('Normal Distribution')
xlim([0,inf])
suptitle('Pre Mouse RAG 10 Experiment 5')

r102 = subplot(2,2,2);
histfit(DataPre10,11,'Negative Binomial')
NegBinPre10 = mle(DataPre10, 'distribution', 'Negative Binomial');
xlabel('Parasite Count')
ylabel('Density')
title('Negative Binomial Distribution')
axis([0,140,0,30])

r103 = subplot(2,2,3);
histfit(DataPre10,11,'Exponential')
ExpoPre10 = mle(DataPre10, 'distribution', 'Exponential');
xlabel('Parasite Count')
ylabel('Density')
title('Exponential Distribution')

r104 = subplot(2,2,4);
histfit(DataPre10,11,'Poisson')
PoiPre10 = mle(DataPre10, 'distribution', 'Poisson');
xlabel('Parasite Count')
ylabel('Density')
title('Poisson Distribution')

linkaxes([r101,r103,r104],'xy')

figure;
R101 = subplot(2,2,1);
histfit(DataPost10,11,'Normal')
NormPost10 = mle(DataPost10, 'distribution', 'Normal');
xlabel('Parasite Count')
ylabel('Density')
title('Normal Distribution')
xlim([0,inf])
suptitle('Post Mouse RAG 10 Experiment 5')

R102 = subplot(2,2,2);
histfit(DataPost10,11,'Negative Binomial')
NegBinPost10 = mle(DataPost10, 'distribution', 'Negative Binomial');
xlabel('Parasite Count')
ylabel('Density')
title('Negative Binomial Distribution')
axis([0,5170,0,40])

R103 = subplot(2,2,3);
histfit(DataPost10,11,'Exponential')
ExpoPost10 = mle(DataPost10, 'distribution', 'Exponential');
xlabel('Parasite Count')
ylabel('Density')
title('Exponential Distribution')

R104 = subplot(2,2,4);
histfit(DataPost10,11,'Poisson')
PoiPost10 = mle(DataPost10, 'distribution', 'Poisson');
xlabel('Parasite Count')
ylabel('Density')
title('Poisson Distribution')

linkaxes([R101,R103,R104],'xy')
