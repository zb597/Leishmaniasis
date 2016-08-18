%RPMI mice

%% Experimental ID F Experiment 2

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

%% Experimental ID J Experiment 2

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

%% Experimental ID 3 Experiment 4

Pre34 = [0	0.6	0	0.4	0	0	0	5.8	0	0	0	0	0	0	0	0	0	0.1	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	36.4	0	0];
Post34 = [1	0.4	25.5	0.3	33.9	0.9	0.4	0.3	14.6	0.1	0.1	0.3	0	0.8	0.1	4460.3	0.8	0.1	0.1	0.7	0.2	0.1	0.1	70.1	0.1	0.5	528.8	9.1	0.2	0.4	0	0	0.3	1.5	0	6.5	0.1	1.5	0.3	0.2];

%data
DataPre34 = round(reshape(Pre34, [numel(Pre34), 1]));
DataPost34 = round(reshape(Post34, [numel(Post34), 1]));

%log likelihood
norPre34 = fitdist(DataPre34, 'Normal');
xPre34 = norPre34.NLogL; %log likelihood (normal dist.)
norPost34 = fitdist(DataPost34, 'Normal');
xPost34 = norPost34.NLogL; %log likelihood (normal dist.)

nbPre34 = fitdist(DataPre34, 'Negative Binomial');
yPre34 = nbPre34.NLogL; %log likelihood (Neg. Bin. dist.)
nbPost34 = fitdist(DataPost34, 'Negative Binomial');
yPost34 = nbPost34.NLogL; %log likelihood (Neg. Bin. dist.)

expPre34 = fitdist(DataPre34, 'Exponential');
zPre34 = expPre34.NLogL; %log likelihood (exponential dist.)
expPost34 = fitdist(DataPost34, 'Exponential');
zPost34 = expPost34.NLogL; %log likelihood (exponential dist.)

poiPre34 = fitdist(DataPre34, 'Poisson');
wPre34 = poiPre34.NLogL; %log likelihood (poisson dist.)
poiPost34 = fitdist(DataPost34, 'Poisson');
wPost34 = poiPost34.NLogL; %log likelihood (poisson dist.)

%Akaike Information Criterion
aicxPre34 = 2*2-2*log(xPre34); %AIC = 2*#parameters - 2ln(ll)
aicyPre34 = 2*2-2*log(yPre34);
aiczPre34 = 2*1-2*log(zPre34);
aicwPre34 = 2*1-2*log(wPre34);
aicPre34 = aicbic([xPre34, yPre34, zPre34, wPre34],[2, 2, 1, 1]); 
aicxPost34 = 2*2-2*log(xPost34); %AIC = 2*#parameters - 2ln(ll)
aicyPost34 = 2*2-2*log(yPost34);
aiczPost34 = 2*1-2*log(zPost34);
aicwPost34 = 2*1-2*log(wPost34);
aicPost34 = aicbic([xPost34, yPost34, zPost34, wPost34],[2, 2, 1, 1]);

%Anderson-Darling test
distnPre34 = makedist('normal','mu',1.075,'sigma',5.74406);
[h,p] = adtest(DataPre34,'Distribution',distnPre34);
distnPost34 = makedist('normal','mu',125,'sigma',707.341);
[h,p] = adtest(DataPost34,'Distribution',distnPost34);

distnbPre34 = makedist('negative binomial','r',0.0193504,'p',0.0176821);
[h,p] = adtest(DataPre34,'Distribution',distnbPre34);
distnbPost34 = makedist('negative binomial','r',0.0594375,'p',0.000460544);
[h,p] = adtest(DataPost34,'Distribution',distnbPost34);

distePre34 = makedist('exponential','mu',1.075);
[h,p] = adtest(DataPre34,'Distribution',distePre34);
distePost34 = makedist('exponential','mu',129);
[h,p] = adtest(DataPost34,'Distribution',distePost34);

distpPre34 = makedist('poisson','lambda',1.075);
[h,p] = adtest(DataPre34,'Distribution',distpPre34);
distpPost34 = makedist('poisson','lambda',129);
[h,p] = adtest(DataPost34,'Distribution',distpPost34);

%fitting distributions
figure;
r341 = subplot(2,2,1);
norPre34 = fitdist(DataPre34, 'Normal');
histfit(DataPre34,11,'Normal')
NormPre34 = mle(DataPre34, 'distribution', 'Normal');
xlabel('Parasite Count')
ylabel('Density')
title('Normal Distribution')
xlim([0,inf])
suptitle('Pre Mouse RAG 3 Experiment 4')

r342 = subplot(2,2,2);
nbPre34 = fitdist(DataPre34, 'Negative Binomial');
histfit(DataPre34,11,'Negative Binomial')
NegBinPre34 = mle(DataPre34, 'distribution', 'Negative Binomial');
xlabel('Parasite Count')
ylabel('Density')
title('Negative Binomial Distribution')
axis([0,35,0,40])

r343 = subplot(2,2,3);
expPre34 = fitdist(DataPre34, 'Exponential');
histfit(DataPre34,11,'Exponential')
ExpoPre34 = mle(DataPre34, 'distribution', 'Exponential');
xlabel('Parasite Count')
ylabel('Density')
title('Exponential Distribution')

r344 = subplot(2,2,4);
poiPre34 = fitdist(DataPre34, 'Poisson');
histfit(DataPre34,11,'Poisson')
PoiPre34 = mle(DataPre34, 'distribution', 'Poisson');
xlabel('Parasite Count')
ylabel('Density')
title('Poisson Distribution')

linkaxes([r341,r343,r344],'xy')

figure;
R341 = subplot(2,2,1);
norPost34 = fitdist(DataPost34, 'Normal');
histfit(DataPost34,11,'Normal')
NormPost34 = mle(DataPost34, 'distribution', 'Normal');
xlabel('Parasite Count')
ylabel('Density')
title('Normal Distribution')
xlim([0,inf])
suptitle('Post Mouse RAG 3 Experiment 4')

R342 = subplot(2,2,2);
nbPost34 = fitdist(DataPost34, 'Negative Binomial');
histfit(DataPost34,11,'Negative Binomial')
NegBinPost34 = mle(DataPost34, 'distribution', 'Negative Binomial');
xlabel('Parasite Count')
ylabel('Density')
title('Negative Binomial Distribution')
axis([0,4305,0,40])

R343 = subplot(2,2,3);
expPost34 = fitdist(DataPost34, 'Exponential');
histfit(DataPost34,11,'Exponential')
ExpoPost34 = mle(DataPost34, 'distribution', 'Exponential');
xlabel('Parasite Count')
ylabel('Density')
title('Exponential Distribution')

R344 = subplot(2,2,4);
poiPost34 = fitdist(DataPost34, 'Poisson');
histfit(DataPost34,11,'Poisson')
PoiPost34 = mle(DataPost34, 'distribution', 'Poisson');
xlabel('Parasite Count')
ylabel('Density')
title('Poisson Distribution')

linkaxes([R341,R343,R344],'xy')

%% Experimental ID 10 Experiment 4

Pre10 = [45917.9	30642.6	3.4	39.8	203.7	4395.4	3471.5	16373.9	16270.9	23.2	7529	21822	1772	105894.9	8905	155.5	161.4	117	1.7	13632.7	57148.3	9698.5	9941.3	2118.8	60.9	59572	0.1	886.4	162	1143.7	4831.8	4164.2	46721.5	20933.6	575.1	1.3	198.9	4753.2	269.4	16621.4];
Post10 = [16948.1	201229.3	56581.7	112934.9	17244.1	71996.4	84816.5	197454.1	58578.8	0	105377.5	2719.7	6702.3	20780.8	334.7	77717.5	219061.9	90332.4	173098.7	41.3	97045.4	99730.3	124903.2	87495.9	5	19.5	34.7	59999.6	29938.8	27444.7	47092.1	3952.7	104650	150732	171.7	0	61254.3	177648.8	43357.5	29698.6];

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
distnPre10 = makedist('normal','mu',12928.4,'sigma',22023.5);
[h,p] = adtest(DataPre10,'Distribution',distnPre10);
distnPost10 = makedist('normal','mu',66478.2,'sigma',64180.1);
[h,p] = adtest(DataPost10,'Distribution',distnPost10);

distnbPre10 = makedist('negative binomial','r',0.288192,'p',0.0000222908);
[h,p] = adtest(DataPre10,'Distribution',distnbPre10);
distnbPost10 = makedist('negative binomial','r',0.345845,'p',0.00000520235);
[h,p] = adtest(DataPost10,'Distribution',distnbPost10);

distePre10 = makedist('exponential','mu',12928.4);
[h,p] = adtest(DataPre10,'Distribution',distePre10);
distePost10 = makedist('exponential','mu',66478.2);
[h,p] = adtest(DataPost10,'Distribution',distePost10);

distpPre10 = makedist('poisson','lambda',12928.4);
[h,p] = adtest(DataPre10,'Distribution',distpPre10);
distpPost10 = makedist('poisson','lambda',66478.2);
[h,p] = adtest(DataPost10,'Distribution',distpPost10);

%fitting distributions
figure;
r101 = subplot(2,2,1);
norPre10 = fitdist(DataPre10, 'Normal');
histfit(DataPre10,11,'Normal')
NormPre10 = mle(DataPre10, 'distribution', 'Normal');
xlabel('Parasite Count')
ylabel('Density')
title('Normal Distribution')
xlim([0,inf])
suptitle('Pre Mouse RAG 10 Experiment 4')

r102 = subplot(2,2,2);
nbPre10 = fitdist(DataPre10, 'Negative Binomial');
histfit(DataPre10,11,'Negative Binomial')
NegBinPre10 = mle(DataPre10, 'distribution', 'Negative Binomial');
xlabel('Parasite Count')
ylabel('Density')
title('Negative Binomial Distribution')
axis([0,101850,0,30])

r103 = subplot(2,2,3);
expPre10 = fitdist(DataPre10, 'Exponential');
histfit(DataPre10,11,'Exponential')
ExpoPre10 = mle(DataPre10, 'distribution', 'Exponential');
xlabel('Parasite Count')
ylabel('Density')
title('Exponential Distribution')

r104 = subplot(2,2,4);
poiPre10 = fitdist(DataPre10, 'Poisson');
histfit(DataPre10,11,'Poisson')
PoiPre10 = mle(DataPre10, 'distribution', 'Poisson');
xlabel('Parasite Count')
ylabel('Density')
title('Poisson Distribution')

linkaxes([r101,r103,r104],'xy')

figure;
R101 = subplot(2,2,1);
norPost10 = fitdist(DataPost10, 'Normal');
histfit(DataPost10,11,'Normal')
NormPost10 = mle(DataPost10, 'distribution', 'Normal');
xlabel('Parasite Count')
ylabel('Density')
title('Normal Distribution')
axis([0,210000,0,15])
suptitle('Post Mouse RAG 10 Experiment 4')

R102 = subplot(2,2,2);
nbPost10 = fitdist(DataPost10, 'Negative Binomial');
histfit(DataPost10,11,'Negative Binomial')
NegBinPost10 = mle(DataPost10, 'distribution', 'Negative Binomial');
xlabel('Parasite Count')
ylabel('Density')
title('Negative Binomial Distribution')
axis([0,210000,0,15])

R103 = subplot(2,2,3);
expPost10 = fitdist(DataPost10, 'Exponential');
histfit(DataPost10,11,'Exponential')
ExpoPost10 = mle(DataPost10, 'distribution', 'Exponential');
xlabel('Parasite Count')
ylabel('Density')
title('Exponential Distribution')
axis([0,210000,0,15])

R104 = subplot(2,2,4);
poiPost10 = fitdist(DataPost10, 'Poisson');
histfit(DataPost10,11,'Poisson')
PoiPost10 = mle(DataPost10, 'distribution', 'Poisson');
xlabel('Parasite Count')
ylabel('Density')
title('Poisson Distribution')
axis([0,210000,0,15])

%% Experimental ID 3 Experiment 5

Pre35 = [1394.6	23.4	16.8	459.2	127.3	3	0.7	0.1	1652.1	0.1	5.6	0.2	1.2	0	0.1	0.1	0.2	284.4	264.9	0.2	1	36.5	0.9	2.3	0.3	72.3	0.1	10.8	1.7	0.8	0.2	0	3.8	11.6	2507.7	7.3	17.8	70.9	15.7	54.2];
Post35 = [74176.4	0.8	0.3	2191.4	562.1	4431.3	0.2	5.8	138.3	9.4	11491.5	1	10.4	1.8	1422	0.5	624.6	2.2	0.1	1918.5	2.8	2840.7	4.9	13.7	22.5	127.8	0	483.7	0	0	1.3	0	403.8	2.1	0	90.4	4	0.3	0.6	5.3];

%data
DataPre35 = round(reshape(Pre35, [numel(Pre35), 1]));
DataPost35= round(reshape(Post35, [numel(Post35), 1]));

%log likelihood
norPre35 = fitdist(DataPre35, 'Normal');
xPre35 = norPre35.NLogL; %log likelihood (normal dist.)
norPost35 = fitdist(DataPost35, 'Normal');
xPost35 = norPost35.NLogL; %log likelihood (normal dist.)

nbPre35 = fitdist(DataPre35, 'Negative Binomial');
yPre35 = nbPre35.NLogL; %log likelihood (Neg. Bin. dist.)
nbPost35 = fitdist(DataPost35, 'Negative Binomial');
yPost35 = nbPost35.NLogL; %log likelihood (Neg. Bin. dist.)

expPre35 = fitdist(DataPre35, 'Exponential');
zPre35 = expPre35.NLogL; %log likelihood (exponential dist.)
expPost35 = fitdist(DataPost35, 'Exponential');
zPost35 = expPost35.NLogL; %log likelihood (exponential dist.)

poiPre35 = fitdist(DataPre35, 'Poisson');
wPre35 = poiPre35.NLogL; %log likelihood (poisson dist.)
poiPost35 = fitdist(DataPost35, 'Poisson');
wPost35 = poiPost35.NLogL; %log likelihood (poisson dist.)

%Akaike Information Criterion
aicxPre35 = 2*2-2*log(xPre35); %AIC = 2*#parameters - 2ln(ll)
aicyPre35 = 2*2-2*log(yPre35);
aiczPre35 = 2*1-2*log(zPre35);
aicwPre35 = 2*1-2*log(wPre35);
aicPre35 = aicbic([xPre35, yPre35, zPre35, wPre35],[2, 2, 1, 1]); 
aicxPost35 = 2*2-2*log(xPost35); %AIC = 2*#parameters - 2ln(ll)
aicyPost35 = 2*2-2*log(yPost35);
aiczPost35 = 2*1-2*log(zPost35);
aicwPost35 = 2*1-2*log(wPost35);
aicPost35 = aicbic([xPost35, yPost35, zPost35, wPost35],[2, 2, 1, 1]);

%Anderson-Darling test
distnPre35 = makedist('normal','mu',176.25,'sigma',509.314);
[h,p] = adtest(DataPre35,'Distribution',distnPre35);
distnPost35 = makedist('normal','mu',2524.82,'sigma',11787.1);
[h,p] = adtest(DataPost35,'Distribution',distnPost35);

distnbPre35 = makedist('negative binomial','r',0.146251,'p',0.000829106);
[h,p] = adtest(DataPre35,'Distribution',distnbPre35);
distnbPost35 = makedist('negative binomial','r',0.11377,'p',0.0000450586);
[h,p] = adtest(DataPost35,'Distribution',distnbPost35);

distePre35 = makedist('exponential','mu',176.25);
[h,p] = adtest(DataPre35,'Distribution',distePre35);
distePost35 = makedist('exponential','mu',2524.82);
[h,p] = adtest(DataPost35,'Distribution',distePost35);

distpPre35 = makedist('poisson','lambda',176.25);
[h,p] = adtest(DataPre35,'Distribution',distpPre35);
distpPost35 = makedist('poisson','lambda',2524.82);
[h,p] = adtest(DataPost35,'Distribution',distpPost35);

%fitting distributions
figure;
r351 = subplot(2,2,1);
norPre35 = fitdist(DataPre35, 'Normal');
histfit(DataPre35,11,'Normal')
NormPre35 = mle(DataPre35, 'distribution', 'Normal');
xlabel('Parasite Count')
ylabel('Density')
title('Normal Distribution')
xlim([0,inf])
suptitle('Pre Mouse RAG 3 Experiment 5')

r352 = subplot(2,2,2);
nbPre35 = fitdist(DataPre35, 'Negative Binomial');
histfit(DataPre35,11,'Negative Binomial')
NegBinPre35 = mle(DataPre35, 'distribution', 'Negative Binomial');
xlabel('Parasite Count')
ylabel('Density')
title('Negative Binomial Distribution')
axis([0,2415,0,40])

r353 = subplot(2,2,3);
expPre35 = fitdist(DataPre35, 'Exponential');
histfit(DataPre35,11,'Exponential')
ExpoPre35 = mle(DataPre35, 'distribution', 'Exponential');
xlabel('Parasite Count')
ylabel('Density')
title('Exponential Distribution')

r354 = subplot(2,2,4);
poiPre35 = fitdist(DataPre35, 'Poisson');
histfit(DataPre35,11,'Poisson')
PoiPre35 = mle(DataPre35, 'distribution', 'Poisson');
xlabel('Parasite Count')
ylabel('Density')
title('Poisson Distribution')

linkaxes([r351,r353,r354],'xy')

figure;
R351 = subplot(2,2,1);
norPost35 = fitdist(DataPost35, 'Normal');
histfit(DataPost35,11,'Normal')
NormPost35 = mle(DataPost35, 'distribution', 'Normal');
xlabel('Parasite Count')
ylabel('Density')
title('Normal Distribution')
xlim([0,inf])
suptitle('Post Mouse RAG 3 Experiment 5')

R352 = subplot(2,2,2);
nbPost35 = fitdist(DataPost35, 'Negative Binomial');
histfit(DataPost35,11,'Negative Binomial')
NegBinPost35 = mle(DataPost35, 'distribution', 'Negative Binomial');
xlabel('Parasite Count')
ylabel('Density')
title('Negative Binomial Distribution')
axis([0,71400,0,40])

R353 = subplot(2,2,3);
expPost35 = fitdist(DataPost35, 'Exponential');
histfit(DataPost35,11,'Exponential')
ExpoPost35 = mle(DataPost35, 'distribution', 'Exponential');
xlabel('Parasite Count')
ylabel('Density')
title('Exponential Distribution')

R354 = subplot(2,2,4);
poiPost35 = fitdist(DataPost35, 'Poisson');
histfit(DataPost35,11,'Poisson')
PoiPost35 = mle(DataPost35, 'distribution', 'Poisson');
xlabel('Parasite Count')
ylabel('Density')
title('Poisson Distribution')

linkaxes([R351,R353,R354],'xy')

%% Experimental ID 4 Experiment 5

Pre4 = [13.3	36.2	2.5	145	6.7	5.3	35.9	1.5	14.7	35.8	46.2	15.9	231	1.4	10	776.8	6.6	11.8	4	129.2	19.4	13.3	7.8	23.9	31.7	56.5	87	138.7	3.2	139.6	6.2	28.8	91.4	25.9	25.5	2.1	5.7	6.6	919.7	52.2];
Post4 = [0.8	1682.6	0	0	322.3	0	188.8	0.1	0.3	911.5	0	0	0	0.3	0.4	4838.1	0.6	0	0	0	0	4.9	0.1	0	0.4	0	0.3	0	0	0.1	0	0	153.4	0	0	0	0	0	0.2	2.5];

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
distnPre4 = makedist('normal','mu',80.45,'sigma',186.366);
[h,p] = adtest(DataPre4,'Distribution',distnPre4);
distnPost4 = makedist('normal','mu',202.675,'sigma',810.037);
[h,p] = adtest(DataPost4,'Distribution',distnPost4);

distnbPre4 = makedist('negative binomial','r',0.48524,'p',0.00599541);
[h,p] = adtest(DataPre4,'Distribution',distnbPre4);
distnbPost4 = makedist('negative binomial','r',0.0329398,'p',0.000162499);
[h,p] = adtest(DataPost4,'Distribution',distnbPost4);

distePre4 = makedist('exponential','mu',80.45);
[h,p] = adtest(DataPre4,'Distribution',distePre4);
distePost4 = makedist('exponential','mu',202.675);
[h,p] = adtest(DataPost4,'Distribution',distePost4);

distpPre4 = makedist('poisson','lambda',80.45);
[h,p] = adtest(DataPre4,'Distribution',distpPre4);
distpPost4 = makedist('poisson','lambda',202.675);
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
suptitle('Pre Mouse RAG 4 Experiment 5')

r42 = subplot(2,2,2);
nbPre4 = fitdist(DataPre4, 'Negative Binomial');
histfit(DataPre4,11,'Negative Binomial')
NegBinPre4 = mle(DataPre4, 'distribution', 'Negative Binomial');
xlabel('Parasite Count')
ylabel('Density')
title('Negative Binomial Distribution')
axis([0,882,0,40])

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
suptitle('Post Mouse RAG 4 Experiment 5')

R42 = subplot(2,2,2);
nbPost4 = fitdist(DataPost4, 'Negative Binomial');
histfit(DataPost4,11,'Negative Binomial')
NegBinPost4 = mle(DataPost4, 'distribution', 'Negative Binomial');
xlabel('Parasite Count')
ylabel('Density')
title('Negative Binomial Distribution')
axis([0,4620,0,40])

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

%% Experimental ID 7 Experiment 5

Pre7 = [0.2	19.9	3.5	0.2	0.8	17.7	105.2	86	117	0.6	86.4	3	20.9	6.2	0	15.1	34.1	3.3	0.3	0	0.8	2	1.4	1.7	8.1	38.8	156.7	0	0.5	0	8.7	1.2	10.2	10.1	0	0.6	0	8.8	0	0];
Post7 = [0.1	6.2	48.7	0.1	3599.4	13.5	0	0	1.1	0	1.2	0	0	0	0.2	0	0	0	0	0	0	0	0	0.3	0	0	1638.6	80.9	0	0.5	0	0	0.2	0	1.3	0	0	0	2	0];

%data
DataPre7 = round(reshape(Pre7, [numel(Pre7), 1]));
DataPost7 = round(reshape(Post7, [numel(Post7), 1]));

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
distnPre7 = makedist('normal','mu',19.275,'sigma',37.1663);
[h,p] = adtest(DataPre7,'Distribution',distnPre7);
distnPost7 = makedist('normal','mu',134.85,'sigma',618.537);
[h,p] = adtest(DataPost7,'Distribution',distnPost7);

distnbPre7 = makedist('negative binomial','r',0.269991,'p',0.0138138);
[h,p] = adtest(DataPre7,'Distribution',distnbPre7);
distnbPost7 = makedist('negative binomial','r',0.0372163,'p',0.000275907);
[h,p] = adtest(DataPost7,'Distribution',distnbPost7);

distePre7 = makedist('exponential','mu',19.275);
[h,p] = adtest(DataPre7,'Distribution',distePre7);
distePost7 = makedist('exponential','mu',134.85);
[h,p] = adtest(DataPost7,'Distribution',distePost7);

distpPre7 = makedist('poisson','lambda',19.275);
[h,p] = adtest(DataPre7,'Distribution',distpPre7);
distpPost7 = makedist('poisson','lambda',134.85);
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
suptitle('Pre Mouse RAG 7 Experiment 5')

r72 = subplot(2,2,2);
nbPre7 = fitdist(DataPre7, 'Negative Binomial');
histfit(DataPre7,11,'Negative Binomial')
NegBinPre7 = mle(DataPre7, 'distribution', 'Negative Binomial');
xlabel('Parasite Count')
ylabel('Density')
title('Negative Binomial Distribution')
axis([0,160,0,30])

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
suptitle('Post Mouse RAG 7 Experiment 5')

R72 = subplot(2,2,2);
nbPost7 = fitdist(DataPost7, 'Negative Binomial');
histfit(DataPost7,11,'Negative Binomial')
NegBinPost7 = mle(DataPost7, 'distribution', 'Negative Binomial');
xlabel('Parasite Count')
ylabel('Density')
title('Negative Binomial Distribution')
axis([0,3465,0,40])

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

%% Experimental ID 9 Experiment 5

Pre9 = [66.1	595.8	6.9	129.4	11.4	1.6	50.9	107.9	19.4	339.3	10	19.9	12.6	11.4	19.9	133.3	52.6	20.2	28.5	5.8	550.1	1675.3	22.3	112.7	382	1302.6	5.3	6.9	4414.4	10	2.4	28.8	503.5	387.7	3	15.9	2003.7	25.3	59.3	319];
Post9 = [0	431.8	29.3	5762.4	0	0.2	513.4	0.8	0.1	0	217	39.5	83.3	284.2	79.4	0	0	339.7	0.2	1.1	1549.5	0.2	0.3	117.9	3.2	288.8	0.6	0.7	207.4	0.1	2.2	210.9	0.3	0	141.3	2.4	63.2	79.4	29.4	0.2];

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
distnPre9 = makedist('normal','mu',336.825,'sigma',800.114);
[h,p] = adtest(DataPre9,'Distribution',distnPre9);
distnPost9 = makedist('normal','mu',261.925,'sigma',930.721);
[h,p] = adtest(DataPost9,'Distribution',distnPost9);

distnbPre9 = makedist('negative binomial','r',0.355801,'p',0.00105522);
[h,p] = adtest(DataPre9,'Distribution',distnbPre9);
distnbPost9 = makedist('negative binomial','r',0.134039,'p',0.000511485);
[h,p] = adtest(DataPost9,'Distribution',distnbPost9);

distePre9 = makedist('exponential','mu',336.825);
[h,p] = adtest(DataPre9,'Distribution',distePre9);
distePost9 = makedist('exponential','mu',261.925);
[h,p] = adtest(DataPost9,'Distribution',distePost9);

distpPre9 = makedist('poisson','lambda',336.825);
[h,p] = adtest(DataPre9,'Distribution',distpPre9);
distpPost9 = makedist('poisson','lambda',261.925);
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
suptitle('Pre Mouse RAG 9 Experiment 5')

r92 = subplot(2,2,2);
nbPre9 = fitdist(DataPre9, 'Negative Binomial');
histfit(DataPre9,11,'Negative Binomial')
NegBinPre9 = mle(DataPre9, 'distribution', 'Negative Binomial');
xlabel('Parasite Count')
ylabel('Density')
title('Negative Binomial Distribution')
axis([0,inf,0,40])

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
suptitle('Post Mouse RAG 9 Experiment 5')

R92 = subplot(2,2,2);
nbPost9 = fitdist(DataPost9, 'Negative Binomial');
histfit(DataPost9,11,'Negative Binomial')
NegBinPost9 = mle(DataPost9, 'distribution', 'Negative Binomial');
xlabel('Parasite Count')
ylabel('Density')
title('Negative Binomial Distribution')
axis([0,5600,0,40])

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

%% Experimental ID 11 Experiment 5

Pre11 = [93.2	14.5	1580.5	710	55.2	1267.4	19.3	29.5	148.7	75.8	66.5	72.1	602.6	13	7.5	259.4	5103.3	2.9	273.8	1265.1	1084.2	748.5	3.9	6.1	8.9	52.6	22.2	849.8	6303.9	105.8	405.3	2809.8	9.5	18.5	4.9	15.4	5.4	1034.1	7.8	652.5];
Post11 = [0.5	0.9	125.4	789.6	682.5	1378.6	44.2	508.9	7.2	16.9	573.6	31.9	55.5	8.1	738.5	0.6	27265.5	2.8	2429.9	104.3	1565.9	4.8	509.3	33.5	6291.7	5.1	586.1	11.2	2825.6	13	0.1	13644.3	5	0.2	0	0.9	1.3	886.5	4.6	2809.3];

%data
DataPre11 = round(reshape(Pre11, [numel(Pre11), 1]));
DataPost11 = round(reshape(Post11, [numel(Post11), 1]));

%log likelihood
norPre11 = fitdist(DataPre11, 'Normal');
xPre11 = norPre11.NLogL; %log likelihood (normal dist.)
norPost11 = fitdist(DataPost11, 'Normal');
xPost11 = norPost11.NLogL; %log likelihood (normal dist.)

nbPre11 = fitdist(DataPre11, 'Negative Binomial');
yPre11 = nbPre11.NLogL; %log likelihood (Neg. Bin. dist.)
nbPost11 = fitdist(DataPost11, 'Negative Binomial');
yPost11 = nbPost11.NLogL; %log likelihood (Neg. Bin. dist.)

expPre11 = fitdist(DataPre11, 'Exponential');
zPre11 = expPre11.NLogL; %log likelihood (exponential dist.)
expPost11 = fitdist(DataPost11, 'Exponential');
zPost11 = expPost11.NLogL; %log likelihood (exponential dist.)

poiPre11 = fitdist(DataPre11, 'Poisson');
wPre11 = poiPre11.NLogL; %log likelihood (poisson dist.)
poiPost11 = fitdist(DataPost11, 'Poisson');
wPost11 = poiPost11.NLogL; %log likelihood (poisson dist.)

%Akaike Information Criterion
aicxPre11 = 2*2-2*log(xPre11); %AIC = 2*#parameters - 2ln(ll)
aicyPre11 = 2*2-2*log(yPre11);
aiczPre11 = 2*1-2*log(zPre11);
aicwPre11 = 2*1-2*log(wPre11);
aicPre11 = aicbic([xPre11, yPre11, zPre11, wPre11],[2, 2, 1, 1]); 
aicxPost11 = 2*2-2*log(xPost11); %AIC = 2*#parameters - 2ln(ll)
aicyPost11 = 2*2-2*log(yPost11);
aiczPost11 = 2*1-2*log(zPost11);
aicwPost11 = 2*1-2*log(wPost11);
aicPost11 = aicbic([xPost11, yPost11, zPost11, wPost11],[2, 2, 1, 1]);

%Anderson-Darling test
distnPre11 = makedist('normal','mu',645.325,'sigma',1319.66);
[h,p] = adtest(DataPre11,'Distribution',distnPre11);
distnPost11 = makedist('normal','mu',1640.21,'sigma',4850.52);
[h,p] = adtest(DataPost11,'Distribution',distnPost11);

distnbPre11 = makedist('negative binomial','r',0.350129,'p',0.000542268);
[h,p] = adtest(DataPre11,'Distribution',distnbPre11);
distnbPost11 = makedist('negative binomial','r',0.205923,'p',0.000125531);
[h,p] = adtest(DataPost11,'Distribution',distnbPost11);

distePre11 = makedist('exponential','mu',645.325);
[h,p] = adtest(DataPre11,'Distribution',distePre11);
distePost11 = makedist('exponential','mu',1640.21);
[h,p] = adtest(DataPost11,'Distribution',distePost11);

distpPre11 = makedist('poisson','lambda',645.325);
[h,p] = adtest(DataPre11,'Distribution',distpPre11);
distpPost11 = makedist('poisson','lambda',1640.21);
[h,p] = adtest(DataPost11,'Distribution',distpPost11);

%fitting distributions
figure;
r111 = subplot(2,2,1);
norPre11 = fitdist(DataPre11, 'Normal');
histfit(DataPre11,11,'Normal')
NormPre11 = mle(DataPre11, 'distribution', 'Normal');
xlabel('Parasite Count')
ylabel('Density')
title('Normal Distribution')
xlim([0,inf])
suptitle('Pre Mouse RAG 11 Experiment 5')

r112 = subplot(2,2,2);
nbPre11 = fitdist(DataPre11, 'Negative Binomial');
histfit(DataPre11,11,'Negative Binomial')
NegBinPre11 = mle(DataPre11, 'distribution', 'Negative Binomial');
xlabel('Parasite Count')
ylabel('Density')
title('Negative Binomial Distribution')
axis([0,6170,0,30])

r113 = subplot(2,2,3);
expPre11 = fitdist(DataPre11, 'Exponential');
histfit(DataPre11,11,'Exponential')
ExpoPre11 = mle(DataPre11, 'distribution', 'Exponential');
xlabel('Parasite Count')
ylabel('Density')
title('Exponential Distribution')

r114 = subplot(2,2,4);
poiPre11 = fitdist(DataPre11, 'Poisson');
histfit(DataPre11,11,'Poisson')
PoiPre11 = mle(DataPre11, 'distribution', 'Poisson');
xlabel('Parasite Count')
ylabel('Density')
title('Poisson Distribution')

linkaxes([r111,r113,r114],'xy')

figure;
R111 = subplot(2,2,1);
norPost11 = fitdist(DataPost11, 'Normal');
histfit(DataPost11,11,'Normal')
NormPost11 = mle(DataPost11, 'distribution', 'Normal');
xlabel('Parasite Count')
ylabel('Density')
title('Normal Distribution')
xlim([0,inf])
suptitle('Post Mouse RAG 11 Experiment 5')

R112 = subplot(2,2,2);
nbPost11 = fitdist(DataPost11, 'Negative Binomial');
histfit(DataPost11,11,'Negative Binomial')
NegBinPost11 = mle(DataPost11, 'distribution', 'Negative Binomial');
xlabel('Parasite Count')
ylabel('Density')
title('Negative Binomial Distribution')
axis([0,28750,0,40])

R113 = subplot(2,2,3);
expPost11 = fitdist(DataPost11, 'Exponential');
histfit(DataPost11,11,'Exponential')
ExpoPost11 = mle(DataPost11, 'distribution', 'Exponential');
xlabel('Parasite Count')
ylabel('Density')
title('Exponential Distribution')

R114 = subplot(2,2,4);
poiPost11 = fitdist(DataPost11, 'Poisson');
histfit(DataPost11,11,'Poisson')
PoiPost11 = mle(DataPost11, 'distribution', 'Poisson');
xlabel('Parasite Count')
ylabel('Density')
title('Poisson Distribution')

linkaxes([R111,R113,R114],'xy')
