%CD4 mice

%% Experimental ID A Experiment 2

MA = [6242.309	40213.92	31013.05	208242.7	60042.1	98883.61	60476.59	249944.2	10892.03	6057.031	112096.5	35391.51	86337.55	45613.8	46020.12	64171.25	56225.69	20209.33	55936.79	42591.23	88499.5	5941.761	5517.263	139639];

%data
DataA = round(reshape(MA, [numel(MA), 1]));

%log likelihood
norA = fitdist(DataA, 'Normal');
xA = norA.NLogL; %log likelihood (normal dist.)

nbA = fitdist(DataA, 'Negative Binomial');
yA = nbA.NLogL; %log likelihood (Neg. Bin. dist.)

expA = fitdist(DataA, 'Exponential');
zA = expA.NLogL; %log likelihood (exponential dist.)

poiA = fitdist(DataA, 'Poisson');
wA = poiA.NLogL; %log likelihood (poisson dist.)

%Akaike Information Criterion
aicxA = 2*2-2*log(xA); %AIC = 2*#parameters - 2ln(ll)
aicyA = 2*2-2*log(yA);
aiczA = 2*1-2*log(zA);
aicwA = 2*1-2*log(wA);
aicA = aicbic([xA, yA, zA, wA],[2, 2, 1, 1]); 

%Anderson-Darling test
distnA = makedist('normal','mu',65675,'sigma',61726.1);
[h,p] = adtest(DataA,'Distribution',distnA);

distnbA = makedist('negative binomial','r',1.19602,'p',0.0000182109);
[h,p] = adtest(DataA,'Distribution',distnbA);

disteA = makedist('exponential','mu',65675);
[h,p] = adtest(DataA,'Distribution',disteA);

distpA = makedist('poisson','lambda',65675);
[h,p] = adtest(DataA,'Distribution',distpA);

%fitting distributions
figure;
a1 = subplot(2,2,1);
histfit(DataA,11,'Normal')
NormA = mle(DataA, 'distribution', 'Normal');
xlabel('Parasite Count')
ylabel('Density')
title('Normal Distribution')
xlim([0,inf])
suptitle('Mouse A')

a2 = subplot(2,2,2);
histfit(DataA,11,'Negative Binomial')
NegBinA = mle(DataA, 'distribution', 'Negative Binomial');
xlabel('Parasite Count')
ylabel('Density')
title('Negative Binomial Distribution')
axis([0,253000,0,6])

a3 = subplot(2,2,3);
histfit(DataA,11,'Exponential')
ExpoA = mle(DataA, 'distribution', 'Exponential');
xlabel('Parasite Count')
ylabel('Density')
title('Exponential Distribution')
axis([0,253000,0,6])

a4 = subplot(2,2,4);
histfit(DataA,11,'Poisson')
PoiA = mle(DataA, 'distribution', 'Poisson');
xlabel('Parasite Count')
ylabel('Density')
title('Poisson Distribution')

linkaxes([a1,a4],'xy')

%% Experimental ID B Experiment 2

MB = [113652.9	23657.08	7039.903	67137.05	64529.57	131615.3	518365.6	541672.1	1800.071	14785.82	1394753	706560.7	130131.8	23352.05	208169.6	608231.4	106441.2	32235.32	44494.23	94896.09	6424.299	5645.563	27523.34	32097.01];

%data
DataB = round(reshape(MB, [numel(MB), 1]));

%log likelihood
norB = fitdist(DataB, 'Normal');
xB = norB.NLogL; %log likelihood (normal dist.)

nbB = fitdist(DataB, 'Negative Binomial');
yB = nbB.NLogL; %log likelihood (Neg. Bin. dist.)

expB = fitdist(DataB, 'Exponential');
zB = expB.NLogL; %log likelihood (exponential dist.)

poiB = fitdist(DataB, 'Poisson');
wB = poiB.NLogL; %log likelihood (poisson dist.)

%Akaike Information Criterion
aicxB = 2*2-2*log(xB); %AIC = 2*#parameters - 2ln(ll)
aicyB = 2*2-2*log(yB);
aiczB = 2*1-2*log(zB);
aicwB = 2*1-2*log(wB);
aicB = aicbic([xB, yB, zB, wB],[2, 2, 1, 1]); 

%Anderson-Darling test
distnB = makedist('normal','mu',204384,'sigma',329680);
[h,p] = adtest(DataB,'Distribution',distnB);

distnbB = makedist('negative binomial','r',0.533042,'p',0.00000260804);
[h,p] = adtest(DataB,'Distribution',distnbB);

disteB = makedist('exponential','mu',204384);
[h,p] = adtest(DataB,'Distribution',disteB);

distpB = makedist('poisson','lambda',204384);
[h,p] = adtest(DataB,'Distribution',distpB);

% fitting distributions
figure;
b1 = subplot(2,2,1);
histfit(DataB,11,'Normal')
NormB = mle(DataB, 'distribution', 'Normal');
xlabel('Parasite Count')
ylabel('Density')
title('Normal Distribution')
xlim([0,inf])
suptitle('Mouse B')

b2 = subplot(2,2,2);
histfit(DataB,11,'Negative Binomial')
NegBinB = mle(DataB, 'distribution', 'Negative Binomial');
xlabel('Parasite Count')
ylabel('Density')
title('Negative Binomial Distribution')
axis([0,1400000,0,16])

b3 = subplot(2,2,3);
histfit(DataB,11,'Exponential')
ExpoB = mle(DataB, 'distribution', 'Exponential');
xlabel('Parasite Count')
ylabel('Density')
title('Exponential Distribution')

b4 = subplot(2,2,4);
histfit(DataB,11,'Poisson')
PoiB = mle(DataB, 'distribution', 'Poisson');
xlabel('Parasite Count')
ylabel('Density')
title('Poisson Distribution')
%axis([0,1430000,0,16])

linkaxes([b1,b3,b4],'xy')

%% Experimental ID 1 Experiment 4

M14 = [18169.2	12718.1	16364.7	34639.3	102337.6	79552.5	115456.1	262826.1	81365.2	110573.1	209569.7	67009.6	66966.9	29948.9	232116.4	766783.9	63812.4	36401	77071.4	56991.7	14903.7	60864.6	8198.5	1496.9];

%data
Data14 = round(reshape(M14, [numel(M14), 1]));

%log likelihood
nor14 = fitdist(Data14, 'Normal');
x14 = nor14.NLogL; %log likelihood (normal dist.)

nb14 = fitdist(Data14, 'Negative Binomial');
y14 = nb14.NLogL; %log likelihood (Neg. Bin. dist.)

exp14 = fitdist(Data14, 'Exponential');
z14 = exp14.NLogL; %log likelihood (exponential dist.)

poi14 = fitdist(Data14, 'Poisson');
w14 = poi14.NLogL; %log likelihood (poisson dist.)

%Akaike Information Criterion
aicx14 = 2*2-2*log(x14); %AIC = 2*#parameters - 2ln(ll)
aicy14 = 2*2-2*log(y14);
aicz14 = 2*1-2*log(z14);
aicw14 = 2*1-2*log(w14);
aic14 = aicbic([x14, y14, z14, w14],[2, 2, 1, 1]); 

%Anderson-Darling test
distn14 = makedist('normal','mu',105256,'sigma',157195);
[h,p] = adtest(Data14,'Distribution',distn14);

distnb14 = makedist('negative binomial','r',0.840364,'p',0.00000798395);
[h,p] = adtest(Data14,'Distribution',distnb14);

diste14 = makedist('exponential','mu',105256);
[h,p] = adtest(Data14,'Distribution',diste14);

distp14 = makedist('poisson','lambda',105256);
[h,p] = adtest(Data14,'Distribution',distp14);

%fitting distributions
figure;
r141 = subplot(2,2,1);
histfit(Data14,11,'Normal')
Norm14 = mle(Data14, 'distribution', 'Normal');
xlabel('Parasite Count')
ylabel('Density')
title('Normal Distribution')
xlim([0,inf])
suptitle('Mouse RAG 1 Experiment 4')

r142 = subplot(2,2,2);
histfit(Data14,11,'Negative Binomial')
NegBin14 = mle(Data14, 'distribution', 'Negative Binomial');
xlabel('Parasite Count')
ylabel('Density')
title('Negative Binomial Distribution')
axis([0,770000,0,15])

r143 = subplot(2,2,3);
histfit(Data14,11,'Exponential')
Expo14 = mle(Data14, 'distribution', 'Exponential');
xlabel('Parasite Count')
ylabel('Density')
title('Exponential Distribution')

r144 = subplot(2,2,4);
histfit(Data14,11,'Poisson')
Poi14 = mle(Data14, 'distribution', 'Poisson');
xlabel('Parasite Count')
ylabel('Density')
title('Poisson Distribution')

linkaxes([r141,r143,r144],'xy')

%% Experimental ID 2 Experiment 4

M2 = [44545.2	530171.5	668891.6	2469231	19219.4	571615.5	1491697	725753	42283.1	263689.3	1826778	318708.8	164573.9	91659.4	462878.5	60530.8	7163.7	47056.6	50569.9	38927.7	3931.7	6870.1	23407.3	26510.7];

%data
Data2 = round(reshape(M2, [numel(M2), 1]));

%log likelihood
nor2 = fitdist(Data2, 'Normal');
x2 = nor2.NLogL; %log likelihood (normal dist.)

nb2 = fitdist(Data2, 'Negative Binomial');
y2 = nb2.NLogL; %log likelihood (Neg. Bin. dist.)

exp2 = fitdist(Data2, 'Exponential');
z2 = exp2.NLogL; %log likelihood (exponential dist.)

poi2 = fitdist(Data2, 'Poisson');
w2 = poi2.NLogL; %log likelihood (poisson dist.)

%Akaike Information Criterion
aicx2 = 2*2-2*log(x2); %AIC = 2*#parameters - 2ln(ll)
aicy2 = 2*2-2*log(y2);
aicz2 = 2*1-2*log(z2);
aicw2 = 2*1-2*log(w2);
aic2 = aicbic([x2, y2, z2, w2],[2, 2, 1, 1]); 

%Anderson-Darling test
distn2 = makedist('normal','mu',414861,'sigma',644273);
[h,p] = adtest(Data2,'Distribution',distn2);

distnb2 = makedist('negative binomial','r',0.486133,'p',0.0000011718);
[h,p] = adtest(Data2,'Distribution',distnb2);

diste2 = makedist('exponential','mu',414861);
[h,p] = adtest(Data2,'Distribution',diste2);

distp2 = makedist('poisson','lambda',414861);
[h,p] = adtest(Data2,'Distribution',distp2);

%fitting distributions
figure;
r21 = subplot(3,1,1);
histfit(Data2,11,'Normal')
Norm2 = mle(Data2, 'distribution', 'Normal');
xlabel('Parasite Count')
ylabel('Density')
title('Normal Distribution')
xlim([0,inf])
suptitle('Mouse RAG 2 Experiment 4')

% r22 = subplot(2,2,2);
% histfit(Data2,11,'Negative Binomial')
% NegBin2 = mle(Data2, 'distribution', 'Negative Binomial');
% xlabel('Parasite Count')
% ylabel('Density')
% title('Negative Binomial Distribution')
% axis([0,2530000,0,15])

r23 = subplot(3,1,2);
histfit(Data2,11,'Exponential')
Expo2 = mle(Data2, 'distribution', 'Exponential');
xlabel('Parasite Count')
ylabel('Density')
title('Exponential Distribution')
axis([0,2530000,0,15])

r24 = subplot(3,1,3);
histfit(Data2,11,'Poisson')
Poi2 = mle(Data2, 'distribution', 'Poisson');
xlabel('Parasite Count')
ylabel('Density')
title('Poisson Distribution')

linkaxes([r21,r24],'xy')

%% Experimental ID 6 Experiment 4

M64 = [598196	1926035	488947.2	173848.7	553694.9	4061678	14662020	3776093	90795.3	3336097	3694323	9669794	293533.2	3504126	10284320	5941550	450755.2	4809105	3664616	1397573	12728.2	49876.7	67301.9	70181.8];

%data
Data64 = round(reshape(M64, [numel(M64), 1]));

%log likelihood
nor64 = fitdist(Data64, 'Normal');
x64 = nor64.NLogL; %log likelihood (normal dist.)

nb64 = fitdist(Data64, 'Negative Binomial');
y64 = nb64.NLogL; %log likelihood (Neg. Bin. dist.)

exp64 = fitdist(Data64, 'Exponential');
z64 = exp64.NLogL; %log likelihood (exponential dist.)

poi64 = fitdist(Data64, 'Poisson');
w64 = poi64.NLogL; %log likelihood (poisson dist.)

%Akaike Information Criterion
aicx64 = 2*2-2*log(x64); %AIC = 2*#parameters - 2ln(ll)
aicy64 = 2*2-2*log(y64);
aicz64 = 2*1-2*log(z64);
aicw64 = 2*1-2*log(w64);
aic64 = aicbic([x64, y64, z64, w64],[2, 2, 1, 1]); 

%Anderson-Darling test
distn64 = makedist('normal','mu',3065720,'sigma',3817910);
[h,p] = adtest(Data64,'Distribution',distn64);

distnb64 = makedist('negative binomial','r',0.535471,'p',0.000000174664);
[h,p] = adtest(Data64,'Distribution',distnb64);

diste64 = makedist('exponential','mu',3065720);
[h,p] = adtest(Data64,'Distribution',diste64);

distp64 = makedist('poisson','lambda',3065720);
[h,p] = adtest(Data64,'Distribution',distp64);

%fitting distributions
figure;
r641 = subplot(3,1,1);
histfit(Data64,11,'Normal')
Norm64 = mle(Data64, 'distribution', 'Normal');
xlabel('Parasite Count')
ylabel('Density')
title('Normal Distribution')
axis([0,15400000,0,12])
suptitle('Mouse RAG 6 Experiment 4')

% r642 = subplot(2,2,2);
% histfit(Data64,11,'Negative Binomial')
% NegBin64 = mle(Data64, 'distribution', 'Negative Binomial');
% xlabel('Parasite Count')
% ylabel('Density')
% title('Negative Binomial Distribution')
% axis([0,15400000,0,12])

r643 = subplot(3,1,2);
histfit(Data64,11,'Exponential')
Expo64 = mle(Data64, 'distribution', 'Exponential');
xlabel('Parasite Count')
ylabel('Density')
title('Exponential Distribution')
axis([0,15400000,0,12])

r644 = subplot(3,1,3);
histfit(Data64,11,'Poisson')
Poi64 = mle(Data64, 'distribution', 'Poisson');
xlabel('Parasite Count')
ylabel('Density')
title('Poisson Distribution')
axis([0,15400000,0,12])

%% Experimental ID 8 Experiment 4

M8 = [3137887	6669272	1751551	2126678	6912118	5712355	15197510	17344610	7932178	4256139	80643520	11343150	8134356	6343389	7060130	6448552	3770326	2546905	6616066	6271906	185902.1	74053.5	248866.6	341170.5];

%data
Data8 = round(reshape(M8, [numel(M8), 1]));

%log likelihood
nor8 = fitdist(Data8, 'Normal');
x8 = nor8.NLogL; %log likelihood (normal dist.)

nb8 = fitdist(Data8, 'Negative Binomial');
y8 = nb8.NLogL; %log likelihood (Neg. Bin. dist.)

exp8 = fitdist(Data8, 'Exponential');
z8 = exp8.NLogL; %log likelihood (exponential dist.)

poi8 = fitdist(Data8, 'Poisson');
w8 = poi8.NLogL; %log likelihood (poisson dist.)

%Akaike Information Criterion
aicx8 = 2*2-2*log(x8); %AIC = 2*#parameters - 2ln(ll)
aicy8 = 2*2-2*log(y8);
aicz8 = 2*1-2*log(z8);
aicw8 = 2*1-2*log(w8);
aic8 = aicbic([x8, y8, z8, w8],[2, 2, 1, 1]); 

% %Anderson-Darling test
% distn8 = makedist('normal','mu',8794520,'sigma',15924900);
% [h,p] = adtest(Data8,'Distribution',distn8);
% 
% distnb8 = makedist('negative binomial','r',0.685018,'p',0.0000000778914);
% [h,p] = adtest(Data8,'Distribution',distnb8);
% 
% diste8 = makedist('exponential','mu',8794520);
% [h,p] = adtest(Data8,'Distribution',diste8);
% 
% distp8 = makedist('poisson','lambda',8794520);
% [h,p] = adtest(Data8,'Distribution',distp8);

%fitting distributions
figure;
r81 = subplot(3,1,1);
histfit(Data8,11,'Normal')
Norm8 = mle(Data8, 'distribution', 'Normal');
xlabel('Parasite Count')
ylabel('Density')
title('Normal Distribution')
xlim([0,inf])
suptitle('Mouse RAG 8 Experiment 4')

% r82 = subplot(2,2,2);
% histfit(Data8,11,'Negative Binomial')
% NegBin8 = mle(Data8, 'distribution', 'Negative Binomial');
% xlabel('Parasite Count')
% ylabel('Density')
% title('Negative Binomial Distribution')
% axis([0,49450,0,30])

r83 = subplot(3,1,2);
histfit(Data8,11,'Exponential')
Expo8 = mle(Data8, 'distribution', 'Exponential');
xlabel('Parasite Count')
ylabel('Density')
title('Exponential Distribution')

r84 = subplot(3,1,3);
histfit(Data8,11,'Poisson')
Poi8 = mle(Data8, 'distribution', 'Poisson');
xlabel('Parasite Count')
ylabel('Density')
title('Poisson Distribution')

linkaxes([r81,r83,r84],'xy')

%% Experimental ID 1 Experiment 5

M15 = [43606.68	39610.4	5482.782	31639.34	32049.93	33995.91	37138.89	172406.8	28486.92	33673.63	50312.51	12606.38	123821	152871.1	83785.74	27238.4	476785.3	79134.71	36514.55	45777.34	13625.87	20583.27	4612.97	140612.8];

%data
Data15 = round(reshape(M15, [numel(M15), 1]));

%log likelihood
nor15 = fitdist(Data15, 'Normal');
x15 = nor15.NLogL; %log likelihood (normal dist.)

nb15 = fitdist(Data15, 'Negative Binomial');
y15 = nb15.NLogL; %log likelihood (Neg. Bin. dist.)

exp15 = fitdist(Data15, 'Exponential');
z15 = exp15.NLogL; %log likelihood (exponential dist.)

poi15 = fitdist(Data15, 'Poisson');
w15 = poi15.NLogL; %log likelihood (poisson dist.)

%Akaike Information Criterion
aicx15 = 2*2-2*log(x15); %AIC = 2*#parameters - 2ln(ll)
aicy15 = 2*2-2*log(y15);
aicz15 = 2*1-2*log(z15);
aicw15 = 2*1-2*log(w15);
aic15 = aicbic([x15, y15, z15, w15],[2, 2, 1, 1]); 

%Anderson-Darling test
distn15 = makedist('normal','mu',71932.3,'sigma',98272.6);
[h,p] = adtest(Data15,'Distribution',distn15);

distnb15 = makedist('negative binomial','r',1.03875,'p',0.0000144405);
[h,p] = adtest(Data15,'Distribution',distnb15);

diste15 = makedist('exponential','mu',71932.3);
[h,p] = adtest(Data15,'Distribution',diste15);

distp15 = makedist('poisson','lambda',71932.3);
[h,p] = adtest(Data15,'Distribution',distp15);

%fitting distributions
figure;
r151 = subplot(2,2,1);
histfit(Data15,11,'Normal')
Norm15 = mle(Data15, 'distribution', 'Normal');
xlabel('Parasite Count')
ylabel('Density')
title('Normal Distribution')
xlim([0,inf])
suptitle('Mouse RAG 1 Experiment 5')

r152 = subplot(2,2,2);
histfit(Data15,11,'Negative Binomial')
NegBin15 = mle(Data15, 'distribution', 'Negative Binomial');
xlabel('Parasite Count')
ylabel('Density')
title('Negative Binomial Distribution')
axis([0,475000,0,20])

r153 = subplot(2,2,3);
histfit(Data15,11,'Exponential')
Expo15 = mle(Data15, 'distribution', 'Exponential');
xlabel('Parasite Count')
ylabel('Density')
title('Exponential Distribution')

r154 = subplot(2,2,4);
histfit(Data15,11,'Poisson')
Poi15 = mle(Data15, 'distribution', 'Poisson');
xlabel('Parasite Count')
ylabel('Density')
title('Poisson Distribution')

linkaxes([r151,r153,r154],'xy')

%% Experimental ID 6 Experiment 5

M65 = [2260.624	9063.668	3060.94	9444.282	4412.582	3038.998	2265.217	970.2748	7934.078	4313.613	7064.852	1514.515	15138.66	8959.973	9887.268	790.0161	7823.779	2604.275	6803.212	2141.352	4012.497	8171.049	906.79	4226.915];

%data
Data65 = round(reshape(M65, [numel(M65), 1]));

%log likelihood
nor65 = fitdist(Data65, 'Normal');
x65 = nor65.NLogL; %log likelihood (normal dist.)

nb65 = fitdist(Data65, 'Negative Binomial');
y65 = nb65.NLogL; %log likelihood (Neg. Bin. dist.)

exp65 = fitdist(Data65, 'Exponential');
z65 = exp65.NLogL; %log likelihood (exponential dist.)

poi65 = fitdist(Data65, 'Poisson');
w65 = poi65.NLogL; %log likelihood (poisson dist.)

%Akaike Information Criterion
aicx65 = 2*2-2*log(x65); %AIC = 2*#parameters - 2ln(ll)
aicy65 = 2*2-2*log(y65);
aicz65 = 2*1-2*log(z65);
aicw65 = 2*1-2*log(w65);
aic65 = aicbic([x65, y65, z65, w65],[2, 2, 1, 1]); 

%Anderson-Darling test
distn65 = makedist('normal','mu',5283.75,'sigma',3686.66);
[h,p] = adtest(Data65,'Distribution',distn65);

distnb65 = makedist('negative binomial','r',1.9186,'p',0.000362981);
[h,p] = adtest(Data65,'Distribution',distnb65);

diste65 = makedist('exponential','mu',5283.75);
[h,p] = adtest(Data65,'Distribution',diste65);

distp65 = makedist('poisson','lambda',5283.75);
[h,p] = adtest(Data65,'Distribution',distp65);

%fitting distributions
figure;
r651 = subplot(2,2,1);
histfit(Data65,11,'Normal')
Norm65 = mle(Data65, 'distribution', 'Normal');
xlabel('Parasite Count')
ylabel('Density')
title('Normal Distribution')
axis([0,15400,0,6])
suptitle('Mouse RAG 6 Experiment 5')

r652 = subplot(2,2,2);
histfit(Data65,11,'Negative Binomial')
NegBin65 = mle(Data65, 'distribution', 'Negative Binomial');
xlabel('Parasite Count')
ylabel('Density')
title('Negative Binomial Distribution')
axis([0,15400,0,6])

r653 = subplot(2,2,3);
histfit(Data65,11,'Exponential')
Expo65 = mle(Data65, 'distribution', 'Exponential');
xlabel('Parasite Count')
ylabel('Density')
title('Exponential Distribution')
axis([0,15400,0,6])

r654 = subplot(2,2,4);
histfit(Data65,11,'Poisson')
Poi65 = mle(Data65, 'distribution', 'Poisson');
xlabel('Parasite Count')
ylabel('Density')
title('Poisson Distribution')
axis([0,15400,0,6])

%% Experimental ID 10 Experiment 5

M10 = [34873.48	178186	24155.73	12365.47	15526.27	22582.66	38169.07	6815.063	4903.425	36461.53	36999.52	27371.52	23405.3	130282.7	31532.11	41119.96	22138.27	42285.54	67661.07	53227.32	6173.713	7604.093	31383.42	45624.58];

%data
Data10 = round(reshape(M10, [numel(M10), 1]));

%log likelihood
nor10 = fitdist(Data10, 'Normal');
x10 = nor10.NLogL; %log likelihood (normal dist.)

nb10 = fitdist(Data10, 'Negative Binomial');
y10 = nb10.NLogL; %log likelihood (Neg. Bin. dist.)

exp10 = fitdist(Data10, 'Exponential');
z10 = exp10.NLogL; %log likelihood (exponential dist.)

poi10 = fitdist(Data10, 'Poisson');
w10 = poi10.NLogL; %log likelihood (poisson dist.)

%Akaike Information Criterion
aicx10 = 2*2-2*log(x10); %AIC = 2*#parameters - 2ln(ll)
aicy10 = 2*2-2*log(y10);
aicz10 = 2*1-2*log(z10);
aicw10 = 2*1-2*log(w10);
aic10 = aicbic([x10, y10, z10, w10],[2, 2, 1, 1]); 

%Anderson-Darling test
distn10 = makedist('normal','mu',39202,'sigma',39344.7);
[h,p] = adtest(Data10,'Distribution',distn10);

distnb10 = makedist('negative binomial','r',1.5164,'p',0.0000386801);
[h,p] = adtest(Data10,'Distribution',distnb10);

diste10 = makedist('exponential','mu',39202);
[h,p] = adtest(Data10,'Distribution',diste10);

distp10 = makedist('poisson','lambda',39202);
[h,p] = adtest(Data10,'Distribution',distp10);

%fitting distributions
figure;
r101 = subplot(2,2,1);
histfit(Data10,11,'Normal')
Norm10 = mle(Data10, 'distribution', 'Normal');
xlabel('Parasite Count')
ylabel('Density')
title('Normal Distribution')
xlim([0,inf])
suptitle('Mouse RAG 10 Experiment 5')

r102 = subplot(2,2,2);
histfit(Data10,11,'Negative Binomial')
NegBin10 = mle(Data10, 'distribution', 'Negative Binomial');
xlabel('Parasite Count')
ylabel('Density')
title('Negative Binomial Distribution')
axis([0,187000,0,10])

r103 = subplot(2,2,3);
histfit(Data10,11,'Exponential')
Expo10 = mle(Data10, 'distribution', 'Exponential');
xlabel('Parasite Count')
ylabel('Density')
title('Exponential Distribution')
axis([0,187000,0,10])

r104 = subplot(2,2,4);
histfit(Data10,11,'Poisson')
Poi10 = mle(Data10, 'distribution', 'Poisson');
xlabel('Parasite Count')
ylabel('Density')
title('Poisson Distribution')

linkaxes([r101,r104],'xy')
