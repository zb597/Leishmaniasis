%CD8 mice

%% Experimental ID C Experiment 2

MC = [86988.09	2483.024	14647.37	48430.08	10765.91	9871.273	9871.273	28048.61	201290.3	13951.6	42517.28	971090.4	2863766	36039.06	42434.35	179498.3	220141	10395.69	40059.79	52496.26	44046.84	6026.122	2963.452	2269.224];

%data
DataC = round(reshape(MC, [numel(MC), 1]));

%log likelihood
norC = fitdist(DataC, 'Normal');
xC = norC.NLogL; %log likelihood (normal dist.)

nbC = fitdist(DataC, 'Negative Binomial');
yC = nbC.NLogL; %log likelihood (Neg. Bin. dist.)

expC = fitdist(DataC, 'Exponential');
zC = expC.NLogL; %log likelihood (exponential dist.)

poiC = fitdist(DataC, 'Poisson');
wC = poiC.NLogL; %log likelihood (poisson dist.)

%Akaike Information Criterion
aicxC = 2*2-2*log(xC); %AIC = 2*#parameters - 2ln(ll)
aicyC = 2*2-2*log(yC);
aiczC = 2*1-2*log(zC);
aicwC = 2*1-2*log(wC);
aicC = aicbic([xC, yC, zC, wC],[2, 2, 1, 1]); 

%Anderson-Darling test
distnC = makedist('normal','mu',205837,'sigma',599685);
[h,p] = adtest(DataC,'Distribution',distnC);

distnbC = makedist('negative binomial','r',0.366567,'p',0.00000178086);
[h,p] = adtest(DataC,'Distribution',distnbC);

disteC = makedist('exponential','mu',205837);
[h,p] = adtest(DataC,'Distribution',disteC);

distpC = makedist('poisson','lambda',205837);
[h,p] = adtest(DataC,'Distribution',distpC);

%fitting distributions
figure;
c1 = subplot(3,1,1);
histfit(DataC,11,'Normal')
NormC = mle(DataC, 'distribution', 'Normal');
xlabel('Parasite Count')
ylabel('Density')
title('Normal Distribution')
xlim([0,inf])
ylim([0,25])
suptitle('Mouse C')

% c2 = subplot(2,2,2);
% histfit(DataC,11,'Negative Binomial')
% NegBinC = mle(DataC, 'distribution', 'Negative Binomial');
% xlabel('Parasite Count')
% ylabel('Density')
% title('Negative Binomial Distribution')
% axis([0,170,0,40])

c3 = subplot(3,1,2);
histfit(DataC,11,'Exponential')
ExpoC = mle(DataC, 'distribution', 'Exponential');
xlabel('Parasite Count')
ylabel('Density')
title('Exponential Distribution')

c4 = subplot(3,1,3);
histfit(DataC,11,'Poisson')
PoiC = mle(DataC, 'distribution', 'Poisson');
xlabel('Parasite Count')
ylabel('Density')
title('Poisson Distribution')

linkaxes([c1,c3,c4],'xy')

%% Experimental ID I Experiment 2

MI = [1637.44	1071.902	1039.782	1847.272	1918.505	696.3065	3563.663	2304.055	2853.535	674.1523	3040.928	322.5037	2841.976	509.4277	514.2527	1855.882	816.8228	1061.694	1126.813	1587.542	695.0186	215.0774	571.5085	7049.146];

%data
DataI = round(reshape(MI, [numel(MI), 1]));

%log likelihood
norI = fitdist(DataI, 'Normal');
xI = norI.NLogL; %log likelihood (normal dist.)

nbI = fitdist(DataI, 'Negative Binomial');
yI = nbI.NLogL; %log likelihood (Neg. Bin. dist.)

expI = fitdist(DataI, 'Exponential');
zI = expI.NLogL; %log likelihood (exponential dist.)

poiI = fitdist(DataI, 'Poisson');
wI = poiI.NLogL; %log likelihood (poisson dist.)

%Akaike Information Criterion
aicxI = 2*2-2*log(xI); %AIC = 2*#parameters - 2ln(ll)
aicyI = 2*2-2*log(yI);
aiczI = 2*1-2*log(zI);
aicwI = 2*1-2*log(wI);
aicI = aicbic([xI, yI, zI, wI],[2, 2, 1, 1]); 

%Anderson-Darling test
distnI = makedist('normal','mu',1659.04,'sigma',1482.8);
[h,p] = adtest(DataI,'Distribution',distnI);

distnbI = makedist('negative binomial','r',1.71863,'p',0.00103484);
[h,p] = adtest(DataI,'Distribution',distnbI);

disteI = makedist('exponential','mu',1659.04);
[h,p] = adtest(DataI,'Distribution',disteI);

distpI = makedist('poisson','lambda',1659.04);
[h,p] = adtest(DataI,'Distribution',distpI);

%fitting distributions
figure;
i1 = subplot(2,2,1);
histfit(DataI,11,'Normal')
NormI = mle(DataI, 'distribution', 'Normal');
xlabel('Parasite Count')
ylabel('Density')
title('Normal Distribution')
xlim([0,inf])
suptitle('Mouse I')

i2 = subplot(2,2,2);
histfit(DataI,11,'Negative Binomial')
NegBinI = mle(DataI, 'distribution', 'Negative Binomial');
xlabel('Parasite Count')
ylabel('Density')
title('Negative Binomial Distribution')
axis([0,7150,0,10])

i3 = subplot(2,2,3);
histfit(DataI,11,'Exponential')
ExpoI = mle(DataI, 'distribution', 'Exponential');
xlabel('Parasite Count')
ylabel('Density')
title('Exponential Distribution')
axis([0,7150,0,10])

i4 = subplot(2,2,4);
histfit(DataI,11,'Poisson')
PoiI = mle(DataI, 'distribution', 'Poisson');
xlabel('Parasite Count')
ylabel('Density')
title('Poisson Distribution')

linkaxes([i1,i4],'xy')

%% Experimental ID 4 Experiment 4

M4 = [98782.1	23432.9	547804.1	218793.6	30498.7	42836.4	132063.2	154165	22458.6	30615.8	68735.9	56174.5	13688	30541.9	123152.8	22664.2	166397.5	24192	19657.9	42091.4	4579.2	852.9	6435.1	4415.2];

%data
Data4 = round(reshape(M4, [numel(M4), 1]));

%log likelihood
nor4 = fitdist(Data4, 'Normal');
x4 = nor4.NLogL; %log likelihood (normal dist.)

nb4 = fitdist(Data4, 'Negative Binomial');
y4 = nb4.NLogL; %log likelihood (Neg. Bin. dist.)

exp4 = fitdist(Data4, 'Exponential');
z4 = exp4.NLogL; %log likelihood (exponential dist.)

poi4 = fitdist(Data4, 'Poisson');
w4 = poi4.NLogL; %log likelihood (poisson dist.)

%Akaike Information Criterion
aicx4 = 2*2-2*log(x4); %AIC = 2*#parameters - 2ln(ll)
aicy4 = 2*2-2*log(y4);
aicz4 = 2*1-2*log(z4);
aicw4 = 2*1-2*log(w4);
aic4 = aicbic([x4, y4, z4, w4],[2, 2, 1, 1]); 

%Anderson-Darling test
distn4 = makedist('normal','mu',78542.9,'sigma',116249);
[h,p] = adtest(Data4,'Distribution',distn4);

distnb4 = makedist('negative binomial','r',0.731727,'p',0.00000931618);
[h,p] = adtest(Data4,'Distribution',distnb4);

diste4 = makedist('exponential','mu',78542.9);
[h,p] = adtest(Data4,'Distribution',diste4);

distp4 = makedist('poisson','lambda',78542.9);
[h,p] = adtest(Data4,'Distribution',distp4);

%fitting distributions
figure;
r41 = subplot(2,2,1);
histfit(Data4,11,'Normal')
Norm4 = mle(Data4, 'distribution', 'Normal');
xlabel('Parasite Count')
ylabel('Density')
title('Normal Distribution')
xlim([0,inf])
suptitle('Mouse RAG 4 Experiment 4')

r42 = subplot(2,2,2);
histfit(Data4,11,'Negative Binomial')
NegBin4 = mle(Data4, 'distribution', 'Negative Binomial');
xlabel('Parasite Count')
ylabel('Density')
title('Negative Binomial Distribution')
axis([0,550000,0,15])

r43 = subplot(2,2,3);
histfit(Data4,11,'Exponential')
Expo4 = mle(Data4, 'distribution', 'Exponential');
xlabel('Parasite Count')
ylabel('Density')
title('Exponential Distribution')

r44 = subplot(2,2,4);
histfit(Data4,11,'Poisson')
Poi4 = mle(Data4, 'distribution', 'Poisson');
xlabel('Parasite Count')
ylabel('Density')
title('Poisson Distribution')

linkaxes([r41,r43,r44],'xy')

%% Experimental ID 5 Experiment 4

M54 = [181818.8	109668.8	412396.3	157043.2	10899.3	37303.5	52511.8	200723.1	75878.3	25746.7	16414.3	29792	96075.4	242080	28506	33780.8	46730.9	4845.2	19985.8	23996	16646.6	14755.1	4314.1	30883.5];

%data
Data54 = round(reshape(M54, [numel(M54), 1]));

%log likelihood
nor54 = fitdist(Data54, 'Normal');
x54 = nor54.NLogL; %log likelihood (normal dist.)

nb54 = fitdist(Data54, 'Negative Binomial');
y54 = nb54.NLogL; %log likelihood (Neg. Bin. dist.)

exp54 = fitdist(Data54, 'Exponential');
z54 = exp54.NLogL; %log likelihood (exponential dist.)

poi54 = fitdist(Data54, 'Poisson');
w54 = poi54.NLogL; %log likelihood (poisson dist.)

%Akaike Information Criterion
aicx54 = 2*2-2*log(x54); %AIC = 2*#parameters - 2ln(ll)
aicy54 = 2*2-2*log(y54);
aicz54 = 2*1-2*log(z54);
aicw54 = 2*1-2*log(w54);
aic54 = aicbic([x54, y54, z54, w54],[2, 2, 1, 1]); 

%Anderson-Darling test
distn54 = makedist('normal','mu',78033.2,'sigma',97857.4);
[h,p] = adtest(Data54,'Distribution',distn54);

distnb54 = makedist('negative binomial','r',0.897635,'p',0.0000115031);
[h,p] = adtest(Data54,'Distribution',distnb54);

diste54 = makedist('exponential','mu',78033.2);
[h,p] = adtest(Data54,'Distribution',diste54);

distp54 = makedist('poisson','lambda',78033.2);
[h,p] = adtest(Data54,'Distribution',distp54);

%fitting distributions
figure;
r541 = subplot(2,2,1);
histfit(Data54,11,'Normal')
Norm54 = mle(Data54, 'distribution', 'Normal');
xlabel('Parasite Count')
ylabel('Density')
title('Normal Distribution')
xlim([0,inf])
suptitle('Mouse RAG 5 Experiment 4')

r542 = subplot(2,2,2);
histfit(Data54,11,'Negative Binomial')
NegBin54 = mle(Data54, 'distribution', 'Negative Binomial');
xlabel('Parasite Count')
ylabel('Density')
title('Negative Binomial Distribution')
axis([0,418000,0,15])

r543 = subplot(2,2,3);
histfit(Data54,11,'Exponential')
Expo54 = mle(Data54, 'distribution', 'Exponential');
xlabel('Parasite Count')
ylabel('Density')
title('Exponential Distribution')
axis([0,418000,0,15])

r544 = subplot(2,2,4);
histfit(Data54,11,'Poisson')
Poi54 = mle(Data54, 'distribution', 'Poisson');
xlabel('Parasite Count')
ylabel('Density')
title('Poisson Distribution')

linkaxes([r541,r544],'xy')

%% Experimental ID 7 Experiment 4

M7 = [867842	552477.6	388486.9	1449490	164748.1	489848.7	374590.3	1272215	65345.5	16331.7	132568.7	218877.4	119276.6	184811.2	1266134	2730995	865704.1	3535952	5601397	4343749	5967.6	10030.2	764994.3	72690];

%data
Data7 = round(reshape(M7, [numel(M7), 1]));

%log likelihood
nor7 = fitdist(Data7, 'Normal');
x7 = nor7.NLogL; %log likelihood (normal dist.)

nb7 = fitdist(Data7, 'Negative Binomial');
y7 = nb7.NLogL; %log likelihood (Neg. Bin. dist.)

exp7 = fitdist(Data7, 'Exponential');
z7 = exp7.NLogL; %log likelihood (exponential dist.)

poi7 = fitdist(Data7, 'Poisson');
w7 = poi7.NLogL; %log likelihood (poisson dist.)

%Akaike Information Criterion
aicx7 = 2*2-2*log(x7); %AIC = 2*#parameters - 2ln(ll)
aicy7 = 2*2-2*log(y7);
aicz7 = 2*1-2*log(z7);
aicw7 = 2*1-2*log(w7);
aic7 = aicbic([x7, y7, z7, w7],[2, 2, 1, 1]); 

%Anderson-Darling test
distn7 = makedist('normal','mu',1062270,'sigma',1496670);
[h,p] = adtest(Data7,'Distribution',distn7);

distnb7 = makedist('negative binomial','r',0.55064,'p',0.000000518361);
[h,p] = adtest(Data7,'Distribution',distnb7);

diste7 = makedist('exponential','mu',1062270);
[h,p] = adtest(Data7,'Distribution',diste7);

distp7 = makedist('poisson','lambda',1062270);
[h,p] = adtest(Data7,'Distribution',distp7);

%fitting distributions
figure;
r71 = subplot(3,1,1);
histfit(Data7,11,'Normal')
Norm7 = mle(Data7, 'distribution', 'Normal');
xlabel('Parasite Count')
ylabel('Density')
title('Normal Distribution')
xlim([0,inf])
suptitle('Mouse RAG 7 Experiment 4')

% r72 = subplot(2,2,2);
% histfit(Data7,11,'Negative Binomial')
% NegBin7 = mle(Data7, 'distribution', 'Negative Binomial');
% xlabel('Parasite Count')
% ylabel('Density')
% title('Negative Binomial Distribution')
% axis([0,5610000,0,20])

r73 = subplot(3,1,2);
histfit(Data7,11,'Exponential')
Expo7 = mle(Data7, 'distribution', 'Exponential');
xlabel('Parasite Count')
ylabel('Density')
title('Exponential Distribution')
axis([0,5610000,0,20])

r74 = subplot(3,1,3);
histfit(Data7,11,'Poisson')
Poi7 = mle(Data7, 'distribution', 'Poisson');
xlabel('Parasite Count')
ylabel('Density')
title('Poisson Distribution')

linkaxes([r71,r74],'xy')

%% Experimental ID 9 Experiment 4

M9 = [46414.6	219530.5	52236.3	128603.4	62248.6	19690.2	22252.3	79189.2	11349.9	14101.1	12007.4	170496.8	55922.1	158562.8	456310.2	499857.7	332860.1	99346.4	216622.8	79926.2	9206.4	12650.2	3191.1	2906.4];

%data
Data9 = round(reshape(M9, [numel(M9), 1]));

%log likelihood
nor9 = fitdist(Data9, 'Normal');
x9 = nor9.NLogL; %log likelihood (normal dist.)

nb9 = fitdist(Data9, 'Negative Binomial');
y9 = nb9.NLogL; %log likelihood (Neg. Bin. dist.)

exp9 = fitdist(Data9, 'Exponential');
z9 = exp9.NLogL; %log likelihood (exponential dist.)

poi9 = fitdist(Data9, 'Poisson');
w9 = poi9.NLogL; %log likelihood (poisson dist.)

%Akaike Information Criterion
aicx9 = 2*2-2*log(x9); %AIC = 2*#parameters - 2ln(ll)
aicy9 = 2*2-2*log(y9);
aicz9 = 2*1-2*log(z9);
aicw9 = 2*1-2*log(w9);
aic9 = aicbic([x9, y9, z9, w9],[2, 2, 1, 1]); 

%Anderson-Darling test
distn9 = makedist('normal','mu',115228,'sigma',140204);
[h,p] = adtest(Data9,'Distribution',distn9);

distnb9 = makedist('negative binomial','r',0.721425,'p',0.00000626079);
[h,p] = adtest(Data9,'Distribution',distnb9);

diste9 = makedist('exponential','mu',115228);
[h,p] = adtest(Data9,'Distribution',diste9);

distp9 = makedist('poisson','lambda',115228);
[h,p] = adtest(Data9,'Distribution',distp9);

%fitting distributions
figure;
r91 = subplot(2,2,1);
histfit(Data9,11,'Normal')
Norm9 = mle(Data9, 'distribution', 'Normal');
xlabel('Parasite Count')
ylabel('Density')
title('Normal Distribution')
axis([0,506000,0,10])
suptitle('Mouse RAG 9 Experiment 4')

r92 = subplot(2,2,2);
histfit(Data9,11,'Negative Binomial')
NegBin9 = mle(Data9, 'distribution', 'Negative Binomial');
xlabel('Parasite Count')
ylabel('Density')
title('Negative Binomial Distribution')
axis([0,506000,0,10])

r93 = subplot(2,2,3);
histfit(Data9,11,'Exponential')
Expo9 = mle(Data9, 'distribution', 'Exponential');
xlabel('Parasite Count')
ylabel('Density')
title('Exponential Distribution')
axis([0,506000,0,10])

r94 = subplot(2,2,4);
histfit(Data9,11,'Poisson')
Poi9 = mle(Data9, 'distribution', 'Poisson');
xlabel('Parasite Count')
ylabel('Density')
title('Poisson Distribution')
axis([0,506000,0,10])

%% Experimental ID 2 Experiment 5

M2 = [53771.23	15734.01	12067.08	35575.91	28530.69	12811.91	10118.98	13885.84	7139.758	16043.8	9844.204	6656.57	14711.64	56160.36	132494	48295.63	87786.71	68436.54	36531.82	72108.79	12688.16	14848.2	8130.676	10732.51];

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
distn2 = makedist('normal','mu',32712.8,'sigma',31852.1);
[h,p] = adtest(Data2,'Distribution',distn2);

distnb2 = makedist('negative binomial','r',1.42469,'p',0.0000435496);
[h,p] = adtest(Data2,'Distribution',distnb2);

diste2 = makedist('exponential','mu',32712.8);
[h,p] = adtest(Data2,'Distribution',diste2);

distp2 = makedist('poisson','lambda',32712.8);
[h,p] = adtest(Data2,'Distribution',distp2);

%fitting distributions
figure;
r21 = subplot(2,2,1);
histfit(Data2,11,'Normal')
Norm2 = mle(Data2, 'distribution', 'Normal');
xlabel('Parasite Count')
ylabel('Density')
title('Normal Distribution')
xlim([0,inf])
suptitle('Mouse RAG 2 Experiment 5')

r22 = subplot(2,2,2);
histfit(Data2,11,'Negative Binomial')
NegBin2 = mle(Data2, 'distribution', 'Negative Binomial');
xlabel('Parasite Count')
ylabel('Density')
title('Negative Binomial Distribution')
axis([0,143000,0,10])

r23 = subplot(2,2,3);
histfit(Data2,11,'Exponential')
Expo2 = mle(Data2, 'distribution', 'Exponential');
xlabel('Parasite Count')
ylabel('Density')
title('Exponential Distribution')
axis([0,143000,0,10])

r24 = subplot(2,2,4);
histfit(Data2,11,'Poisson')
Poi2 = mle(Data2, 'distribution', 'Poisson');
xlabel('Parasite Count')
ylabel('Density')
title('Poisson Distribution')

linkaxes([r21,r24],'xy')

%% Experimental ID 5 Experiment 5

M55 = [2099.204	849.3257	3339.803	2380.775	3682.561	1291.839	5476.861	1360.699	1888.724	4222.9	2481.203	1476.645	3769.24	1287.345	1493.629	1151.788	2676.89	3870.346	3575.469	6870.077	901.0186	1427.349	1226.367	1000.991];

%data
Data55 = round(reshape(M55, [numel(M55), 1]));

%log likelihood
nor55 = fitdist(Data55, 'Normal');
x55 = nor55.NLogL; %log likelihood (normal dist.)

nb55 = fitdist(Data55, 'Negative Binomial');
y55 = nb55.NLogL; %log likelihood (Neg. Bin. dist.)

exp55 = fitdist(Data55, 'Exponential');
z55 = exp55.NLogL; %log likelihood (exponential dist.)

poi55 = fitdist(Data55, 'Poisson');
w55 = poi55.NLogL; %log likelihood (poisson dist.)

%Akaike Information Criterion
aicx55 = 2*2-2*log(x55); %AIC = 2*#parameters - 2ln(ll)
aicy55 = 2*2-2*log(y55);
aicz55 = 2*1-2*log(z55);
aicw55 = 2*1-2*log(w55);
aic55 = aicbic([x55, y55, z55, w55],[2, 2, 1, 1]); 

%Anderson-Darling test
distn55 = makedist('normal','mu',2491.71,'sigma',1569.58);
[h,p] = adtest(Data55,'Distribution',distn55);

distnb55 = makedist('negative binomial','r',3.0176,'p',0.00120959);
[h,p] = adtest(Data55,'Distribution',distnb55);

diste55 = makedist('exponential','mu',2491.71);
[h,p] = adtest(Data55,'Distribution',diste55);

distp55 = makedist('poisson','lambda',2491.71);
[h,p] = adtest(Data55,'Distribution',distp55);

%fitting distributions
figure;
r551 = subplot(2,2,1);
histfit(Data55,11,'Normal')
Norm55 = mle(Data55, 'distribution', 'Normal');
xlabel('Parasite Count')
ylabel('Density')
title('Normal Distribution')
axis([0,6880,0,10])
suptitle('Mouse RAG 5 Experiment 5')

r552 = subplot(2,2,2);
histfit(Data55,11,'Negative Binomial')
NegBin55 = mle(Data55, 'distribution', 'Negative Binomial');
xlabel('Parasite Count')
ylabel('Density')
title('Negative Binomial Distribution')
axis([0,6880,0,10])

r553 = subplot(2,2,3);
histfit(Data55,11,'Exponential')
Expo55 = mle(Data55, 'distribution', 'Exponential');
xlabel('Parasite Count')
ylabel('Density')
title('Exponential Distribution')
axis([0,6880,0,10])

r554 = subplot(2,2,4);
histfit(Data55,11,'Poisson')
Poi55 = mle(Data55, 'distribution', 'Poisson');
xlabel('Parasite Count')
ylabel('Density')
title('Poisson Distribution')
axis([0,6880,0,10])

%% Experimental ID 8 Experiment 5

M8 = [28611.06	3129.208	3214.273	4881.184	3551.307	17673.08	47348.29	15514.55	860.3441	3260.179	64844.03	3.775513	17518.26	19609.43	3.540197	7105.975	6733.547	9861.456	11472.41	11072.78	4183.498	14446.12	9036.298	11701.71];

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

%Anderson-Darling test
distn8 = makedist('normal','mu',13151.4,'sigma',15206.6);
[h,p] = adtest(Data8,'Distribution',distn8);

distnb8 = makedist('negative binomial','r',0.614947,'p',0.0000467568);
[h,p] = adtest(Data8,'Distribution',distnb8);

diste8 = makedist('exponential','mu',13151.4);
[h,p] = adtest(Data8,'Distribution',diste8);

distp8 = makedist('poisson','lambda',13151.4);
[h,p] = adtest(Data8,'Distribution',distp8);

%fitting distributions
figure;
r81 = subplot(2,2,1);
histfit(Data8,11,'Normal')
Norm8 = mle(Data8, 'distribution', 'Normal');
xlabel('Parasite Count')
ylabel('Density')
title('Normal Distribution')
xlim([0,inf])
suptitle('Mouse RAG 8 Experiment 5')

r82 = subplot(2,2,2);
histfit(Data8,11,'Negative Binomial')
NegBin8 = mle(Data8, 'distribution', 'Negative Binomial');
xlabel('Parasite Count')
ylabel('Density')
title('Negative Binomial Distribution')
axis([0,64900,0,10])

r83 = subplot(2,2,3);
histfit(Data8,11,'Exponential')
Expo8 = mle(Data8, 'distribution', 'Exponential');
xlabel('Parasite Count')
ylabel('Density')
title('Exponential Distribution')
axis([0,64900,0,10])

r84 = subplot(2,2,4);
histfit(Data8,11,'Poisson')
Poi8 = mle(Data8, 'distribution', 'Poisson');
xlabel('Parasite Count')
ylabel('Density')
title('Poisson Distribution')

linkaxes([r81,r84],'xy')
