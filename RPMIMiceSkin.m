%RPMI mice

%% Experimental ID F Experiment 2

F = [1918.327	5368.238	14975.54	3805.637	4451.461	17625.02	27372.53	17554.89	250.9739	11183.28	7324	7664.599	1906.456	5972.289	5060.821	11847.68	7574.186	717.4968	1669.183	10734.92	43123.45	2413.856	731.7327	42438.04];

%data
DataF = round(reshape(F, [numel(F), 1]));

%log likelihood
norF = fitdist(DataF, 'Normal');
xF = norF.NLogL; %log likelihood (normal dist.)

nbF = fitdist(DataF, 'Negative Binomial');
yF = nbF.NLogL; %log likelihood (Neg. Bin. dist.)

expF = fitdist(DataF, 'Exponential');
zF = expF.NLogL; %log likelihood (exponential dist.)

poiF = fitdist(DataF, 'Poisson');
wF = poiF.NLogL; %log likelihood (poisson dist.)

%Akaike Information Criterion
aicxF = 2*2-2*log(xF); %AIC = 2*#parameters - 2ln(ll)
aicyF = 2*2-2*log(yF);
aiczF = 2*1-2*log(zF);
aicwF = 2*1-2*log(wF);
aicF = aicbic([xF, yF, zF, wF],[2, 2, 1, 1]); 

%Anderson-Darling test
distnF = makedist('normal','mu',10570.2,'sigma',11911.5);
[h,p] = adtest(DataF,'Distribution',distnF);

distnbF = makedist('negative binomial','r',0.90809,'p',0.0000859033);
[h,p] = adtest(DataF,'Distribution',distnbF);

disteF = makedist('exponential','mu',10570.2);
[h,p] = adtest(DataF,'Distribution',disteF);

distpF = makedist('poisson','lambda',10570.2);
[h,p] = adtest(DataF,'Distribution',distpF);

%fitting distributions
figure;
f1 = subplot(2,2,1);
histfit(DataF,11,'Normal')
NormF = mle(DataF, 'distribution', 'Normal');
xlabel('Parasite Count')
ylabel('Density')
title('Normal Distribution')
axis([0,44000,0,10])
suptitle('Mouse F')

f2 = subplot(2,2,2);
histfit(DataF,11,'Negative Binomial')
NegBinF = mle(DataF, 'distribution', 'Negative Binomial');
xlabel('Parasite Count')
ylabel('Density')
title('Negative Binomial Distribution')
axis([0,44000,0,10])

f3 = subplot(2,2,3);
histfit(DataF,11,'Exponential')
ExpoF = mle(DataF, 'distribution', 'Exponential');
xlabel('Parasite Count')
ylabel('Density')
title('Exponential Distribution')
axis([0,44000,0,10])

f4 = subplot(2,2,4);
histfit(DataF,11,'Poisson')
PoiF = mle(DataF, 'distribution', 'Poisson');
xlabel('Parasite Count')
ylabel('Density')
title('Poisson Distribution')
axis([0,44000,0,10])

%% Experimental ID J Experiment 2

J = [15997.09	3017.408	5034.854	9239.502	2739.726	513.6808	1732.656	1191.046	1036.687	3348.615	1548.505	860.771	7333.507	1555.587	2160.678	3284.401	3297.141	922.8627	4877.168	6981.653	555.0819	763.1017	27997.12	45156.43];

%data
DataJ = round(reshape(J, [numel(J), 1]));

%log likelihood
norJ = fitdist(DataJ, 'Normal');
xJ = norJ.NLogL; %log likelihood (normal dist.)

nbJ = fitdist(DataJ, 'Negative Binomial');
yJ = nbJ.NLogL; %log likelihood (Neg. Bin. dist.)

expJ = fitdist(DataJ, 'Exponential');
zJ = expJ.NLogL; %log likelihood (exponential dist.)

poiJ = fitdist(DataJ, 'Poisson');
wJ = poiJ.NLogL; %log likelihood (poisson dist.)

%Akaike Information Criterion
aicxJ = 2*2-2*log(xJ); %AIC = 2*#parameters - 2ln(ll)
aicyJ = 2*2-2*log(yJ);
aiczJ = 2*1-2*log(zJ);
aicwJ = 2*1-2*log(wJ);
aicJ = aicbic([xJ, yJ, zJ, wJ],[2, 2, 1, 1]); 

%Anderson-Darling test
distnJ = makedist('normal','mu',6297.83,'sigma',10281.2);
[h,p] = adtest(DataJ,'Distribution',distnJ);

distnbJ = makedist('negative binomial','r',0.772314,'p',0.000122617);
[h,p] = adtest(DataJ,'Distribution',distnbJ);

disteJ = makedist('exponential','mu',6297.83);
[h,p] = adtest(DataJ,'Distribution',disteJ);

distpJ = makedist('poisson','lambda',6297.83);
[h,p] = adtest(DataJ,'Distribution',distpJ);

%fitting distributions
figure;
j1 = subplot(2,2,1);
histfit(DataJ,11,'Normal')
NormJ = mle(DataJ, 'distribution', 'Normal');
xlabel('Parasite Count')
ylabel('Density')
title('Normal Distribution')
xlim([0,inf])
suptitle('Mouse J')

j2 = subplot(2,2,2);
histfit(DataJ,11,'Negative Binomial')
NegBinJ = mle(DataJ, 'distribution', 'Negative Binomial');
xlabel('Parasite Count')
ylabel('Density')
title('Negative Binomial Distribution')
axis([0,46100,0,20])

j3 = subplot(2,2,3);
histfit(DataJ,11,'Exponential')
ExpoJ = mle(DataJ, 'distribution', 'Exponential');
xlabel('Parasite Count')
ylabel('Density')
title('Exponential Distribution')

j4 = subplot(2,2,4);
histfit(DataJ,11,'Poisson')
PoiJ = mle(DataJ, 'distribution', 'Poisson');
xlabel('Parasite Count')
ylabel('Density')
title('Poisson Distribution')

linkaxes([j1,j3,j4],'xy')

%% Experimental ID 3 Experiment 4

M34 = [79325.4	8075.2	12407.7	1621.6	12520.3	14747.8	6971.3	4854.6	2778.3	1358.6	14703.1	5527.3	3123.6	3629.4	5529.5	4353.3	6067.5	2363.5	4163.9	2053.3	5635.2	3689.9	8057.8	12618.4];

%data
Data34 = round(reshape(M34, [numel(M34), 1]));

%log likelihood
nor34 = fitdist(Data34, 'Normal');
x34 = nor34.NLogL; %log likelihood (normal dist.)

nb34 = fitdist(Data34, 'Negative Binomial');
y34 = nb34.NLogL; %log likelihood (Neg. Bin. dist.)

exp34 = fitdist(Data34, 'Exponential');
z34 = exp34.NLogL; %log likelihood (exponential dist.)

poi34 = fitdist(Data34, 'Poisson');
w34 = poi34.NLogL; %log likelihood (poisson dist.)

%Akaike Information Criterion
aicx34 = 2*2-2*log(x34); %AIC = 2*#parameters - 2ln(ll)
aicy34 = 2*2-2*log(y34);
aicz34 = 2*1-2*log(z34);
aicw34 = 2*1-2*log(w34);
aic34 = aicbic([x34, y34, z34, w34],[2, 2, 1, 1]); 

%Anderson-Darling test
distn34 = makedist('normal','mu',9424.04,'sigma',15451.7);
[h,p] = adtest(Data34,'Distribution',distn34);

distnb34 = makedist('negative binomial','r',1.15228,'p',0.000122256);
[h,p] = adtest(Data34,'Distribution',distnb34);

diste34 = makedist('exponential','mu',9424.04);
[h,p] = adtest(Data34,'Distribution',diste34);

distp34 = makedist('poisson','lambda',9424.04);
[h,p] = adtest(Data34,'Distribution',distp34);

%fitting distributions
figure;
r341 = subplot(2,2,1);
histfit(Data34,11,'Normal')
Norm34 = mle(Data34, 'distribution', 'Normal');
xlabel('Parasite Count')
ylabel('Density')
title('Normal Distribution')
xlim([0,inf])
suptitle('Mouse RAG 3 Experiment 4')

r342 = subplot(2,2,2);
histfit(Data34,11,'Negative Binomial')
NegBin34 = mle(Data34, 'distribution', 'Negative Binomial');
xlabel('Parasite Count')
ylabel('Density')
title('Negative Binomial Distribution')
axis([0,80300,0,20])

r343 = subplot(2,2,3);
histfit(Data34,11,'Exponential')
Expo34 = mle(Data34, 'distribution', 'Exponential');
xlabel('Parasite Count')
ylabel('Density')
title('Exponential Distribution')

r344 = subplot(2,2,4);
histfit(Data34,11,'Poisson')
Poi34 = mle(Data34, 'distribution', 'Poisson');
xlabel('Parasite Count')
ylabel('Density')
title('Poisson Distribution')

linkaxes([r341,r343,r344],'xy')

%% Experimental ID 10 Experiment 4

M10 = [4243310	5552395	4087400	4627846	5546594	17248300	17070900	9188658	3691940	25452820	34302760	23473170	7529109	9002134	4844421	3201581	3572419	5662281	9909828	4986509	1713397	193055.3	593888.8	833192.3];

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
distn10 = makedist('normal','mu',8605330,'sigma',8713650);
[h,p] = adtest(Data10,'Distribution',distn10);

% distnb10 = makedist('negative binomial','r',1.08474,'p',0.000000126054);
% [h,p] = adtest(Data10,'Distribution',distnb10);

diste10 = makedist('exponential','mu',8605330);
[h,p] = adtest(Data10,'Distribution',diste10);

distp10 = makedist('poisson','lambda',8605330);
[h,p] = adtest(Data10,'Distribution',distp10);

%fitting distributions
figure;
r101 = subplot(3,1,1);
histfit(Data10,11,'Normal')
Norm10 = mle(Data10, 'distribution', 'Normal');
xlabel('Parasite Count')
ylabel('Density')
title('Normal Distribution')
axis([0,35200000,0,15])
suptitle('Mouse RAG 10 Experiment 4')

% r102 = subplot(2,2,2);
% histfit(Data10,11,'Negative Binomial')
% NegBin10 = mle(Data10, 'distribution', 'Negative Binomial');
% xlabel('Parasite Count')
% ylabel('Density')
% title('Negative Binomial Distribution')
% axis([0,35200000,0,15])

r103 = subplot(3,1,2);
histfit(Data10,11,'Exponential')
Expo10 = mle(Data10, 'distribution', 'Exponential');
xlabel('Parasite Count')
ylabel('Density')
title('Exponential Distribution')
axis([0,35200000,0,15])

r104 = subplot(3,1,3);
histfit(Data10,11,'Poisson')
Poi10 = mle(Data10, 'distribution', 'Poisson');
xlabel('Parasite Count')
ylabel('Density')
title('Poisson Distribution')
axis([0,35200000,0,15])

%% Experimental ID 3 Experiment 5

M35 = [32331.17	11941.36	30553.34	23279.64	148774.2	335843	42148.09	38192.13	12640.68	31603.28	103105.3	68400.34	76280.03	174505.1	40626.54	108187.5	23332.35	43019.24	114974.5	98845.32	6164.186	13454.91	26539.47	11808.43];

%data
Data35 = round(reshape(M35, [numel(M35), 1]));

%log likelihood
nor35 = fitdist(Data35, 'Normal');
x35 = nor35.NLogL; %log likelihood (normal dist.)

nb35 = fitdist(Data35, 'Negative Binomial');
y35 = nb35.NLogL; %log likelihood (Neg. Bin. dist.)

exp35 = fitdist(Data35, 'Exponential');
z35 = exp35.NLogL; %log likelihood (exponential dist.)

poi35 = fitdist(Data35, 'Poisson');
w35 = poi35.NLogL; %log likelihood (poisson dist.)

%Akaike Information Criterion
aicx35 = 2*2-2*log(x35); %AIC = 2*#parameters - 2ln(ll)
aicy35 = 2*2-2*log(y35);
aicz35 = 2*1-2*log(z35);
aicw35 = 2*1-2*log(w35);
aic35 = aicbic([x35, y35, z35, w35],[2, 2, 1, 1]); 

%Anderson-Darling test
distn35 = makedist('normal','mu',67356.2,'sigma',73567.2);
[h,p] = adtest(Data35,'Distribution',distn35);

distnb35 = makedist('negative binomial','r',1.21461,'p',0.0000180324);
[h,p] = adtest(Data35,'Distribution',distnb35);

diste35 = makedist('exponential','mu',67356.2);
[h,p] = adtest(Data35,'Distribution',diste35);

distp35 = makedist('poisson','lambda',67356.2);
[h,p] = adtest(Data35,'Distribution',distp35);

%fitting distributions
figure;
r351 = subplot(2,2,1);
histfit(Data35,11,'Normal')
Norm35 = mle(Data35, 'distribution', 'Normal');
xlabel('Parasite Count')
ylabel('Density')
title('Normal Distribution')
xlim([0,inf])
suptitle('Mouse RAG 3 Experiment 5')

r352 = subplot(2,2,2);
histfit(Data35,11,'Negative Binomial')
NegBin35 = mle(Data35, 'distribution', 'Negative Binomial');
xlabel('Parasite Count')
ylabel('Density')
title('Negative Binomial Distribution')
axis([0,341000,0,10])

r353 = subplot(2,2,3);
histfit(Data35,11,'Exponential')
Expo35 = mle(Data35, 'distribution', 'Exponential');
xlabel('Parasite Count')
ylabel('Density')
title('Exponential Distribution')
axis([0,341000,0,10])

r354 = subplot(2,2,4);
histfit(Data35,11,'Poisson')
Poi35 = mle(Data35, 'distribution', 'Poisson');
xlabel('Parasite Count')
ylabel('Density')
title('Poisson Distribution')

linkaxes([r351,r354],'xy')

%% Experimental ID 4 Experiment 5

M4 = [6816.114	9865.926	8074.566	2144.396	3790.682	4670.691	287609.7	17519.6	675.1464	2809.901	4721.435	17965.03	10342.02	15399.22	1696.771	16272.53	3538.723	1732.837	4501.902	630.3101	6376.592	79432.82	1688.43	4939.643];

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
distn4 = makedist('normal','mu',21384,'sigma',58858.9);
[h,p] = adtest(Data4,'Distribution',distn4);

distnb4 = makedist('negative binomial','r',0.501679,'p',0.0000234599);
[h,p] = adtest(Data4,'Distribution',distnb4);

diste4 = makedist('exponential','mu',21384);
[h,p] = adtest(Data4,'Distribution',diste4);

distp4 = makedist('poisson','lambda',21384);
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
suptitle('Mouse RAG 4 Experiment 5')

r42 = subplot(2,2,2);
histfit(Data4,11,'Negative Binomial')
NegBin4 = mle(Data4, 'distribution', 'Negative Binomial');
xlabel('Parasite Count')
ylabel('Density')
title('Negative Binomial Distribution')
axis([0,297000,0,30])

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

%% Experimental ID 7 Experiment 5

M7 = [8176.658	4045.592	11200.34	27846.53	7319.354	3848.755	2961.245	10853.75	13593.43	82828.38	5796.723	14738.95	44362.54	2780.33	11222.72	18019.24	31724.3	4964.867	7930.374	8426.884	40682.51	12424.21	6102.104	2921.042];

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
distn7 = makedist('normal','mu',16032.1,'sigma',18360.1);
[h,p] = adtest(Data7,'Distribution',distn7);

distnb7 = makedist('negative binomial','r',1.27581,'p',0.0000795722);
[h,p] = adtest(Data7,'Distribution',distnb7);

diste7 = makedist('exponential','mu',16032.1);
[h,p] = adtest(Data7,'Distribution',diste7);

distp7 = makedist('poisson','lambda',16032.1 );
[h,p] = adtest(Data7,'Distribution',distp7);

%fitting distributions
figure;
r71 = subplot(2,2,1);
histfit(Data7,11,'Normal')
Norm7 = mle(Data7, 'distribution', 'Normal');
xlabel('Parasite Count')
ylabel('Density')
title('Normal Distribution')
xlim([0,inf])
suptitle('Mouse RAG 7 Experiment 5')

r72 = subplot(2,2,2);
histfit(Data7,11,'Negative Binomial')
NegBin7 = mle(Data7, 'distribution', 'Negative Binomial');
xlabel('Parasite Count')
ylabel('Density')
title('Negative Binomial Distribution')
axis([0,83600,0,10])

r73 = subplot(2,2,3);
histfit(Data7,11,'Exponential')
Expo7 = mle(Data7, 'distribution', 'Exponential');
xlabel('Parasite Count')
ylabel('Density')
title('Exponential Distribution')

r74 = subplot(2,2,4);
histfit(Data7,11,'Poisson')
Poi7 = mle(Data7, 'distribution', 'Poisson');
xlabel('Parasite Count')
ylabel('Density')
title('Poisson Distribution')

linkaxes([r71,r73,r74],'xy')

%% Experimental ID 9 Experiment 5

M9 = [94325.42	62468.18	71573.48	59310.95	74108.88	103655.9	51661.73	193353.3	83582.15	84971.95	27241.44	274852.2	253882.5	50821.75	153677.9	296190.4	44208.68	353988.8	340442.5	128945	28428.64	28618.89	6126.26	7119.881];

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
distn9 = makedist('normal','mu',119732,'sigma',107149);
[h,p] = adtest(Data9,'Distribution',distn9);

distnb9 = makedist('negative binomial','r',1.24623,'p',0.0000104085);
[h,p] = adtest(Data9,'Distribution',distnb9);

diste9 = makedist('exponential','mu',119732);
[h,p] = adtest(Data9,'Distribution',diste9);

distp9 = makedist('poisson','lambda',119732);
[h,p] = adtest(Data9,'Distribution',distp9);

%fitting distributions
figure;
r91 = subplot(2,2,1);
histfit(Data9,11,'Normal')
Norm9 = mle(Data9, 'distribution', 'Normal');
xlabel('Parasite Count')
ylabel('Density')
title('Normal Distribution')
axis([0,363000,0,6])
suptitle('Mouse RAG 9 Experiment 5')

r92 = subplot(2,2,2);
histfit(Data9,11,'Negative Binomial')
NegBin9 = mle(Data9, 'distribution', 'Negative Binomial');
xlabel('Parasite Count')
ylabel('Density')
title('Negative Binomial Distribution')
axis([0,363000,0,6])

r93 = subplot(2,2,3);
histfit(Data9,11,'Exponential')
Expo9 = mle(Data9, 'distribution', 'Exponential');
xlabel('Parasite Count')
ylabel('Density')
title('Exponential Distribution')
axis([0,363000,0,6])

r94 = subplot(2,2,4);
histfit(Data9,11,'Poisson')
Poi9 = mle(Data9, 'distribution', 'Poisson');
xlabel('Parasite Count')
ylabel('Density')
title('Poisson Distribution')
axis([0,363000,0,6])

%% Experimental ID 11 Experiment 5

M11 = [43841.88	29074.17	13066.45	81767.38	119878.5	33898.7	165113.4	100725.1	73412.82	147603.8	58429.43	247062.3	122996.5	214556.3	62988.18	116733.1	114132	179094.4	68986.91	118917.3	35624.51	3932.968	12931.35	70823.13];

%data
Data11 = round(reshape(M11, [numel(M11), 1]));

%log likelihood
nor11 = fitdist(Data11, 'Normal');
x11 = nor11.NLogL; %log likelihood (normal dist.)

nb11 = fitdist(Data11, 'Negative Binomial');
y11 = nb11.NLogL; %log likelihood (Neg. Bin. dist.)

exp11 = fitdist(Data11, 'Exponential');
z11 = exp11.NLogL; %log likelihood (exponential dist.)

poi11 = fitdist(Data11, 'Poisson');
w11 = poi11.NLogL; %log likelihood (poisson dist.)

%Akaike Information Criterion
aicx11 = 2*2-2*log(x11); %AIC = 2*#parameters - 2ln(ll)
aicy11 = 2*2-2*log(y11);
aicz11 = 2*1-2*log(z11);
aicw11 = 2*1-2*log(w11);
aicPre11 = aicbic([x11, y11, z11, w11],[2, 2, 1, 1]); 

%Anderson-Darling test
distn11 = makedist('normal','mu',93149.5,'sigma',64323.7);
[h,p] = adtest(Data11,'Distribution',distn11);

distnb11 = makedist('negative binomial','r',1.64396,'p',0.0000176483);
[h,p] = adtest(Data11,'Distribution',distnb11);

diste11 = makedist('exponential','mu',93149.5);
[h,p] = adtest(Data11,'Distribution',diste11);

distp11 = makedist('poisson','lambda',93149.5);
[h,p] = adtest(Data11,'Distribution',distp11);

%fitting distributions
figure;
r111 = subplot(2,2,1);
histfit(Data11,11,'Normal')
Norm11 = mle(Data11, 'distribution', 'Normal');
xlabel('Parasite Count')
ylabel('Density')
title('Normal Distribution')
axis([0,253000,0,5])
suptitle('Mouse RAG 11 Experiment 5')

r112 = subplot(2,2,2);
histfit(Data11,11,'Negative Binomial')
NegBin11 = mle(Data11, 'distribution', 'Negative Binomial');
xlabel('Parasite Count')
ylabel('Density')
title('Negative Binomial Distribution')
axis([0,253000,0,5])

r113 = subplot(2,2,3);
histfit(Data11,11,'Exponential')
Expo11 = mle(Data11, 'distribution', 'Exponential');
xlabel('Parasite Count')
ylabel('Density')
title('Exponential Distribution')
axis([0,253000,0,5])

r114 = subplot(2,2,4);
histfit(Data11,11,'Poisson')
Poi11 = mle(Data11, 'distribution', 'Poisson');
xlabel('Parasite Count')
ylabel('Density')
title('Poisson Distribution')
axis([0,253000,0,5])
