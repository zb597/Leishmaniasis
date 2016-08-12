%CD8 mice

%% Untreated Rag 1

Rag1 = [1734.9	4606.1	4255.6	5866.0	8419.2	1115.6	915.1	813.4	1905.9	7063.1	915.2	1245.8	1602.3	8424.4	1362.0	2523.7	1140.7	8686.8	12205.9	2051.1	3327.6	1340.7	980.6	173.5];

%data
Data1 = round(reshape(Rag1, [numel(Rag1), 1]));

%log likelihood
nor1 = fitdist(Data1, 'Normal');
x1 = nor1.NLogL; %log likelihood (normal dist.)

nb1 = fitdist(Data1, 'Negative Binomial');
y1 = nb1.NLogL; %log likelihood (Neg. Bin. dist.)

exp1 = fitdist(Data1, 'Exponential');
z1 = exp1.NLogL; %log likelihood (exponential dist.)

poi1 = fitdist(Data1, 'Poisson');
w1 = poi1.NLogL; %log likelihood (poisson dist.)

%Akaike Information Criterion
aicx1 = 2*2-2*log(x1); %AIC = 2*#parameters - 2ln(ll)
aicy1 = 2*2-2*log(y1);
aicz1 = 2*1-2*log(z1);
aicw1 = 2*1-2*log(w1);
aic1 = aicbic([x1, y1, z1, w1],[2, 2, 1, 1]); 

%Anderson-Darling test
distn1 = makedist('normal','mu',3444.88,'sigma',3277.39);
[h,p] = adtest(Data1,'Distribution',distn1);

distnb1 = makedist('negative binomial','r',1.26968,'p',0.000368435);
[h,p] = adtest(Data1,'Distribution',distnb1);

diste1 = makedist('exponential','mu',3444.88);
[h,p] = adtest(Data1,'Distribution',diste1);

distp1 = makedist('poisson','lambda',3444.88);
[h,p] = adtest(Data1,'Distribution',distp1);

%fitting distributions
figure;
c1 = subplot(2,2,1);
histfit(Data1,11,'Normal')
Norm1 = mle(Data1, 'distribution', 'Normal');
xlabel('Parasite Count')
ylabel('Density')
title('Normal Distribution')
xlim([0,inf])
suptitle('Rag 1')

c2 = subplot(2,2,2);
histfit(Data1,11,'Negative Binomial')
NegBin1 = mle(Data1, 'distribution', 'Negative Binomial');
xlabel('Parasite Count')
ylabel('Density')
title('Negative Binomial Distribution')
axis([0,13200,0,10])

c3 = subplot(2,2,3);
histfit(Data1,11,'Exponential')
Expo1 = mle(Data1, 'distribution', 'Exponential');
xlabel('Parasite Count')
ylabel('Density')
title('Exponential Distribution')
axis([0,13200,0,10])

c4 = subplot(2,2,4);
histfit(Data1,11,'Poisson')
Poi1 = mle(Data1, 'distribution', 'Poisson');
xlabel('Parasite Count')
ylabel('Density')
title('Poisson Distribution')

linkaxes([c1,c4],'xy')

%% Untreated Rag 2

Rag2 = [5581.7	1476.5	756.9	2162.6	16668.0	422.8	1995.1	11084.2	2355.6	6956.3	221.8	2708.2	1839.7	1675.0	2105.4	737.8	249.1	2177.8	470.5	167.1	3667.1	1235.0	656.6	3548.2];

%data
Data2 = round(reshape(Rag2, [numel(Rag2), 1]));

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
distn2 = makedist('normal','mu',2955.04,'sigma',3841.35);
[h,p] = adtest(Data2,'Distribution',distn2);

distnb2 = makedist('negative binomial','r',0.915594,'p',0.000309745);
[h,p] = adtest(Data2,'Distribution',distnb2);

diste2 = makedist('exponential','mu',2955.04);
[h,p] = adtest(Data2,'Distribution',diste2);

distp2 = makedist('poisson','lambda',2955.04);
[h,p] = adtest(Data2,'Distribution',distp2);

%fitting distributions
figure;
i1 = subplot(2,2,1);
histfit(Data2,11,'Normal')
Norm2 = mle(Data2, 'distribution', 'Normal');
xlabel('Parasite Count')
ylabel('Density')
title('Normal Distribution')
xlim([0,inf])
suptitle('Rag 2')

i2 = subplot(2,2,2);
histfit(Data2,11,'Negative Binomial')
NegBin2 = mle(Data2, 'distribution', 'Negative Binomial');
xlabel('Parasite Count')
ylabel('Density')
title('Negative Binomial Distribution')
axis([0,17600,0,10])

i3 = subplot(2,2,3);
histfit(Data2,11,'Exponential')
Expo2 = mle(Data2, 'distribution', 'Exponential');
xlabel('Parasite Count')
ylabel('Density')
title('Exponential Distribution')
axis([0,17600,0,10])

i4 = subplot(2,2,4);
histfit(Data2,11,'Poisson')
Poi2 = mle(Data2, 'distribution', 'Poisson');
xlabel('Parasite Count')
ylabel('Density')
title('Poisson Distribution')

linkaxes([i1,i4],'xy')

%% Untreated Rag 3

Rag3 = [505.5	173.2	1884.3	5646.4	1433.9	1635.7	800.3	1670.8	16082.9	3466.9	1376.7	524.0	1751.0	4847.2	3003.9	1316.0	1356.7	2526.7	2671.0	5941.0	7748.1	6120.4	3849.7	714.5];

%data
Data3 = round(reshape(Rag3, [numel(Rag3), 1]));

%log likelihood
nor3 = fitdist(Data3, 'Normal');
x3 = nor3.NLogL; %log likelihood (normal dist.)

nb3 = fitdist(Data3, 'Negative Binomial');
y3 = nb3.NLogL; %log likelihood (Neg. Bin. dist.)

exp3 = fitdist(Data3, 'Exponential');
z3 = exp3.NLogL; %log likelihood (exponential dist.)

poi3 = fitdist(Data3, 'Poisson');
w3 = poi3.NLogL; %log likelihood (poisson dist.)

%Akaike Information Criterion
aicx3 = 2*2-2*log(x3); %AIC = 2*#parameters - 2ln(ll)
aicy3 = 2*2-2*log(y3);
aicz3 = 2*1-2*log(z3);
aicw3 = 2*1-2*log(w3);
aic3 = aicbic([x3, y3, z3, w3],[2, 2, 1, 1]); 

%Anderson-Darling test
distn3 = makedist('normal','mu',3210.33,'sigma',3427.34);
[h,p] = adtest(Data3,'Distribution',distn3);

distnb3 = makedist('negative binomial','r',1.24736,'p',0.000388394);
[h,p] = adtest(Data3,'Distribution',distnb3);

diste3 = makedist('exponential','mu',3210.33);
[h,p] = adtest(Data3,'Distribution',diste3);

distp3 = makedist('poisson','lambda',3210.33);
[h,p] = adtest(Data3,'Distribution',distp3);

%fitting distributions
figure;
r31 = subplot(2,2,1)
histfit(Data3,11,'Normal')
Norm3 = mle(Data3, 'distribution', 'Normal');
xlabel('Parasite Count')
ylabel('Density')
title('Normal Distribution')
xlim([0,inf])
suptitle('Rag 3')

r32 = subplot(2,2,2);
histfit(Data3,11,'Negative Binomial')
NegBin3 = mle(Data3, 'distribution', 'Negative Binomial');
xlabel('Parasite Count')
ylabel('Density')
title('Negative Binomial Distribution')
axis([0,16500,0,10])

r33 = subplot(2,2,3);
histfit(Data3,11,'Exponential')
Expo3 = mle(Data3, 'distribution', 'Exponential');
xlabel('Parasite Count')
ylabel('Density')
title('Exponential Distribution')
axis([0,16500,0,10])

r34 = subplot(2,2,4);
histfit(Data3,11,'Poisson')
Poi3 = mle(Data3, 'distribution', 'Poisson');
xlabel('Parasite Count')
ylabel('Density')
title('Poisson Distribution')

linkaxes([r31,r34],'xy')

%% Untreated Rag 4

Rag4 = [1082.3	841.0	787.7	680.5	1473.9	66.4	578.4	521.7	755.7	529.4	781.4	781.1	1660.2	1170.1	1677.1	567.9	1291.2	1171.7	1797.6	700.2	2031.9	1318.6	1630.4	1447.0];

%data
Data4 = round(reshape(Rag4, [numel(Rag4), 1]));

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
distn4 = makedist('normal','mu',1055.96,'sigma',497.755);
[h,p] = adtest(Data4,'Distribution',distn4);

distnb4 = makedist('negative binomial','r',3.29603,'p',0.00311165);
[h,p] = adtest(Data4,'Distribution',distnb4);

diste4 = makedist('exponential','mu',1055.96);
[h,p] = adtest(Data4,'Distribution',diste4);

distp4 = makedist('poisson','lambda',1055.96);
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
suptitle('Rag 4')

r42 = subplot(2,2,2);
histfit(Data4,11,'Negative Binomial')
NegBin4 = mle(Data4, 'distribution', 'Negative Binomial');
xlabel('Parasite Count')
ylabel('Density')
title('Negative Binomial Distribution')
axis([0,2090,0,4])

r43 = subplot(2,2,3);
histfit(Data4,11,'Exponential')
Expo4 = mle(Data4, 'distribution', 'Exponential');
xlabel('Parasite Count')
ylabel('Density')
title('Exponential Distribution')
axis([0,2090,0,4])

r44 = subplot(2,2,4);
histfit(Data4,11,'Poisson')
Poi4 = mle(Data4, 'distribution', 'Poisson');
xlabel('Parasite Count')
ylabel('Density')
title('Poisson Distribution')

linkaxes([r41,r44],'xy')

%% Untreated Rag 5

Rag5 = [868.3	303.7	2543.6	248.4	12076.7	98.8	246.4	413.2	1042.1	685.7	792.8	1656.2	1659.3	717.2	2846.7	2748.0	623.6	564.2	917.4	2703.1	1084.0	1620.8	402.9	129.7];

%data
Data5 = round(reshape(Rag5, [numel(Rag5), 1]));

%log likelihood
nor5 = fitdist(Data5, 'Normal');
x5 = nor5.NLogL; %log likelihood (normal dist.)

nb5 = fitdist(Data5, 'Negative Binomial');
y5 = nb5.NLogL; %log likelihood (Neg. Bin. dist.)

exp5 = fitdist(Data5, 'Exponential');
z5 = exp5.NLogL; %log likelihood (exponential dist.)

poi5 = fitdist(Data5, 'Poisson');
w5 = poi5.NLogL; %log likelihood (poisson dist.)

%Akaike Information Criterion
aicx5 = 2*2-2*log(x5); %AIC = 2*#parameters - 2ln(ll)
aicy5 = 2*2-2*log(y5);
aicz5 = 2*1-2*log(z5);
aicw5 = 2*1-2*log(w5);
aic5 = aicbic([x5, y5, z5, w5],[2, 2, 1, 1]); 

%Anderson-Darling test
distn5 = makedist('normal','mu',1541.38,'sigma',2405.98);
[h,p] = adtest(Data5,'Distribution',distn5);

distnb5 = makedist('negative binomial','r',0.952277,'p',0.000617428);
[h,p] = adtest(Data5,'Distribution',distnb5);

diste5 = makedist('exponential','mu',1541.38);
[h,p] = adtest(Data5,'Distribution',diste5);

distp5 = makedist('poisson','lambda',1541.38);
[h,p] = adtest(Data5,'Distribution',distp5);

%fitting distributions
figure;
r51 = subplot(2,2,1);
histfit(Data5,11,'Normal')
Norm5 = mle(Data5, 'distribution', 'Normal');
xlabel('Parasite Count')
ylabel('Density')
title('Normal Distribution')
xlim([0,inf])
suptitle('Rag 5')

r52 = subplot(2,2,2);
histfit(Data5,11,'Negative Binomial')
NegBin5 = mle(Data5, 'distribution', 'Negative Binomial');
xlabel('Parasite Count')
ylabel('Density')
title('Negative Binomial Distribution')
axis([0,12100,0,20])

r53 = subplot(2,2,3);
histfit(Data5,11,'Exponential')
Expo5 = mle(Data5, 'distribution', 'Exponential');
xlabel('Parasite Count')
ylabel('Density')
title('Exponential Distribution')
axis([0,12100,0,20])

r54 = subplot(2,2,4);
histfit(Data5,11,'Poisson')
Poi5 = mle(Data5, 'distribution', 'Poisson');
xlabel('Parasite Count')
ylabel('Density')
title('Poisson Distribution')

linkaxes([r51,r54],'xy')

%% Untreated Rag 6

Rag6 = [18815.3	61973.8	37314.1	2225.0	144193.4	28475.9	122334.3	15205.1	6420.3	4582.2	11234.9	47808.0	15934.2	3369.7	22128.1	5687.2	7504.0	1921.5	3141.7	1870.7	34220.9	5318.9	899.5	48252.0];

%data
Data6 = round(reshape(Rag6, [numel(Rag6), 1]));

%log likelihood
nor6 = fitdist(Data6, 'Normal');
x6 = nor6.NLogL; %log likelihood (normal dist.)

nb6 = fitdist(Data6, 'Negative Binomial');
y6 = nb6.NLogL; %log likelihood (Neg. Bin. dist.)

exp6 = fitdist(Data6, 'Exponential');
z6 = exp6.NLogL; %log likelihood (exponential dist.)

poi6 = fitdist(Data6, 'Poisson');
w6 = poi6.NLogL; %log likelihood (poisson dist.)

%Akaike Information Criterion
aicx6 = 2*2-2*log(x6); %AIC = 2*#parameters - 2ln(ll)
aicy6 = 2*2-2*log(y6);
aicz6 = 2*1-2*log(z6);
aicw6 = 2*1-2*log(w6);
aic6 = aicbic([x6, y6, z6, w6],[2, 2, 1, 1]); 

%Anderson-Darling test
distn6 = makedist('normal','mu',27118,'sigma',37069.4);
[h,p] = adtest(Data6,'Distribution',distn6);

distnb6 = makedist('negative binomial','r',0.722932,'p',0.0000266581);
[h,p] = adtest(Data6,'Distribution',distnb6);

diste6 = makedist('exponential','mu',27118);
[h,p] = adtest(Data6,'Distribution',diste6);

distp6 = makedist('poisson','lambda',27118);
[h,p] = adtest(Data6,'Distribution',distp6);

%fitting distributions
figure;
r61 = subplot(2,2,1);
histfit(Data6,11,'Normal')
Norm6 = mle(Data6, 'distribution', 'Normal');
xlabel('Parasite Count')
ylabel('Density')
title('Normal Distribution')
axis([0,154000,0,12])
suptitle('Rag 6')

r62 = subplot(2,2,2);
histfit(Data6,11,'Negative Binomial')
NegBin6 = mle(Data6, 'distribution', 'Negative Binomial');
xlabel('Parasite Count')
ylabel('Density')
title('Negative Binomial Distribution')
axis([0,154000,0,12])

r63 = subplot(2,2,3);
histfit(Data6,11,'Exponential')
Expo6 = mle(Data6, 'distribution', 'Exponential');
xlabel('Parasite Count')
ylabel('Density')
title('Exponential Distribution')
axis([0,154000,0,12])

r64 = subplot(2,2,4);
histfit(Data6,11,'Poisson')
Poi6 = mle(Data6, 'distribution', 'Poisson');
xlabel('Parasite Count')
ylabel('Density')
title('Poisson Distribution')
axis([0,154000,0,12])

%% Untreated Rag 7

Rag7 = [81094.5	397965.5	95943.1	189913.2	54714.1	53272.7	820945.0	1489117.0	76049.7	133271.6	374300.6	943975.9	19880.9	78146.6	369051.2	237573.3	23747.6	262029.2	433661.6	1001859.9	1429.9	15965.9	15016.3	69149.0];

%data
Data7 = round(reshape(Rag7, [numel(Rag7), 1]));

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
distn7 = makedist('normal','mu',301587,'sigma',387002);
[h,p] = adtest(Data7,'Distribution',distn7);

distnb7 = makedist('negative binomial','r',0.656076,'p',0.00000217541);
[h,p] = adtest(Data7,'Distribution',distnb7);

diste7 = makedist('exponential','mu',301587);
[h,p] = adtest(Data7,'Distribution',diste7);

distp7 = makedist('poisson','lambda',301587);
[h,p] = adtest(Data7,'Distribution',distp7);

%fitting distributions
figure;
r71 = subplot(2,2,1);
histfit(Data7,11,'Normal')
Norm7 = mle(Data7, 'distribution', 'Normal');
xlabel('Parasite Count')
ylabel('Density')
title('Normal Distribution')
axis([0,1540000,0,15])
suptitle('Rag 7')

r72 = subplot(2,2,2);
histfit(Data7,11,'Negative Binomial')
NegBin7 = mle(Data7, 'distribution', 'Negative Binomial');
xlabel('Parasite Count')
ylabel('Density')
title('Negative Binomial Distribution')
axis([0,1540000,0,15])

r73 = subplot(2,2,3);
histfit(Data7,11,'Exponential')
Expo7 = mle(Data7, 'distribution', 'Exponential');
xlabel('Parasite Count')
ylabel('Density')
title('Exponential Distribution')
axis([0,1540000,0,15])

r74 = subplot(2,2,4);
histfit(Data7,11,'Poisson')
Poi7 = mle(Data7, 'distribution', 'Poisson');
xlabel('Parasite Count')
ylabel('Density')
title('Poisson Distribution')
axis([0,1540000,0,15])

%% Untreated Rag 8

Rag8 = [266994.7	612286.0	587857.0	578527.3	37708.9	192540.2	204166.2	349983.6	4876.1	61548.5	228418.5	54736.7	77599.3	107644.2	159764.8	34644.1	37963.3	15894.1	154566.9	132574.3	5059.9	2870.5	157096.4	98691.8];

%data
Data8 = round(reshape(Rag8, [numel(Rag8), 1]));

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
distn8 = makedist('normal','mu',173501,'sigma',185049);
[h,p] = adtest(Data8,'Distribution',distn8);

distnb8 = makedist('negative binomial','r',0.812449,'p',0.00000468266);
[h,p] = adtest(Data8,'Distribution',distnb8);

diste8 = makedist('exponential','mu',173501);
[h,p] = adtest(Data8,'Distribution',diste8);

distp8 = makedist('poisson','lambda',173501);
[h,p] = adtest(Data8,'Distribution',distp8);

%fitting distributions
figure;
r81 = subplot(2,2,1);
histfit(Data8,11,'Normal')
Norm8 = mle(Data8, 'distribution', 'Normal');
xlabel('Parasite Count')
ylabel('Density')
title('Normal Distribution')
axis([0,616000,0,8])
suptitle('Rag 8')

r82 = subplot(2,2,2);
histfit(Data8,11,'Negative Binomial')
NegBin8 = mle(Data8, 'distribution', 'Negative Binomial');
xlabel('Parasite Count')
ylabel('Density')
title('Negative Binomial Distribution')
axis([0,616000,0,8])

r83 = subplot(2,2,3);
histfit(Data8,11,'Exponential')
Expo8 = mle(Data8, 'distribution', 'Exponential');
xlabel('Parasite Count')
ylabel('Density')
title('Exponential Distribution')
axis([0,616000,0,8])

r84 = subplot(2,2,4);
histfit(Data8,11,'Poisson')
Poi8 = mle(Data8, 'distribution', 'Poisson');
xlabel('Parasite Count')
ylabel('Density')
title('Poisson Distribution')
axis([0,616000,0,8])

%% Untreated Rag 9

Rag9 = [12413.2	10824.3	3444.4	8597.4	703.2	20057.6	82067.7	12291.2	16895.1	22202.2	26330.9	2902.4	63653.7	41884.9	77456.6	6736.4	5400.5	24222.7	35832.7	16875.7	318.0	572.1	996.2	2235.6];

%data
Data9 = round(reshape(Rag9, [numel(Rag9), 1]));

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
distn9 = makedist('normal','mu',20621.5,'sigma',23740.1);
[h,p] = adtest(Data9,'Distribution',distn9);

distnb9 = makedist('negative binomial','r',0.722228,'p',0.0000350219);
[h,p] = adtest(Data9,'Distribution',distnb9);

diste9 = makedist('exponential','mu',20621.5);
[h,p] = adtest(Data9,'Distribution',diste9);

distp9 = makedist('poisson','lambda',20621.5);
[h,p] = adtest(Data9,'Distribution',distp9);

%fitting distributions
figure;
r91 = subplot(2,2,1);
histfit(Data9,11,'Normal')
Norm9 = mle(Data9, 'distribution', 'Normal');
xlabel('Parasite Count')
ylabel('Density')
title('Normal Distribution')
xlim([0,inf])
suptitle('Rag 9')

r92 = subplot(2,2,2);
histfit(Data9,11,'Negative Binomial')
NegBin9 = mle(Data9, 'distribution', 'Negative Binomial');
xlabel('Parasite Count')
ylabel('Density')
title('Negative Binomial Distribution')
axis([0,82500,0,10])

r93 = subplot(2,2,3);
histfit(Data9,11,'Exponential')
Expo9 = mle(Data9, 'distribution', 'Exponential');
xlabel('Parasite Count')
ylabel('Density')
title('Exponential Distribution')
axis([0,82500,0,10])

r94 = subplot(2,2,4);
histfit(Data9,11,'Poisson')
Poi9 = mle(Data9, 'distribution', 'Poisson');
xlabel('Parasite Count')
ylabel('Density')
title('Poisson Distribution')

linkaxes([r91,r94],'xy')
