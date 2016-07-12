%% mouse A %%
A = [6242.308594 40213.92188 31013.05078 208242.6563 60042.10156 98883.60938 60476.59375 249944.2188 10892.03223 6057.03125	112096.5313	35391.50781	86337.54688	45613.80078	46020.12109	64171.25 56225.69141 20209.32813 55936.78906 42591.23438 88499.5 5941.76123	5517.262695	139638.9688];
DataA = round(reshape(A, [numel(A), 1]));

norA = fitdist(DataA, 'Normal');
xA = norA.NLogL; %log likelihood (normal dist.)

nbA = fitdist(DataA, 'Negative Binomial');
yA = nbA.NLogL; %log likelihood (Neg. Bin. dist.)

expA = fitdist(DataA, 'Exponential');
zA = expA.NLogL; %log likelihood (exponential dist.)

poiA = fitdist(DataA, 'Poisson');
wA = poiA.NLogL; %log likelihood (poisson dist.)

aicxA = 2*2-2*log(xA); %AIC = 2*#parameters - 2ln(ll)
aicyA = 2*2-2*log(yA);
aiczA = 2*1-2*log(zA);
aicwA = 2*1-2*log(wA);
aicA = aicbic([xA, yA, zA, wA],[2, 2, 1, 1]); %Akaike Information Criterion

%% mouse B %%
B = [113652.8672 23657.07813 7039.90332	67137.04688	64529.57031	131615.2969	518365.5938	541672.125	1800.070679	14785.81934	1394753	706560.6875	130131.8047	23352.04688	208169.5938	608231.375	106441.1953	32235.32422	44494.23438	94896.08594	6424.298828	5645.563477	27523.33594	32097.01172];
DataB = round(reshape(B, [numel(B), 1]));

norB = fitdist(DataB, 'Normal');
xB = norB.NLogL; %log likelihood (normal dist.)

nbB = fitdist(DataB, 'Negative Binomial');
yB = nbB.NLogL; %log likelihood (Neg. Bin. dist.)

expB = fitdist(DataB, 'Exponential');
zB = expB.NLogL;%log likelihood (exponential dist.)

poiB = fitdist(DataB, 'Poisson');
wB = poiB.NLogL; %log likelihood (poisson dist.)

aicxB = 2*2-2*log(xB); 
aicyB = 2*2-2*log(yB);
aiczB = 2*1-2*log(zB);
aicwB = 2*1-2*log(wB);
aicB = aicbic([xB, yB, zB, wB],[2, 2, 1, 1]); %Akaike Information Criterion

%% mouse C %%
C = [86988.09375 2483.02417	14647.36816	48430.07813	10765.91016	9871.273438	9871.273438	28048.60938	201290.25	13951.59961	42517.27734	971090.375	2863765.75	36039.0625	42434.35156	179498.3125	220141.0469	10395.6875	40059.78906	52496.25781	44046.84375	6026.121582	2963.452393	2269.223633];
DataC = round(reshape(C, [numel(C), 1]));

norC = fitdist(DataC, 'Normal');
xC = norC.NLogL; %log likelihood (normal dist.)

nbC = fitdist(DataC, 'Negative Binomial');
yC = nbC.NLogL; %log likelihood (Neg. Bin. dist.)

expC = fitdist(DataC, 'Exponential');
zC = expC.NLogL; %log likelihood (exponential dist.)

poiC = fitdist(DataC, 'Poisson');
wC = poiC.NLogL; %log likelihood (poisson dist.)

aicxC = 2*2-2*log(xC); 
aicyC = 2*2-2*log(yC);
aiczC = 2*1-2*log(zC);
aicwC = 2*1-2*log(wC);
aicC = aicbic([xC, yC, zC, wC],[2, 2, 1, 1]); %Akaike Information Criterion

%% mouse D %%
D = [3636.382568 9097.695313 3001.909668 23878.79297 918.046875	4182.071777	1520.401367	1117.206055	2155.180176	5595.969727	1236.572754	17303.19531	290.6083984	1136.777832	13201.63477	1235.423096	1470.355835	4630.273926	6996.893555	2384.648926	22086.73828	3055.782959	370.6409302	16257.79492];
DataD = round(reshape(D, [numel(D), 1]));

norD = fitdist(DataD, 'Normal');
xD = norD.NLogL; %log likelihood (normal dist.)

nbD = fitdist(DataD, 'Negative Binomial');
yD = nbD.NLogL; %log likelihood (Neg. Bin. dist.)

expD = fitdist(DataD, 'Exponential');
zD = expD.NLogL; %log likelihood (exponential dist.)

poiD = fitdist(DataD, 'Poisson');
wD = poiD.NLogL; %log likelihood (poisson dist.)

aicxD = 2*2-2*log(xD); 
aicyD = 2*2-2*log(yD);
aiczD = 2*1-2*log(zD);
aicwD = 2*1-2*log(wD);
aicD = aicbic([xD, yD, zD, wD],[2, 2, 1, 1]); %Akaike Information Criterion

%% mouse E %%
E = [15409.66016 38189.30078 7407.382813 12525.3877	14765.42969	48938.80859	11563.5166	16625.5625	6094.720703	3692.590332	41796.57813	14807.84375	9591.882813	3691.366699	39645.85156	9301.40625	29549.90625	12585.88477	3411.382324	41044.53906	2181.679932	5054.815918	25538.76953	3942.522705];
DataE = round(reshape(E, [numel(E), 1]));

norE = fitdist(DataE, 'Normal');
xE = norE.NLogL; %log likelihood (normal dist.)

nbE = fitdist(DataE, 'Negative Binomial');
yE = nbE.NLogL; %log likelihood (Neg. Bin. dist.)

expE = fitdist(DataE, 'Exponential');
zE = expE.NLogL; %log likelihood (exponential dist.)

poiE = fitdist(DataE, 'Poisson');
wE = poiE.NLogL; %log likelihood (poisson dist.)

aicxE = 2*2-2*log(xE); 
aicyE = 2*2-2*log(yE);
aiczE = 2*1-2*log(zE);
aicwE = 2*1-2*log(wE);
aicE = aicbic([xE, yE, zE, wE],[2, 2, 1, 1]); %Akaike Information Criterion

%% mouse F %%
F = [1918.327393 5368.237793 14975.54102 3805.636719 4451.460938 17625.02344 27372.53125 17554.89453 250.9739227 11183.27734 7323.999512 7664.599121 1906.455933 5972.289063 5060.821289 11847.68164 7574.185547 717.4968262 1669.182861 10734.91504 43123.44922 2413.856445 731.732666	42438.03906];
DataF = round(reshape(F, [numel(F), 1]));

norF = fitdist(DataF, 'Normal');
xF = norF.NLogL; %log likelihood (normal dist.)

nbF = fitdist(DataF, 'Negative Binomial');
yF = nbF.NLogL; %log likelihood (Neg. Bin. dist.)

expF = fitdist(DataF, 'Exponential');
zF = expF.NLogL; %log likelihood (exponential dist.)

poiF = fitdist(DataF, 'Poisson');
wF = poiF.NLogL; %log likelihood (poisson dist.)

aicxF = 2*2-2*log(xF); 
aicyF = 2*2-2*log(yF);
aiczF = 2*1-2*log(zF);
aicwF = 2*1-2*log(wF);
aicF = aicbic([xF, yF, zF, wF],[2, 2, 1, 1]); %Akaike Information Criterion

%% mouse G %%
G = [5879.079102 22170.35547 6642.881836 175939.375	29015.39844	8576.245117	15615.9082	20745.52734	1759.902832	2210.585449	8686.671875	7552.5625	4470.313965	4072.462891	10441.4375	5817.073242	12959.82227	4377.237305	11262.03516	14850.03125	1024.440186	2940.30249	6098.556641	6186.19043];
DataG = round(reshape(G, [numel(G), 1]));

norG = fitdist(DataG, 'Normal');
xG = norG.NLogL; %log likelihood (normal dist.)

nbG = fitdist(DataG, 'Negative Binomial');
yG = nbG.NLogL; %log likelihood (Neg. Bin. dist.)

expG = fitdist(DataG, 'Exponential');
zG = expG.NLogL; %log likelihood (exponential dist.)

poiG = fitdist(DataG, 'Poisson');
wG = poiG.NLogL; %log likelihood (poisson dist.)

aicxG = 2*2-2*log(xG); 
aicyG = 2*2-2*log(yG);
aiczG = 2*1-2*log(zG);
aicwG = 2*1-2*log(wG);
aicG = aicbic([xG, yG, zG, wG],[2, 2, 1, 1]); %Akaike Information Criterion

%% mouse H %%
H = [7838.777344 4216.459473 4139.382813 13826.62305 1222.324097 14711.82324 14711.82324 27832.98828 35677 2495.137695 5517.897461 13922.68457 9070.527344 2650.269531 8553.928711 17811.33398 7258.19043 15696.14648 2292.208496 7610.755859 50638.11719 11762.1748 1077.937134 32804.46484];
DataH = round(reshape(H, [numel(H), 1]));

norH = fitdist(DataH, 'Normal');
xH = norH.NLogL; %log likelihood (normal dist.)

nbH = fitdist(DataH, 'Negative Binomial');
yH = nbH.NLogL; %log likelihood (Neg. Bin. dist.)

expH = fitdist(DataH, 'Exponential');
zH = expH.NLogL; %log likelihood (exponential dist.)

poiH = fitdist(DataH, 'Poisson');
wH = poiH.NLogL; %log likelihood (poisson dist.)

aicxH = 2*2-2*log(xH);
aicyH = 2*2-2*log(yH);
aiczH = 2*1-2*log(zH);
aicwH = 2*1-2*log(wH);
aicH = aicbic([xH, yH, zH, wH],[2, 2, 1, 1]); %Akaike Information Criterion

%% mouse I %%
I = [1637.440186 1071.902344 1039.78186	1847.271973	1918.505127	696.3065186	3563.663086	2304.054688	2853.534668	674.1523438	3040.928223	322.5037231	2841.975586	509.4277344	514.2527466	1855.881958	816.8228149	1061.693604	1126.812866	1587.541748	695.0185547	215.0773773	571.5085449	7049.145996];
DataI = round(reshape(I, [numel(I), 1]));

norI = fitdist(DataI, 'Normal');
xI = norI.NLogL; %log likelihood (normal dist.)

nbI = fitdist(DataI, 'Negative Binomial');
yI = nbI.NLogL; %log likelihood (Neg. Bin. dist.)

expI = fitdist(DataI, 'Exponential');
zI = expI.NLogL; %log likelihood (exponential dist.)

poiI = fitdist(DataI, 'Poisson');
wI = poiI.NLogL; %log likelihood (poisson dist.)

aicxI = 2*2-2*log(xI);
aicyI = 2*2-2*log(yI);
aiczI = 2*1-2*log(zI);
aicwI = 2*1-2*log(wI);
aicI = aicbic([xI, yI, zI, wI],[2, 2, 1, 1]); %Akaike Information Criterion

%% mouse J %%
J = [15997.08789 3017.408203 5034.854004 9239.501953 2739.726074 513.6808472 1732.656006 1191.046387 1036.687012 3348.614746 1548.505249	860.7709961	7333.506836	1555.587158	2160.678223	3284.401367	3297.141357	922.8627319	4877.167969	6981.65332	555.0819092	763.1016846	27997.11719	45156.42969];
DataJ = round(reshape(J, [numel(J), 1]));

norJ = fitdist(DataJ, 'Normal');
xJ = norJ.NLogL; %log likelihood (normal dist.)

nbJ = fitdist(DataJ, 'Negative Binomial');
yJ = nbJ.NLogL; %log likelihood (Neg. Bin. dist.)

expJ = fitdist(DataJ, 'Exponential');
zJ = expJ.NLogL; %log likelihood (exponential dist.)

poiJ = fitdist(DataJ, 'Poisson');
wJ = poiJ.NLogL; %log likelihood (poisson dist.)

aicxJ = 2*2-2*log(xJ);
aicyJ = 2*2-2*log(yJ);
aiczJ = 2*1-2*log(zJ);
aicwJ = 2*1-2*log(wJ);
aicJ = aicbic([xJ, yJ, zJ, wJ],[2, 2, 1, 1]); %Akaike Information Criterion
