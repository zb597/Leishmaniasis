Biop = [7838.777344 4216.459473 4139.382813 13826.62305 1222.324097 14711.82324 14711.82324 27832.98828 35677 2495.137695 5517.897461 13922.68457 9070.527344 2650.269531 8553.928711 17811.33398 7258.19043 15696.14648 2292.208496 7610.755859 50638.11719 11762.1748 1077.937134 32804.46484];
Data = round(reshape(Biop, [numel(Biop), 1]));

%% maximum likelihood estimates %%
% nor = fitdist(Data, 'Normal');
% x = nor.NLogL;
% 
% nb = fitdist(Data, 'Negative Binomial');
% y = nb.NLogL;
% 
% exp = fitdist(Data, 'Exponential');
% z = exp.NLogL;
% 
% poi = fitdist(Data, 'Poisson');
% w = poi.NLogL; %log likelihood
% 
% aic = aicbic([x, y, z, w],[2, 2, 1, 1]);

%%

% nor = fitdist(Data, 'Normal');
% x = 0:1:50000;
% y = pdf(nor, x);
% plot(x, y)
% hold on
% histfit(Data,10,'Normal')
% axis([0,60000,0,8])

% histogram(Data,10)
% Norm = mle(Data, 'distribution', 'Normal')

% nb = fitdist(Data, 'Negative Binomial');
% x = 0:1:50000;
% y = pdf(nb, x);
% plot(x, y)
% hold on
% histfit(Data,10,'Negative Binomial')
% axis([0,60000,0,8])
% 
% histogram(Data,10)
% NegBin = mle(Data, 'distribution', 'Negative Binomial')

% exp = fitdist(Data, 'Exponential');
% x = 0:1:60000;
% y = pdf(exp, x);
% plot(x, y)
% hold on
% histfit(Data,10,'Exponential')
% axis([0,60000,0,10])
% 
% histogram(Data,10)
% Expo = mle(Data, 'distribution', 'Exponential')

% poi = fitdist(Data, 'Poisson');
% x = 0:1:60000;
% y = pdf(poi,x);
% plot(x, y)
% hold on
% histfit(Data,10,'Poisson')
% axis([0,60000,0,10])
% 
% histogram(Data,10)
% Poi = mle(Data, 'distribution', 'Poisson')

%%
% distn = makedist('normal','mu',13056,'sigma',12129);
% [h,p] = adtest(Data,'Distribution',distn);
% 
% distnb = makedist('negative binomial','r',1.2702,'p',0.0001);
% [h,p] = adtest(Data,'Distribution',distnb);
% 
% diste = makedist('exponential','mu',1.3056);
% [h,p] = adtest(Data,'Distribution',diste);
% 
% distp = makedist('poisson','lambda',1.3056);
% [h,p] = adtest(Data,'Distribution',distp);

%%
% subplot(2,2,1)
% N = fitdist(Data,'Normal');
% x_values1 = -25000:1:50000;
% y1 = pdf(N,x_values1);
% plot(x_values1,y1)
% 
% subplot(2,2,2)
% E = fitdist(Data,'Exponential');
% x_values2 = 0:1:50000;
% y2 = pdf(E,x_values2);
% plot(x_values2,y2)
% 
% subplot(2,2,3)
% P = fitdist(Data,'Poisson');
% x_values3 = 10000:1:11000;
% y3 = pdf(P,x_values3);
% plot(x_values3,y3)
% 
% subplot(2,2,4)
% NB = fitdist(Data,'Negative Binomial');
% x_values4 = 0:1:50000;
% y4 = pdf(NB,x_values4);
% plot(x_values4,y4)