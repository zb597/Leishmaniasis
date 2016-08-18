SkinBiopMeans = [10570	6298	9424	8605300	67356	21384	16032	119730	93150	65675	204380	105260	414860	3065700	8794500	71932	5284	39202	205840	1659	78543	78033	1062300	115230	32713	2492	13151];
SandFlyMeans = [0.325	1.95	129	66478	2524.8	202.675	134.85	261.925	1599.1	7671.4	1733.7	64.625	795.525	21763	16172	320.15	114	161.1	3604.8	15.9	410.825	794.475	46686	172.9	140	18.8	0.4];

%%
figure;
scatter(SkinBiopMeans,SandFlyMeans,'x')
xlabel('Means of Skin Biopsy')
ylabel('Means of Sand Fly ')
title('Means of parasite numbers in a skin biopsy against in a sand fly vector')
hold on
for i=1:27
    scatter(SkinBiopMeans(i),SandFlyMeans(i),'r+')
    %label the points with the corresponding 'x' value
    text(SkinBiopMeans(i),SandFlyMeans(i),num2str(i));
end
%xlim([0,50000])

%%
figure;
plot(SkinBiopMeans,'+')
xlabel('Mouse')
ylabel('Mean number of parasites')
hold on
plot(SandFlyMeans,'*')
legend('Skin Biopsies','Sandfly Vector')
mouse = {'F';'J';'M34';'M104';'M35';'M45';'M75';'M95';'M115';'A';'B';'M14';'M24';'M64';'M84';'M15';'M65';'M105';'C';'I';'M44';'M54';'M74';'M94';'M25';'M55';'M85'};
ax = gca;
set(ax,'XTickLabel',mouse)
title('Mean numbers of parasites on different scales')

%%
figure;
subplot(2,1,1)
plot(SkinBiopMeans,'+')
% mouse = ['F';'J';'M34';'M104';'M35';'M45';'M75';'M95';'M115';'A';'B';'M14';'M24';'M64';'M84';'M15';'M65';'M105';'C';'I';'M44';'M54';'M74';'M94';'M25';'M55';'M85'];
% ax = gca;
% set(ax,'XTicklabel',mouse)
xlabel('Mouse')
ylabel('Skin biopsy')
subplot(2,1,2)
plot(SandFlyMeans,'*')
% mouse = ['A';'B';'C';'D';'E';'F';'G';'H';'I';'J'];
% ax = gca;
% set(ax,'XTicklabel',mouse)
xlabel('Mouse')
ylabel('Sandfly vector')
suptitle('Mean number of parasites')
