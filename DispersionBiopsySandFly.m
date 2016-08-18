%% Dispersion parameters
%macro scale dispersion parameter
MSDP = [0.90809 0.772314 1.15228 1.08474 1.21461 0.501679 1.27581 1.24623 1.64396 1.19602 0.533042 0.840364 0.486133 0.535471 0.685018 1.03875 1.9186 1.5164 0.366567 0.0949068 0.731727 0.897635 0.55064 0.721425 1.42469 3.0176 0.614947];
%parasites in a sand fly dispersion parameter
SFDP = [0.0175009 0.041603 0.0594375 0.345845 0.11377 0.0329398 0.0372163 0.134039 0.205923 0.109374 0.0475553 0.17425 0.0794412 0.322317 0.542371 0.125262 0.0482573 0.0408038 0.0657911 0.0492123 0.123687 0.107996 0.402308 0.0948693 0.106422 0.079551 0.138334];
scatter(MSDP,SFDP)
xlabel('Skin Biopsy Data')
ylabel('Sand Fly Data')
title('Dispersion Parameter')

%% Sand fly means against dispersion parameter for skin data
SandFlyMeans = [0.325	1.95	129	66478	2524.8	202.675	134.85	261.925	1599.1	7671.4	1733.7	64.625	795.525	21763	16172	320.15	114	161.1	3604.8	15.9	410.825	794.475	46686	172.9	140	18.8	0.4];
MSDP = [0.90809 0.772314 1.15228 1.08474 1.21461 0.501679 1.27581 1.24623 1.64396 1.19602 0.533042 0.840364 0.486133 0.535471 0.685018 1.03875 1.9186 1.5164 0.366567 0.0949068 0.731727 0.897635 0.55064 0.721425 1.42469 3.0176 0.614947];
scatter(MSDP,SandFlyMeans)
xlabel('Dispersion Parameter For Skin Biopsy Data')
ylabel('Mean Numbers of Parasites in a Sand Fly')
title('Dispersion of Parasites in the Skin Against the Average Number of Parasites in a Sand Fly')

%% Skin biopsy means against dispersion parameter for sand fly
SkinBiopMeans = [10570	6298	9424	8605300	67356	21384	16032	119730	93150	65675	204380	105260	414860	3065700	8794500	71932	5284	39202	205840	1659	78543	78033	1062300	115230	32713	2492	13151];
SFDP = [0.0175009 0.041603 0.0594375 0.345845 0.11377 0.0329398 0.0372163 0.134039 0.205923 0.109374 0.0475553 0.17425 0.0794412 0.322317 0.542371 0.125262 0.0482573 0.0408038 0.0657911 0.0492123 0.123687 0.107996 0.402308 0.0948693 0.106422 0.079551 0.138334];
scatter(SFDP,SkinBiopMeans)
xlabel('Dispersion Parameter for Sand Fly Data')
ylabel('Means Numbers of Parasites in the Skin')
title('Dispersion of Parasites in the Sand Fly Against the Average Number of Parasites in the Skin')