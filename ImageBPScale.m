i = imread('MouseC3scan2_19.tif'); % Read image
% figure(1);
% hold on; 

figure
imshow(i);
region = [2100	1100 5699 2415]; %This draws a 1uL size blood pool sampling site on your image 
% region = [2200 1000 1690 362]; %This draws a prob scale size blood pool
% sampling site on your image
% region = [1100 2674 1690 362]; %Prob
% region = [900 4348 1690 362];  %Prob
% region = [700 6022 1690 362];  %Prob
% region = [600 7696 1690 362];  %Prob
% region = [9000 300 362 1690];  %Prob
rectangle('Position', region, 'EdgeColor', 'r', 'LineWidth', 2);

cropped = imcrop(i, region);

red = cropped(:,:,1); % Red channel

t = 50; %This is the threshold, above which is considered 'red enough to be a parasite!'

% Use logical indexing to find pixels with >t intensity.
parasites = red>t;
figure
imshow(parasites); %This shows only the cropped area's 'red' pixels, i.e. parasites

NR = nnz(parasites); %Count number of just red pixels. This, in theory, could be equated to number of parasites... somehow!
title([num2str(NR) ' of ' num2str(numel(red)) ' pixels are red (' num2str(100*NR/numel(red)) '%, t = ' num2str(t) ')']);
