function [dataCell] = identifyNuclei(channel, diskElementOpen)
%identifyNuclei identifies cell nuclei in a single channel image
%   Inputs
%       channel: single channel image containing nuclei.
%       diskElementOpen: disk element used to remove noise
%   Outputs (as elements in the cell array dataCell):
%       numNuclei: the number of nuclei detected
%       chThreshClean: the input channel with background removed.
%       cellLabels:  a label matrix containing all the detected nuclei.

% Threshold the image to remove the background.
chAverage = mean(channel, 'all');
chSD = std(single(chAverage), 0, 'all');
chThresh = channel > (chAverage + chSD);
% Open the image to remove noise.
chThreshClean = imopen(chThresh, diskElementOpen); % remove small features

% Watershed to split nuclei that are overlapping and then count them.
distTrans = -bwdist(~chThreshClean); % distance transform
mask = imextendedmin(distTrans, 2); % removes noise from distance transform
distTrans2 = imimposemin(distTrans, mask);
basinMap = watershed(distTrans2);
cellArea = imdilate(chThreshClean, strel('disk', 3));
cellLabels = basinMap;
cellLabels(~cellArea) = 0;
numNuclei = max(cellLabels(:));

dataCell = {numNuclei, chThreshClean, cellLabels};
end

