function [dataCells] = identifyFoci(ch, diskElementTH, diskElementOpen)
%identifyFoci identifies and counts bright spots in an image
%   Inputs:
%       ch: a one-channel image (i.e., a matrix if pixel intensities).
%       diskElementTH: element used with a top hat filter to remove
%                      background.
%       diskElementOpen: element used to remove noise.
%   Outputs (as elements of a cell array):
%       nFoci: the number of spots detected.
%       chThresh: the filtered and thresholded input channel.
%       chBinary: binarized image where foci have value of true
%       fociLabels: a label matrix containing all the detected foci.

% Adaptive background removal using a top hat filter.
chFiltered = imtophat(ch, diskElementTH);
% Threshold image so only bright spots remain.
chAverage = mean(chFiltered, 'all');
chSD = std(single(chFiltered), 0, 'all'); % std doesn't accept uint16 data
chBinary = chFiltered > (chAverage + chSD);
chBinaryClean = imopen(chBinary, diskElementOpen);
chThresh = chFiltered;
chThresh(~chBinaryClean) = 0;

% Count foci.
fociLabels = watershed(~chBinaryClean);
fociLabels(~chBinaryClean) = 0;
nFoci = max(fociLabels(:));

dataCells = {nFoci, chThresh, chBinaryClean, fociLabels};
end

