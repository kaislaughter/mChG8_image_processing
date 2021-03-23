function [dataCells] = identifyFoci(chCorr, thresh, diskElementOpen)
%identifyFoci identifies and counts bright spots in an image
%   Inputs:
%       ch: a one-channel image (i.e., a matrix if pixel intensities) that
%           has the background removed and crosstalk corrections applied.
%       thresh: a numeric threshold for the channel
%       diskElementOpen: element used to remove noise.
%   Outputs (as elements of a cell array):
%       nFoci: the number of spots detected.
%       chThresh: the filtered and thresholded input channel.
%       chBinary: binarized image where foci have value of true
%       fociLabels: a label matrix containing all the detected foci.

% Threshold image so only bright spots remain.
chBinary = chCorr > thresh;
chBinaryClean = imopen(chBinary, diskElementOpen);
chThresh = chCorr;
chThresh(~chBinaryClean) = 0;

% Count foci.
fociLabels = watershed(~chBinaryClean);
fociLabels(~chBinaryClean) = 0;
nFoci = max(fociLabels(:));

dataCells = {nFoci, chThresh, chBinaryClean, fociLabels};
end

