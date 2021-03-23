function [dataCells] = Manders_ED(ch1, ch2)
%Manders_ED calculates the Manders correlation coefficient between two
%thresholded images.
%   Both ch1 and ch2 should already be background-subtracted, for example
%   using a top hat filter and despeckling/opening.
products = (double(ch1).*double(ch2));
ch1sq = double(ch1).^2;
ch2sq = double(ch2).^2;

% Calculate coefficients
rP = corr2(ch1, ch2);
rOverlap = sum(products(:)) / (sum(ch1sq(:))*sum(ch2sq(:)));
MOC = sum(products(:)) / sqrt(sum(ch1sq(:)) * sum(ch2sq(:)));


ch1Sum = sum(ch1(:));
ch2Sum = sum(ch2(:));

% Calculate fractional overlaps
ch1OverlapRaw = ch1(ch2>0);
ch2OverlapRaw = ch2(ch1>0);
ch1Overlap = sum(ch1OverlapRaw(:)) / ch1Sum;
ch2Overlap = sum(ch2OverlapRaw(:)) / ch2Sum;

% Construct an array of results to return
dataCells = {rP, rOverlap, MOC, ch1Overlap, ch2Overlap};
end

