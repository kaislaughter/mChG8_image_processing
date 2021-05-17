function [SNR] = Gal8SNRHelper(config, PCData, NCData, testThresh)
%Gal8SNRHelper calculates the ratio between number of foci in positive and
%negative controls.
%   Images in "PC" and "NC" subfolders are loaded and the number of foci
%   are counted after thresholding with the supplied threshold. The ratio
%   between positive and negative controls is returned and is intended to
%   be used by an optimizer function.

    % Calculate signal from the positive controls
    numPC = length(PCData);
    signalPC = double(0);
    for i = 1:numPC
        result = identifyFoci(PCData{i, 1}, testThresh,...
            config('SE_GAL8_OP'));
        signalPC = signalPC + double(result{1}) / PCData{i, 2};
    end
    signalPCMean = signalPC / numPC;
    
    % Calculate signal from the negative controls
    numNC = length(NCData);
    signalPC = double(0);
    for i = 1:numNC
        result = identifyFoci(NCData{i, 1}, testThresh,...
            config('SE_GAL8_OP'));
        signalNC = signalNC + double(result{1}) / NCData{i, 2};
    end
    signalNCMean = signalNC / numNC;
    
    % Calculate SNR
    SNR = signalPCMean / signalNCMean;
end

