function [sensitivityIndex] = Gal8SNRHelper(config, PCData, NCData,...
    testThresh)
%Gal8SNRHelper calculates the ratio between number of foci in positive and
%negative controls.
%   Images in "PC" and "NC" subfolders are loaded and the number of foci
%   are counted after thresholding with the supplied threshold. The ratio
%   between positive and negative controls is returned and is intended to
%   be used by an optimizer function.

    % Calculate signal from the positive controls
    numPC = size(PCData, 1);
    signalsPC = zeros(numPC, 1);
    for i = 1:numPC
        result = identifyFoci(PCData{i, 1}, testThresh,...
            config('SE_GAL8_OP'));
        signalsPC(i) = double(result{1}) / PCData{i, 2};
    end
    signalPCMean = mean(signalsPC);
    signalPCSD = std(signalsPC);
    
    % Calculate signal from the negative controls
    numNC = size(NCData, 1);
    signalNC = double(0);
    signalsNC = zeros(numPC, 1);
    for i = 1:numNC
        result = identifyFoci(NCData{i, 1}, testThresh,...
            config('SE_GAL8_OP'));
        signalNC = signalNC + double(result{1}) / NCData{i, 2};
        signalsNC(i) = double(result{1}) / NCData{i, 2};
    end
    signalNCMean = mean(signalsNC);
    signalNCSD = std(signalsNC);
    
    % Calculate sensitivity index
    sensitivityIndex = 4 * (signalPCMean - signalNCMean)...
        / sqrt(signalPCSD.^2 + signalNCSD.^2);
    % Effectively discard any result with infinite sensitivity.
    if sensitivityIndex == Inf
        sensitivityIndex = -sensitivityIndex;
    end
    fprintf('Sensitivity with threshold %i:\t%f\n', testThresh,...
        sensitivityIndex);
end

