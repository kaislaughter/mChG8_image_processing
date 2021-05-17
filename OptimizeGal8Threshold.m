function [gal8ThreshOptimal] = OptimizeGal8Threshold(config, baseDir)
%OptimizeGal8Threshold determines the best threshold to use for the gal8
%channel to achieve maximum discrimination of real disrupted endosomes from
%noise.
%   Detailed explanation goes here

    % Initialize paths and variables.
    listingPC = dir(strcat(baseDir, 'PC/', imageType));
    listingNC = dir(strcat(baseDir, 'NC/', imageType));
    numPC = length(listingPC);
    numNC = length(listingNC);
    
    % Iterate over positive and negative control images and pre-process
    % in preparation for optimization.
    %
    % We load the images into memory and count nuclei beforehand so the
    % optimizer function (which is called many times) doesn't repeat these
    % expensive steps.
    PCData = cell(numPC, 2);
    for i = 1:numPC
        data = bfopen(strcat(workingdir, listingPC(i,1).name));
        series1 = data{1, 1};
        % Store the tophat-transformed gal8 channel.
        gal8 = series1{config('GAL8_CHANNEL'), 1};
        PCData{i, 1} = imtophat(gal8, config('SE_GAL8_TH'));
        % Store the number of nuclei.
        nuc = series1{config('NUC_CHANNEL'), 1};
        result = identifyNuclei(nuc, config('SE_NUC_OP'));
        PCData{i, 2} = result{1};
    end
    NCData = cell(numNC, 2);
    for i = 1:numNC
        data = bfopen(strcat(workingdir, listingNC(i,1).name));
        series1 = data{1, 1};
        % Store the tophat-transformed gal8 channel.
        gal8 = series1{config('GAL8_CHANNEL'), 1};
        NCData{i, 1} = imtophat(gal8, config('SE_GAL8_TH'));
        % Store the number of nuclei.
        nuc = series1{config('NUC_CHANNEL'), 1};
        result = identifyNuclei(nuc, config('SE_NUC_OP'));
        NCData{i, 2} = result{1};
    end
    
    % Run an optimization routine to maximize the signal-to-noise ratio.
    optimizerFunction = @(x) -Gal8SNRHelper(config, PCData, NCData, x);
    % Assume the optimum is somewhere between 0.5x and 2x the supplied
    % value.
    lowBound = round(config('GAL8_THRESHOLD') / 2);
    highBound = config('GAL8_THRESHOLD') * 2;
    options = optimoptions('ga', 'Display', 'iter');
    % Note: the input variable is constrained to integers only. Pixel
    % intensities are stored as integers so there is no sense searching
    % through different fractional thresholds.
    result = ga(optimizerFunction, 1, [], [], [], [], lowBound,...
        highBound, [], 1, options);
    gal8ThreshOptimal = result(1);
    
    % Print the results for debugging.
    if gal8ThreshOptimal == lowBound
        warning('SNR maximized at the lower bound. Try decreasing the default threshold.');
    elseif gal8ThreshOptimal == highBound
        warning('SNR maximized at the upper bound. Try increasing the default threshold.');
    else
        disp(strcat('Optimal threshold identified as ',...
            num2str(gal8ThreshOptimal)));
    end
end

