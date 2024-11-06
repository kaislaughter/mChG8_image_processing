   %% Actin intensity and NP uptake quantification script

% Adapted in 2024 by Kai Slaughter and Eric Donders under the supervision
% of Molly Shoichet in the University of Toronto.

% Original by Kameron V Kilchrist in the Biomedical Engineering Department
% at Vanderbilt University 2015 - 2018 in the course of his PhD studies
% advised by Craig Duvall (https://my.vanderbilt.edu/duvall/)

% University Email: kameron.v.kilchrist@vanderbilt.edu 
% Permanent Email: kameronkilchrist@gmail.com 
% Web: https://kameronkilchrist.com

% This code may be reused at will. 
% FigShare DOI: 10.6084/m9.figshare.7066472
% Web: http://doi.org/10.6084/m9.figshare.7066472
% Licensed under a Creative Commons Attribution 4.0 International License.

% Derivative code should be published at GitHub or FigShare. Derivative
% scientific works should cite Kilchrist, Dimobi, ..., Duvall; "Gal8
% Visualization of Endosome Disruption Predicts Carrier-mediated Biologic
% Drug Intracellular Bioavailability"

%% Actin Quantification and export
% This will measure the total actin intensity and count cell
% nuclei in each frame.

% This code assumes you have 2 (or 3) frame CZI images consisting of
% a nuclear stain (default channel 1, e.g., DAPI or Hoechest), an actin 
% fluorescent stain (default channel 3, e.g., phalloidin-488) and 
% fluorescently labelled nanoparticles (default channel 2, e.g. Cy5 or DiD)

% The adapted script used images captured with a 40x objective and 14 bit
% depth.

clc, clear, close all;
warning('OFF', 'MATLAB:xlswrite:AddSheet'); % Disable new sheet warning.

%% Setup
disp('Choose input directory')
workingdir = [uigetdir(), filesep]; % Prompts user for input directory
disp('Choose output directory')
exportdir = [uigetdir(), filesep]; % Prompts user for output directory

% Load the configuration data from the output directory
% Note: a single Map object called "config" is created in this script.
run(strcat(exportdir, 'Config.m'));

% Fill in defaults for config keys that may be absent from old versions.
optionalKeys = {...
    'EXPORT_NUC_MAP',...
    'EXPORT_NP_ANNOTATIONS',...
    'EXPORT_UNFILTERED_COMPOSITE',...
    'EXPORT_COMPOSITE',...
    'EXPORT_INDIVIDUAL_CHANNELS',...
    };
% List of optional keys not already present in the config Map.
extraKeys = optionalKeys(~isKey(config, optionalKeys));
% Set these missing optional keys to false (the default)
for i = 1:length(extraKeys)
    config(extraKeys{i}) = false;
end


listing = dir(strcat(workingdir, config('IMAGETYPE')));
numImages = length(listing);

se8 = strel('disk', 8);
se10 = strel('disk', 10);

%% Analysis
% Validate configuration values.
% Note, groupTitles must be assigned as its own variable because temporary
% objects can't be indexed.
groupTitles = config('GROUP_TITLES');
if ~config('TEST_MODE')
    if numImages ~= length(groupTitles) * config('TECH_REPLICATES')
        error(['Mismatch between expected '...
            num2str(length(groupTitles) * config('TECH_REPLICATES'))...
            ' and found ' num2str(numImages) ' number of images!']);
    end
end

% Initialize cell variable for all data.
outputHeaders = {...
    'Run #',...
    'Group #',...
    'Well #',...
    'Filename',...
    'Group',...
    '# Cells',...
    'Actin sum',...
    'Actin/cell',...
    'NP signal sum',...
    'NP signal/cell',...
    '# NPs',...
    '# NPs/cell'};
outputArray = cell(numImages, length(outputHeaders));
figureArray = zeros(numImages, 1);

%%
for i = 1:numImages % Iterate over all images.
    %% Image loading and cleanup
    % Select the ith image
    imageTitle = listing(i,1).name(1:end-4);
    currFile = strcat(workingdir, listing(i,1).name);
    fn = listing(i,1).name;
    
    if config('TEST_MODE')
        groupNum = 1;
        wellNum = 1;
    else
        % Extract the well number from the filename.
        wellNum = extractBetween(fn, '(', ')');
        wellNum = str2double(wellNum{1});
        groupNum = config('PLATE_COLUMNS') * fix((wellNum - 1)...
            / (config('PLATE_COLUMNS') * config('TECH_REPLICATES')))...
            + rem((wellNum - 1), config('PLATE_COLUMNS')) + 1;
    end
    
    if i == 1
    % Note that Bio-Formats toolbox must be downloaded and placed in the
        % Get the base run name that all images are based on.
        runName = extractBefore(fn, '(');
    end
    
    % Load the image.
    % working directory or MATLAB path.
    data = bfopen(currFile);
    % The image is opened as a stack so just take the first one.
    series1 = data{1, 1};
    % Split the multi-channel image into its three components.
    nuc = series1{config('NUC_CHANNEL'), 1}; % Nuclei image
    if config('ACTIN_CHANNEL')
        actin = series1{config('ACTIN_CHANNEL'), 1};
        actinTH = imtophat(actin, config('SE_ACTIN_TH'));
    end
    if config('NP_CHANNEL')
        uptake = series1{config('NP_CHANNEL'), 1};
        NPTH = imtophat(uptake, config('SE_NP_TH'));
    end
    
    % Clean up images.
    
    NPCorr = NPTH;
    actinCorr = actinTH;
    
    %% Feature identification
    % Identify nuclei.
    disp('Identifying nuclei...');
    nucCells = identifyNuclei(nuc, config('SE_NUC_OP'));
    numNuclei = nucCells{1};
    nucThresh = nucCells{2};
    nucBinaryClean = nucCells{3};
    nucLabels = nucCells{4};
    
    % Quantify actin intensity
    if config('ACTIN_CHANNEL')
        disp('Quantifying actin...');
        actinClean = imopen(actinCorr, config('SE_ACTIN_OP'));
        actinThresh = actin;
        actinThresh(actinClean < config('ACTIN_THRESHOLD')) = 0;
        actinSum = sum(actinThresh, 'all');
    else
        actinSum = NaN;
        % Make an empty channel for the composite image.
        actinThresh = zeros(size(nucThresh));
    end
   
    % Identify endocytosed NPs.
    if config('NP_CHANNEL')
        disp('Identifying endocytosed NPs...');
        NPResults = identifyFoci(NPCorr, config('NP_THRESHOLD'),...
            config('SE_NP_OP'));
        numNPs = NPResults{1};
        NPThresh = NPResults{2};
        NPBinaryClean = NPResults{3};
        NPLabels = NPResults{4};
        NPSum = sum(NPThresh, 'all');
    else
        numNPs = NaN;
        NPSum = NaN;
        % Make an empty channel for the composite image.
        NPThresh = zeros(size(nucThresh));
    end
    
    %% Save the results
    % Assign the results to a row in the output array
    outputArray(i,:) = {...
        config('RUN'),...
        groupNum,...
        wellNum,...
        imageTitle,...
        groupTitles{groupNum},...
        double(numNuclei),...
        actinSum,...
        uint64(actinSum) / double(numNuclei),...
        NPSum,...
        uint64(NPSum)/double(numNuclei),...
        double(numNPs),...
        double(numNPs)/double(numNuclei)};
    clc;

    disp([outputHeaders;outputArray(1:i,:)]);
      
    % Export the specified images.
    exportbase = strcat(exportdir, imageTitle, '_', config('RUN'), '_');
    if config('EXPORT_NUC_MAP')
        % Generate and save "sanity check" rainbow map for nuclei.
        disp('Exporting nuclei map...');
        nucMap = label2rgb(nucLabels, 'jet', 'w', 'shuffle');
        imwrite(nucMap, strcat(exportbase, 'nucmap', '.png'),...
            config('FILETYPE'));
    end
    if config('EXPORT_INDIVIDUAL_CHANNELS')
        disp('Exporting individual channels...');
        zero_channel = zeros(size(nucThresh));
        if config('NP_CHANNEL')
            img = cat(3, NPThresh.*config('NP_BRIGHTEN'),...
                zero_channel, NPThresh.*config('NP_BRIGHTEN'));
            imwrite(img, strcat(...
                exportbase, 'NPs.', config('FILETYPE')),...
                config('FILETYPE'));
        end
        if config('ACTIN_CHANNEL')
            img = cat(3, zero_channel,...
                actinThresh.*config('ACTIN_BRIGHTEN'), zero_channel);
            imwrite(img, strcat(...
                exportbase, 'Actin.', config('FILETYPE')),...
                config('FILETYPE'));
        end
        img = cat(3, zero_channel, zero_channel,...
            nucThresh.*config('NUC_BRIGHTEN'));
        imwrite(img, strcat(...
            exportbase, 'Nuclei.', config('FILETYPE')), config('FILETYPE'));
    end
    if config('EXPORT_NP_ANNOTATIONS') && config('NP_CHANNEL')
        disp('Exporting NP annotations...');
        NPCircles = xor(imdilate(NPBinaryClean, se10),...
            imdilate(NPBinaryClean, se8));
        NPCircled = cat(3, NPThresh.*config('NP_BRIGHTEN'),...
            NPCircles.*2^16, nucThresh.*config('NUC_BRIGHTEN'));
        imwrite(NPCircled, strcat(...
            exportbase, 'composite_NP_circled', '.png'),...
            config('FILETYPE'));
    end
    if config('EXPORT_UNFILTERED_COMPOSITE')
        % Generate output unfiltered composite images
        disp('Exporting unfiltered composite image...');
            comp = cat(3, uptake.*config('NP_BRIGHTEN'),...
                actin.*config('ACTIN_BRIGHTEN'),...
                nuc.*config('NUC_BRIGHTEN')...
                + uptake.*config('NP_BRIGHTEN'));
        imwrite(comp, strcat(exportbase, 'unfiltered_composite', '.png'),...
            config('FILETYPE'));
    end
    if config('EXPORT_COMPOSITE')
        % Generate output composite images
        disp('Exporting composite image...');
            comp = cat(3, NPThresh.*config('NP_BRIGHTEN'),...
                actinThresh.*config('ACTIN_BRIGHTEN'),...
                nucThresh.*config('NUC_BRIGHTEN')...
                + NPThresh.*config('NP_BRIGHTEN'));
        imwrite(comp, strcat(exportbase, 'composite', '.png'),...
            config('FILETYPE'));
    end
end

%%
% Export an Excel sheet of your data
disp('Saving results...');
exportDirSplit = strsplit(exportdir, filesep);
exportFilename = strcat(...
    exportdir, runName, '_', exportDirSplit{length(exportDirSplit) - 1});
exportExcelFile = strcat(exportFilename, '.xlsx');
if ~config('TEST_MODE') 
    outputArray = sortrows(outputArray, [1, 2, 3]);
end
writetable(cell2table(outputArray, 'VariableNames', outputHeaders),...
    exportExcelFile);
for i = 1:length(outputHeaders)
    exportGroupedData(outputHeaders{i}, outputArray(:, i),...
        groupTitles, exportExcelFile, i + 1);
end
disp('Finished processing!');