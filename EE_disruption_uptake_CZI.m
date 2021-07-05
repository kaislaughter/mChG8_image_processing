%% Gal8 recruitment and NP uptake quantification script

% Adapted in 2021 by Eric Donders and Kai Slaughter under the supervision
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

%% HTS Gal8 Quantification and export
% This will measure the total Gal8 recruited to puncta and count cell
% nuclei in each frame.

% This code assumes you have 2 (or 3) frame CZI images consisting of
% a nuclear stain (default channel 2, e.g., DAPI or Hoechest), a Gal8 
% fluorescent stain (default channel 1, e.g., Gal8-mCherry) and 
% fluorescently labelled nanoparticles (default channel 3, e.g. Cy5 or DiD)

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

listing = dir(strcat(workingdir, config('IMAGETYPE')));
numImages = length(listing);

se8 = strel('disk', 8);
se10 = strel('disk', 10);

%% Gal8 threshold optimization
% Note: the Global Optimization Toolbox must be installed to use this
% feature.
if config('OPTIMIZE_GAL8_THRESHOLD')
    if ~contains(struct2array(ver), 'Global Optimization Toolbox')
        error(['The Global Optimization Toolbox is required but not '...
            'installed.']);
    end
    optimalGal8Thresh = OptimizeGal8Threshold(config, workingdir);
    fprintf('Gal8 threshold optimized from %i to %i\n',...
        config('GAL8_THRESHOLD'), optimalGal8Thresh);
    config('GAL8_THRESHOLD') = optimalGal8Thresh;
end

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
    'Gal8 sum',...
    'Gal8/cell',...
    '# Foci',...
    '# Foci/cell',...
    'NP signal sum',...
    'NP signal/cell',...
    '# NPs',...
    '# NPs/cell',...
    'PL signal sum',...
    'PL signal/cell',...
    '# PL foci',...
    '# PL foci/cell',...
    'Gal8-NP correlation coefficient',...
    'Gal8-NP Manders overlap coefficient',...
    'Fraction Gal8 signal overlapping with NPs',...
    'Fraction NP signal overlapping with Gal8',...
    '# Number of overlaps betweeen Gal8 and NPs',...
    'Fraction Gal8 foci overlapping with NPs',...
    'Fraction NPs overlapping with Gal8',...
    'Gal8-PL correlation coefficient',...
    'Gal8-PL Manders overlap coefficient',...
    'Fraction Gal8 signal overlapping with PL',...
    'Fraction PL signal overlapping with Gal8',...
    '# Number of overlaps betweeen Gal8 and PL',...
    'Fraction Gal8 foci overlapping with PL',...
    'Fraction PL overlapping with Gal8'...
    'NP-PL correlation coefficient',...
    'NP-PL Manders overlap coefficient',...
    'Fraction NPs signal overlapping with PL',...
    'Fraction PL signal overlapping with NPs',...
    '# Number of overlaps betweeen NPs and PL',...
    'Fraction NPs foci overlapping with PL',...
    'Fraction PL overlapping with NPs'};
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
    if config('GAL8_CHANNEL')
        gal8 = series1{config('GAL8_CHANNEL'), 1};
        gal8TH = imtophat(gal8, config('SE_GAL8_TH'));
    end
    if config('NP_CHANNEL')
        uptake = series1{config('NP_CHANNEL'), 1};
        NPTH = imtophat(uptake, config('SE_NP_TH'));
    end
    if config('PL_CHANNEL')
        PL = series1{config('PL_CHANNEL'), 1};
        PLTH = imtophat(PL, config('SE_PL_TH'));
    end
    
    % Clean up images.
    if config('GAL8_CHANNEL') && config('NP_CHANNEL')
        gal8Corr = gal8TH - config('NP_GAL8_CROSSTALK') * NPTH;
        NPCorr = NPTH - config('GAL8_NP_CROSSTALK')...
            * gal8TH;     
    elseif ~config('GAL8_CHANNEL')
        NPCorr = NPTH;
    else
        gal8Corr = gal8TH;
    end
    
    %% Feature identification
    % Identify nuclei.
    disp('Identifying nuclei...');
    nucCells = identifyNuclei(nuc, config('SE_NUC_OP'));
    numNuclei = nucCells{1};
    nucThresh = nucCells{2};
    nucBinaryClean = nucCells{3};
    nucLabels = nucCells{4};
    
    % Identify gal8 foci.
    if config('GAL8_CHANNEL')
        disp('Identifying Gal8 foci...');
        fociResults = identifyFoci(gal8Corr, config('GAL8_THRESHOLD'),...
            config('SE_GAL8_OP'));
        numFoci = fociResults{1};
        gal8Thresh = fociResults{2};
        gal8BinaryClean = fociResults{3};
        gal8Labels = fociResults{4};
        gal8Sum = sum(gal8Thresh, 'all');
    else
        numFoci = NaN;
        gal8Sum = NaN;
        % Make an empty channel for the composite image.
        gal8Thresh = zeros(size(nucThresh));
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
    
    % Identify phospholipidosis vesicles
    if config('PL_CHANNEL')
        disp('Identifying phospholipidosis vesicles...');
        PLResults = identifyFoci(PLTH, config('PL_THRESHOLD'),...
            config('SE_PL_OP'));
        numPLs = PLResults{1};
        PLThresh = PLResults{2};
        PLBinaryClean = PLResults{3};
        PLLabels = PLResults{4};
        PLSum = sum(PLThresh, 'all');
    end        
    
    %% Correlation calculations
    % Calculate correlation between gal8 foci and NPs.
    if config('GAL8_CHANNEL') && config('NP_CHANNEL')
        disp('Calculating colocalization of Gal8 and NPs...');
        gal8NPmandersCells = Manders_ED(gal8Thresh, NPThresh);
        % Identify number of overlapping regions from binary images
        overlapBinaryClean = gal8BinaryClean & NPBinaryClean;
        overlapMap = watershed(~overlapBinaryClean);
        gal8NPnumOverlap = max(overlapMap(:));
    else
        % We can't calculate overlap without both channels.
        gal8NPmandersCells = {NaN, NaN, NaN, NaN, NaN};
        gal8NPnumOverlap = NaN;
    end
    gal8NPrP = gal8NPmandersCells{1};
    gal8NPrOverlap = gal8NPmandersCells{2};
    gal8NPMOC = gal8NPmandersCells{3};
    gal8NPch1Overlap = gal8NPmandersCells{4};
    gal8NPch2Overlap = gal8NPmandersCells{5};
    
    % Calculate correlation between gal8 and phospholipidosis.
    if config('GAL8_CHANNEL') && config('PL_CHANNEL')
        disp('Calculating colocalization of Gal8 and phospholipidosis...');
        mandersCells = Manders_ED(gal8Thresh, PLThresh);
        % Identify number of overlapping regions from binary images
        overlapBinaryClean = gal8BinaryClean & PLBinaryClean;
        overlapMap = watershed(~overlapBinaryClean);
        gal8PLnumOverlap = max(overlapMap(:));
    else
        % We can't calculate overlap without both channels.
        gal8NPmandersCells = {NaN, NaN, NaN, NaN, NaN};
        gal8PLnumOverlap = NaN;
    end
    gal8PLrP = mandersCells{1};
    gal8PLrOverlap = mandersCells{2};
    gal8PLMOC = mandersCells{3};
    gal8PLch1Overlap = mandersCells{4};
    gal8PLch2Overlap = mandersCells{5};
    
    % Calculate correlation between NPs and phospholipidosis.
    if config('NP_CHANNEL') && config('PL_CHANNEL')
        disp('Calculating colocalization of NPs and phospholipidosis...');
        mandersCells = Manders_ED(NPThresh, PLThresh);
        % Identify number of overlapping regions from binary images
        overlapBinaryClean = NPBinaryClean & PLBinaryClean;
        overlapMap = watershed(~overlapBinaryClean);
        NPPLnumOverlap = max(overlapMap(:));
    else
        % We can't calculate overlap without both channels.
        gal8NPmandersCells = {NaN, NaN, NaN, NaN, NaN};
        NPPLnumOverlap = NaN;
    end
    NPPLrP = mandersCells{1};
    NPPLrOverlap = mandersCells{2};
    NPPLMOC = mandersCells{3};
    NPPLch1Overlap = mandersCells{4};
    NPPLch2Overlap = mandersCells{5};
    
    %% Save the results
    % Assign the results to a row in the output array
    outputArray(i,:) = {...
        config('RUN'),...
        groupNum,...
        wellNum,...
        imageTitle,...
        groupTitles{groupNum},...
        double(numNuclei),...
        gal8Sum,...
        uint64(gal8Sum) / double(numNuclei),...
        double(numFoci),...
        double(numFoci)/double(numNuclei),...
        NPSum,...
        uint64(NPSum)/double(numNuclei),...
        double(numNPs),...
        double(numNPs)/double(numNuclei),...
        PLSum,...
        uint64(PLSum) / double(numNuclei),...
        double(numPLs),...
        double(numPLs) / double(numNuclei),...
        gal8NPrP,...
        gal8NPMOC,...
        gal8NPch1Overlap,...
        gal8NPch2Overlap,...
        gal8NPnumOverlap,...
        double(gal8NPnumOverlap)/double(numFoci),...
        double(gal8NPnumOverlap)/double(numNPs),...
        gal8PLrP,...
        gal8PLMOC,...
        gal8PLch1Overlap,...
        gal8PLch2Overlap,...
        gal8PLnumOverlap,...
        double(gal8PLnumOverlap)/double(numFoci),...
        double(gal8PLnumOverlap)/double(numNPs),...
        NPPLrP,...
        NPPLMOC,...
        NPPLch1Overlap,...
        NPPLch2Overlap,...
        NPPLnumOverlap,...
        double(NPPLnumOverlap)/double(numFoci),...
        double(NPPLnumOverlap)/double(numNPs)};
    clc;
    disp(outputArray(1:i,:));
      
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
        if config('GAL8_CHANNEL')
            img = cat(3, zero_channel,...
                gal8Thresh.*config('GAL8_BRIGHTEN'), zero_channel);
            imwrite(img, strcat(...
                exportbase, 'Gal8_foci.', config('FILETYPE')),...
                config('FILETYPE'));
        end
        if config('PL_CHANNEL')
            img = cat(3, zero_channel, PLThresh.*config('PL_BRIGHTEN'),...
                PLThresh.*config('PL_BRIGHTEN'));
            imwrite(img, strcat(...
                exportbase, 'PL_foci.', config('FILETYPE')),...
                config('FILETYPE'));
        end
        img = cat(3, zero_channel, zero_channel,...
            nucThresh.*config('NUC_BRIGHTEN'));
        imwrite(img, strcat(...
            exportbase, 'Nuclei.', config('FILETYPE')), config('FILETYPE'));
    end
    if config('EXPORT_GAL8_ANNOTATIONS') && config('GAL8_CHANNEL')
        disp('Exporting Gal8 annotations...');
        gal8Circles = xor(imdilate(gal8BinaryClean, se10),...
            imdilate(gal8BinaryClean, se8));
        gal8Circled = cat(3, gal8Circles.*2^16,...
            gal8Thresh.*config('GAL8_BRIGHTEN'),...
            nucThresh.*config('NUC_BRIGHTEN'));
        imwrite(gal8Circled, strcat(...
            exportbase, 'composite_Gal8_circled', '.png'),...
            config('FILETYPE'));
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
    if config('EXPORT_PL_ANNOTATIONS') && config('PL_CHANNEL')
        disp('Exporting PL annotations...');
        PLCircles = xor(imdilate(PLBinaryClean, se10),...
            imdilate(PLBinaryClean, se8));
        PLCircled = cat(3, PLThresh.*config('PL_BRIGHTEN'),...
            PLCircles.*2^16, nucThresh.*config('NUC_BRIGHTEN'));
        imwrite(PLCircled, strcat(...
            exportbase, 'composite_PL_circled', '.png'),...
            config('FILETYPE'));
    end
    if config('EXPORT_OVERLAP_ANNOTATIONS') && config('GAL8_CHANNEL')...
            && config('NP_CHANNEL')
        disp('Exporting colocalization annotations...');
        overlapCircles = xor(imdilate(overlapBinaryClean, se10),...
            imdilate(overlapBinaryClean, se8));
        overlapCircled = cat(3,...
            NPThresh.*config('NP_BRIGHTEN')...
            + uint16(overlapCircles.*2^16),...
            gal8Thresh.*config('GAL8_BRIGHTEN')...
            + uint16(overlapCircles.*2^16),...
            nucThresh.*config('NUC_BRIGHTEN')...
            + NPThresh.*config('NP_BRIGHTEN')...
            + uint16(overlapCircles.*2^16));
        imwrite(overlapCircled, strcat(...
            exportbase, 'composite_overlap_circled', '.png'),...
            config('FILETYPE'));
    end
    if config('EXPORT_COMPOSITE')
        % Generate output composite images
        disp('Exporting composite image...');
        comp = cat(3, NPThresh.*config('NP_BRIGHTEN'),...
            gal8Thresh.*config('GAL8_BRIGHTEN'),...
            nucThresh.*config('NUC_BRIGHTEN')...
            + NPThresh.*config('NP_BRIGHTEN'));
        imwrite(comp, strcat(exportbase, 'composite', '.png'),...
            config('FILETYPE'));
    end
    if config('EXPORT_CORRELATION_PLOTS') && config('GAL8_CHANNEL')...
            && config('NP_CHANNEL')
        % Plot the intensity of Gal8 versus the intensity of NPs.
        disp('Exporting a correlation plot...');
        fig = figure(i);
        scatter(NPCorr(:), gal8Corr(:), '.');
        xlim([0 2000]);
        ylim([0 2000]);
        xline(config('NP_THRESHOLD'));
        yline(config('GAL8_THRESHOLD'));
        xlabel('NP pixel intensity');
        ylabel('Gal8 pixel intensity');
        titleStr = [groupTitles{groupNum} ' (well ' num2str(wellNum) ')'];
        title(titleStr);
        saveas(fig, strcat(exportdir, titleStr), config('FILETYPE'));
        close;
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