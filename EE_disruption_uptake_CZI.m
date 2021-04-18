%% Gal8 recruitment and NP uptake quantification script

% Adapted in 2021 by Eric Donders and Kai Slaughter under the supervision
% of Molly Shoichet in the University of Toront.

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

% This code assumes you have 2 (or 3) frame CZI images consisting of a
% Gal8 fluorescent stain in channel 1 (e.g., Gal8-mCherry)and a nuclear
% stain in channel 2 (e.g., DAPI or Hoechest). MATLAB plays more nicely
% with 16 bit data. Original script used images that were 2048 px squared
% with a 20x objective, so if you use larger or smaller images or higher
% mag objectives, you will likely need to edit this.

% The adapted script used images captured with a 40x objective and 14 bit
% depth.

% Usage instructions:
% 1. Install the Bio-Formats MATLAB toolbox by following the instructions
%    on their site: https://downloads.openmicroscopy.org/bio-formats/5.3.4/
% 2. Put the .CZI images that you want to process in their own folder
% 3. Copy Config.m into the same folder and edit the values according to
%    your needs.
% 4. Make an empty folder for the output data table and optional images
% 5. Run this script and wait for it to complete.
clc, clear, close all;
warning('OFF', 'MATLAB:xlswrite:AddSheet'); % Disable new sheet warning.

%% Setup
disp('Choose input directory')
workingdir = [uigetdir(), '/']; % Prompts user for input directory
disp('Choose output directory')
exportdir = [uigetdir(), '/']; % Prompts user for output directory
listing = dir(strcat(workingdir, '*.CZI'));
numImages = length(listing);

% Load the configuration data from the input directory.
run(strcat(workingdir, 'Config.m'));

%% Analysis
% Validate configuration values.
if ~TEST_MODE
    if numImages ~= length(GROUP_TITLES) * TECH_REPLICATES
        error(['Mismatch between expected '...
            num2str(length(GROUP_TITLES) * TECH_REPLICATES)...
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
    '# NPs',...
    '# NPs/cell',...
    'Gal8-NP correlation coefficient',...
    'Manders overlap coefficient',...
    'Fraction Gal8 signal overlapping with NPs',...
    'Fraction NP signal overlapping with Gal8',...
    '# Number of overlaps betweeen Gal8 and NPs',...
    'Fraction Gal8 foci overlapping with NPs',...
    'Fraction NPs overlapping with Gal8'};
outputArray = cell(numImages, length(outputHeaders));
figureArray = zeros(numImages, 1);

for i = 1:numImages % Iterate over all images.
    % Select the jth image
    imageTitle = listing(i,1).name(1:end-4);
    currfile = strcat(workingdir, listing(i,1).name);
    fn = listing(i,1).name;
    
    if TEST_MODE
        groupNum = 1;
        wellNum = 1;
    else
        % Extract the well number from the filename.
        wellNum = extractBetween(fn, '(', ')');
        wellNum = str2double(wellNum{1});
        groupNum = PLATE_COLUMNS * fix((wellNum - 1)...
            / (PLATE_COLUMNS * TECH_REPLICATES))...
            + rem((wellNum - 1), PLATE_COLUMNS) + 1;
    end
    
    if i == 1
        % Get the base run name that all images are based on.
        runName = extractBefore(fn, '(');
    end
    
    % Load the image.
    % Note that Bio-Formats toolbox must be downloaded and placed in the
    % working directory or MATLAB path.
    data = bfopen(currfile);
    % The image is opened as a stack so just take the first one.
    series1 = data{1, 1};
    % Split the multi-channel image into its three components.
    nuc = series1{NUC_CHANNEL, 1}; % Nuclei image
    if GAL8_CHANNEL
        gal8 = series1{GAL8_CHANNEL, 1};
        gal8TH = imtophat(gal8, SE_GAL8_TH);
    end
    if NP_CHANNEL
        uptake = series1{NP_CHANNEL, 1};
        NPTH = imtophat(uptake, SE_NP_TH);
    end
    
    % Clean up images.
    if GAL8_CHANNEL && NP_CHANNEL
        gal8Corr = gal8TH - NP_GAL8_CROSSTALK * NPTH;
        NPCorr = NPTH - GAL8_NP_CROSSTALK * gal8TH;
    elseif ~GAL8_CHANNEL
        NPCorr = NPTH;
    else
        gal8Corr = gal8TH;
    end
    
    % Identify nuclei.
    disp('Identifying nuclei...');
    nucCells = identifyNuclei(nuc, SE_NUC_OP);
    numNuclei = nucCells{1};
    nucThresh = nucCells{2};
    nucBinaryClean = nucCells{3};
    nucLabels = nucCells{4};
    
    % Identify gal8 foci.
    if GAL8_CHANNEL
        disp('Identifying Gal8 foci...');
        fociResults = identifyFoci(gal8Corr, GAL8_THRESHOLD, SE_GAL8_OP);
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
    if NP_CHANNEL
        disp('Identifying endocytosed NPs...');
        % Correct for some Gal8 fluorescence bleeding into NP channel.
        NPResults = identifyFoci(NPCorr, NP_THRESHOLD, SE_NP_OP);
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
    
    % Calculate correlation between gal8 foci and NPs.
    if GAL8_CHANNEL && NP_CHANNEL
        disp('Calculating colocalization of Gal8 and NPs...');
        mandersCells = Manders_ED(gal8Thresh, NPThresh);
        % Identify number of overlapping regions from binary images
        overlapBinaryClean = gal8BinaryClean & NPBinaryClean;
        overlapMap = watershed(~overlapBinaryClean);
        numOverlap = max(overlapMap(:));
    else
        % We can't calculate overlap without both channels.
        mandersCells = {NaN, NaN, NaN, NaN, NaN};
        numOverlap = NaN;
    end
    rP = mandersCells{1};
    rOverlap = mandersCells{2};
    MOC = mandersCells{3};
    ch1Overlap = mandersCells{4};
    ch2Overlap = mandersCells{5};
    
    % Save the results in the output array.
    outputArray(i,:) = {...
        RUN,...
        groupNum,...
        wellNum,...
        imageTitle,...
        GROUP_TITLES{groupNum},...
        double(numNuclei),...
        gal8Sum,...
        uint64(gal8Sum) / uint64(numNuclei),...
        double(numFoci),...
        double(numFoci)/double(numNuclei),...
        double(NPSum),...
        double(numNPs),...
        double(numNPs)/double(numNuclei),...
        rP,...
        MOC,...
        ch1Overlap,...
        ch2Overlap,...
        numOverlap,...
        double(numOverlap)/double(numFoci),...
        double(numOverlap)/double(numNPs)};
    clc;
    disp(outputArray(1:i,:));
      
    % Export the specified images.
    exportbase = strcat(exportdir, imageTitle, '_', RUN, '_');
    if EXPORT_NUC_MAP
        % Generate and save "sanity check" rainbow map for nuclei.
        disp('Exporting nuclei map...');
        nucMap = label2rgb(nucLabels, 'jet', 'w', 'shuffle');
        imwrite(nucMap, strcat(exportbase, 'nucmap', '.png'), FILETYPE);
    end
    if EXPORT_GAL8_ANNOTATIONS && GAL8_CHANNEL
        disp('Exporting Gal8 annotations...');
        gal8Circles = xor(imdilate(gal8BinaryClean, se10),...
            imdilate(gal8BinaryClean, se8));
        gal8Circled = cat(3, gal8Circles.*2^16,...
            gal8Thresh.*GAL8_BRIGHTEN, nucThresh.*NUC_BRIGHTEN);
        imwrite(gal8Circled, strcat(...
            exportbase, 'composite_Gal8_circled', '.png'), FILETYPE);
    end
    if EXPORT_NP_ANNOTATIONS && NP_CHANNEL
        disp('Exporting NP annotations...');
        NPCircles = xor(imdilate(NPBinaryClean, se10),...
            imdilate(NPBinaryClean, se8));
        NPCircled = cat(3, NPThresh.*NP_BRIGHTEN, NPCircles.*2^16,...
            nucThresh.*NUC_BRIGHTEN);
        imwrite(NPCircled, strcat(...
            exportbase, 'composite_NP_circled', '.png'), FILETYPE);
    end
    if EXPORT_OVERLAP_ANNOTATIONS && GAL8_CHANNEL && NP_CHANNEL
        disp('Exporting colocalization annotations...');
        overlapCircles = xor(imdilate(overlapBinaryClean, se10),...
            imdilate(overlapBinaryClean, se8));
        overlapCircled = cat(3,...
            NPThresh.*NP_BRIGHTEN + uint16(overlapCircles.*2^16),...
            gal8Thresh.*GAL8_BRIGHTEN + uint16(overlapCircles.*2^16),...
            nucThresh.*NUC_BRIGHTEN + uint16(overlapCircles.*2^16));
        imwrite(overlapCircled, strcat(...
            exportbase, 'composite_overlap_circled', '.png'), FILETYPE);
    end
    if EXPORT_COMPOSITE
        % Generate output composite images
        disp('Exporting composite image...');
        comp = cat(3, NPThresh.*NP_BRIGHTEN, gal8Thresh.*GAL8_BRIGHTEN,...
            nucThresh.*NUC_BRIGHTEN);
        imwrite(comp, strcat(exportbase, 'composite', '.png'), FILETYPE);
    end
    if EXPORT_CORRELATION_PLOTS && GAL8_CHANNEL && NP_CHANNEL
        % Plot the intensity of Gal8 versus the intensity of NPs.
        disp('Exporting a correlation plot...');
        fig = figure(i);
        scatter(NPCorr(:), gal8Corr(:), '.');
        xlim([0 2000]);
        ylim([0 2000]);
        xline(NP_THRESHOLD);
        yline(GAL8_THRESHOLD);
        xlabel('NP pixel intensity');
        ylabel('Gal8 pixel intensity');
        titleStr = [GROUP_TITLES{groupNum} ' (well ' num2str(wellNum) ')'];
        title(titleStr);
        saveas(fig, strcat(exportdir, titleStr), FILETYPE);
        close;
    end
end

% Export an Excel sheet of your data
disp('Saving results...');
exportFilename = strcat(exportdir, runName, '_Output');
exportExcelFile = strcat(exportFilename, '.xlsx');
if ~TEST_MODE 
    outputArray = sortrows(outputArray, [1, 2, 3]);
end
writetable(cell2table(outputArray, 'VariableNames', outputHeaders),...
    exportExcelFile);
for i = 1:length(outputHeaders)
    exportGroupedData(outputHeaders{i}, outputArray(:, i),...
        GROUP_TITLES, exportExcelFile, i + 1);
end
disp('Finished processing!');