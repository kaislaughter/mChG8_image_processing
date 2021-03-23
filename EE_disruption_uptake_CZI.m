%% Gal8 Recruitment MATLAB Program

% Written by Kameron V Kilchrist in the Biomedical Engineering Department
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
clc, clear, close all;
warning('OFF', 'MATLAB:xlswrite:AddSheet'); % Disable new sheet warning.

%% Parameters
% The Run number is used to track multiple runs of the software, and is used in
% export file names and in the DataCells array. Note: it is a character
% array / string!
RUN = '1';

% Choose which images will be exported and in what format.
FILETYPE = 'png';
EXPORT_NUC_MAP = false;
EXPORT_GAL8_ANNOTATIONS = true;
EXPORT_NP_ANNOTATIONS = true;
EXPORT_OVERLAP_ANNOTATIONS = true;
EXPORT_COMPOSITE = true;
EXPORT_CORRELATION_PLOTS = true;
GAL8_BRIGHTEN = 200;  % Signal multiplier for display only.
NP_BRIGHTEN = 200;  % Signal multiplier for display only.
NUC_BRIGHTEN = 50;  % Signal multiplier for display only.

% Define thresholds and correction factors for image processing
GAL8_THRESHOLD = 100;
NP_THRESHOLD = 80;
GAL8_NP_CROSSTALK = 0.27;
NP_GAL8_CROSSTALK = 0.01;

% Suppress sorting to test individual images
TEST_MODE = false;

% Define technical replicates and plate size. The script assumes technical
% replicates are placed in the same column.
% Note that these values are for annotation purposes only; not configuring
% them will not affect quantification in any way.

if TEST_MODE
    TECH_REPLICATES = 1;
    PLATE_COLUMNS = 1;
    GROUP_TITLES = {'Group 1'};
else
    TECH_REPLICATES = 3;
    PLATE_COLUMNS = 10;
    GROUP_TITLES = {...
    'Group 1',...
    'Group 2',...
    'Group 3',...
    'Group 4',...
    'Group 5',...
    'Group 6',...
    'Group 7',...
    'Group 8',...
    'Group 9',...
    'Group 10',...
    'Group 11',...
    'Group 12',...
    'Group 13',...
    'Group 14',...
    'Group 15',...
    'Group 16',...
    'Group 17',...
    'Group 18',...
    'Group 19',...
    'Group 20'};
end

% For convenience, a number of sizes of disk shaped structural elements are
% generated for data exploration.
se1 = strel('disk',1);
se2 = strel('disk',2);
se3 = strel('disk',3);
se4 = strel('disk',4);
se5 = strel('disk',5);
se6 = strel('disk',6);
se7 = strel('disk',7);
se8 = strel('disk',8);
se9 = strel('disk',9);
se10 = strel('disk',10);
se20 = strel('disk', 20);
se25 = strel('disk',25);
se100 = strel('disk',100);

% Define structural elements to be used in processing images
SE_GAL8_TH = se20; % Gal8 tophat
SE_GAL8_OP = se2; % Gal8 open
SE_NUC_OP = se25; % Nucleus open
SE_NP_TH = se20; % Uptake tophat
SE_NP_OP = se2; % Uptake open

%% Analysis
disp('Choose input directory')
workingdir = [uigetdir(), '/']; % Prompts user for input directory
disp('Choose output directory')
exportdir = [uigetdir(), '/']; % Prompts user for output directory
listing = dir(strcat(workingdir, '*.CZI'));
numImages = length(listing);

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
    gal8 = series1{1, 1}; % Gal-8 image is the first channel.
    nuc = series1{2, 1}; % Nuclei image is the second channel.
    uptake = series1{3, 1}; % Cy5 (or DiD) image is the third channel.
    
    % Clean up images.
    gal8TH = imtophat(gal8, SE_GAL8_TH);
    NPTH = imtophat(uptake, SE_NP_TH);
    gal8Corr = gal8TH - NP_GAL8_CROSSTALK * NPTH;
    NPCorr = NPTH - GAL8_NP_CROSSTALK * gal8TH;
    
    % Identify gal8 foci.
    disp('Identifying Gal8 foci...');
    fociResults = identifyFoci(gal8Corr, GAL8_THRESHOLD, SE_GAL8_OP);
    numFoci = fociResults{1};
    gal8Thresh = fociResults{2};
    gal8BinaryClean = fociResults{3};
    gal8Labels = fociResults{4};
    gal8Sum = sum(gal8Thresh, 'all');
    
    % Identify endocytosed NPs.
    disp('Identifying endocytosed NPs...');
    % Correct for some Gal8 fluorescence bleeding into NP channel.
    NPResults = identifyFoci(NPCorr, NP_THRESHOLD, SE_NP_OP);
    numNPs = NPResults{1};
    NPThresh = NPResults{2};
    NPBinaryClean = NPResults{3};
    NPLabels = NPResults{4};
    NPSum = sum(NPThresh, 'all');
    
    % Identify nuclei.
    disp('Identifying nuclei...');
    nucCells = identifyNuclei(nuc, SE_NUC_OP);
    numNuclei = nucCells{1};
    nucThresh = nucCells{2};
    nucBinaryClean = nucCells{3};
    nucLabels = nucCells{4};
    
    % Calculate correlation between gal8 foci and NPs.
    disp('Calculating colocalization of Gal8 and NPs...');
    mandersCells = Manders_ED(gal8Thresh, NPThresh);
    rP = mandersCells{1};
    rOverlap = mandersCells{2};
    MOC = mandersCells{3};
    ch1Overlap = mandersCells{4};
    ch2Overlap = mandersCells{5};
    
    % Identify number of overlapping regions from binary images
    overlapBinaryClean = gal8BinaryClean & NPBinaryClean;
    overlapMap = watershed(~overlapBinaryClean);
    numOverlap = max(overlapMap(:));
    
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
    if EXPORT_GAL8_ANNOTATIONS
        disp('Exporting Gal8 annotations...');
        gal8Circles = xor(imdilate(gal8BinaryClean, se10),...
            imdilate(gal8BinaryClean, se8));
        gal8Circled = cat(3, gal8Circles.*2^16,...
            gal8Thresh.*GAL8_BRIGHTEN, nucThresh.*NUC_BRIGHTEN);
        imwrite(gal8Circled, strcat(...
            exportbase, 'composite_Gal8_circled', '.png'), FILETYPE);
    end
    if EXPORT_NP_ANNOTATIONS
        disp('Exporting NP annotations...');
        NPCircles = xor(imdilate(NPBinaryClean, se10),...
            imdilate(NPBinaryClean, se8));
        NPCircled = cat(3, NPThresh.*NP_BRIGHTEN, NPCircles.*2^16,...
            nucThresh.*NUC_BRIGHTEN);
        imwrite(NPCircled, strcat(...
            exportbase, 'composite_NP_circled', '.png'), FILETYPE);
    end
    if EXPORT_OVERLAP_ANNOTATIONS
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
    if EXPORT_CORRELATION_PLOTS
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
