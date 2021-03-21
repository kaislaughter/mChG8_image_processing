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
clc, clear;

%% Parameters
DEBUG = 0;

% Choose which images will be exported and in what format.
FILETYPE = 'png';
EXPORT_NUC_MAP = true;
EXPORT_GAL8_ANNOTATIONS = true;
EXPORT_NP_ANNOTATIONS = true;
EXPORT_COMPOSITE = true;

% The Run number is used to track multiple runs of the software, and is used in
% export file names and in the DataCells array. Note: it is a character
% array / string!
RUN = '1'; 

% Define technical replicates and plate size. The script assumes technical
% replicates are placed in the same column.
% Note that these values are for annotation purposes only; not configuring
% them will not affect quantification in any way.
TECH_REPLICATES = 3;
PLATE_COLUMNS = 10;

% For convenience, a number of sizes of disk shaped structural elements are
% generated for data exploration.
se1=strel('disk',1);
se2=strel('disk',2);
se3=strel('disk',3);
se4=strel('disk',4);
se5=strel('disk',5);
se6=strel('disk',6);
se7=strel('disk',7);
se8=strel('disk',8);
se9=strel('disk',9);
se10=strel('disk',10);
se25=strel('disk',25);
se100=strel('disk',100);

% Define structural elements to be used in processing images
SE_GAL8_TH = se10; % Gal8 tophat
SE_GAL8_OP = se3; % Gal8 open
SE_NUC_OP = se25; % Nucleus open
SE_NP_TH = se10; % Uptake tophat
SE_NP_OP = se3; % Uptake open

%% Analysis
disp('Choose input directory')
workingdir = [uigetdir(), '/']; % Prompts user for input directory
disp('Choose output directory')
exportdir = [uigetdir(), '/']; % Prompts user for output directory
listing = dir(strcat(workingdir, '*.CZI'));
num_images = length(listing);

% Initialize cell variable for all data.
outputHeaders = {'Run', 'Filename', 'Well #', 'Group #', '# Cells',...
    'Gal8 sum', 'Gal8/cell', '# Foci', '# Foci/cell', 'NP signal sum',...
    '# NPs', '# NPs/cell', 'Gal8-NP correlation coefficient',...
    'Overlap of Gal8 onto NPs', 'Overlap of NPs onto Gal8',...
    'Fraction Gal8 channel overlapping with NPs',...
    'Fraction NP channel overlapping with Gal8'};
outputArray = cell(num_images + 1, length(outputHeaders));
outputArray(1,:) = outputHeaders;

for j = 1:num_images % Iterate over all images.
    % Select the jth image
    title = listing(j,1).name(1:end-4);
    currfile = strcat(workingdir, listing(j,1).name);
    fn = listing(j,1).name;
    % Extract the well number from the filename.
    wellNum = extractBetween(fn, '(', ')');
    wellNum = wellNum{1};
    groupNum = 10 * fix(wellNum / (PLATE_COLUMNS * TECH_REPLICATES))...
        + rem(wellNum, PLATE_COLUMNS);
    
    
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
    
    % Identify gal8 foci.
    fociCells = identifyFoci(gal8, SE_GAL8_TH, SE_GAL8_OP);
    numFoci = fociCells{1};
    gal8Thresh = fociCells{2};
    gal8BinaryClean = fociCells{3};
    gal8Labels = fociCells{4};
    focimapc=(label2rgb(gal8Labels, 'jet', 'w', 'shuffle')); % rainbow map of foci
    gal8Sum = sum(gal8Thresh, 'all');
    
    % Identify endocytosed NPs.
    NPCells = identifyFoci(uptake, SE_NP_TH, SE_NP_OP);
    numNPs = NPCells{1};
    NPThresh = NPCells{2};
    NPBinaryClean = NPCells{3};
    NPLabels = NPCells{4};
    npmapc=(label2rgb(NPLabels, 'jet', 'w', 'shuffle')); % rainbow map of foci
    NPSum = sum(NPThresh, 'all');
    
    % Identify nuclei.
    nucCells = identifyNuclei(nuc, SE_NUC_OP);
    numNuclei = nucCells{1};
    nucThresh = nucCells{2};
    nucBinaryClean = nucCells{3};
    nucLabels = nucCells{4};
    
    % Calculate correlation between gal8 foci and NPs.
    mandersCells = Manders_ED(gal8Thresh, NPThresh);
    rP = mandersCells{1};
    rOverlap = mandersCells{2};
    rch1 = mandersCells{3};
    rch2 = mandersCells{4};
    ch1Overlap = mandersCells{5};
    ch2Overlap = mandersCells{6};
    
    % Save the results in the output array.
    clc;
    outputArray(j+1,:) = {RUN, title, wellNum, groupNum, numNuclei,...
        gal8Sum, uint64(gal8Sum) / uint64(numNuclei), numFoci,...
        double(numFoci)/double(numNuclei), NPSum, numNPs,...
        double(numNPs)/double(numNuclei), rP, rch1,...
        rch2, ch1Overlap, ch2Overlap}
      
    % Export the specified images.
    exportbase = strcat(exportdir, title, '_', RUN, '_');
    if EXPORT_NUC_MAP
        % Generate and save "sanity check" rainbow map for nuclei.
        nucMap = label2rgb(nucLabels, 'jet', 'w', 'shuffle');
        imwrite(nucMap, strcat(exportbase, 'nucmap', '.png'), FILETYPE);
    end
    if EXPORT_GAL8_ANNOTATIONS
        gal8Circles = xor(imdilate(gal8BinaryClean, se10),...
            imdilate(gal8BinaryClean, se8));
        gal8Circled = cat(3, gal8Circles.*2^16, gal8Thresh.*100, nucThresh.*50);
        imwrite(gal8Circled, strcat(...
            exportbase, 'composite_Gal8_circled', '.png'), FILETYPE);
    end
    if EXPORT_NP_ANNOTATIONS
        NPCircles = xor(imdilate(NPBinaryClean, se10),...
            imdilate(NPBinaryClean, se8));
        NPCircled = cat(3, NPThresh.*100, NPCircles.*2^16, nucThresh.*50);
        imwrite(NPCircled, strcat(...
            exportbase, 'composite_NP_circled', '.png'), FILETYPE);
    end
    if EXPORT_COMPOSITE
        % Generate output composite images
        comp = cat(3, NPThresh.*100, gal8Thresh.*100, nucThresh.*50);
        imwrite(comp, strcat(exportbase, 'composite', '.png'), FILETYPE);
    end
end

% Export an Excel sheet of your data
writetable(cell2table(outputArray),strcat(exportdir,'EE_disruption_uptake_CZI_output.xlsx')); 
