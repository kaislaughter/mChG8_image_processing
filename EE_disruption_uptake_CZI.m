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

% You *must* edit the workingdir and exportdir variables for this code to
% work. You will need to edit the following variables to fit your data:
% galectin8_threshold, nuclear_threshold, These and other parts of code
% that may need to be modified are tagged with "% EditHere" comments which
% can be easily found with Control-F text search

clc, clear, clear all;
display('Choose input directory')
workingdir=[uigetdir(),'/']; % Prompts user for input directory
display('Choose output directory')
exportdir=[uigetdir(),'/']; % Prompts user for output directory
filetype='png';
listing=dir(strcat(workingdir,'*.CZI'));
num_images = length(listing)

%% Parameters
debug=0;

% The Run number is used to track multiple runs of the software, and is used in
% export file names and in the DataCells array. Note: it is a character
% array / string!
run='1'; 

tech_replicates = 3;
plate_columns = 10;

%erosiondisk=strel('disk', 10); % EditHere
%analysisdisk=strel('disk', 20); % EditHere
%tophatdisk=strel('disk',30); % EditHere

% For conveience, a number of sizes of disk shaped structural elements are
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

% These are values of Gal8 and Nuclear stain above background. 
galectin8_threshold=100;
nuclear_threshold=400;
uptake_threshold=80;

% Define structural elements to be used in processing images
se_gal8th = se25; % Gal8 tophat
se_gal8op = se3; % Gal8 open
se_nucop = se25; % Nucleus open
se_uptaketh = se25; % Uptake tophat
se_uptakeop = se3; % Uptake open

%% Analysis

% Initialize cell variable for all data.
outputHeaders = {'Run','Filename','Well #','# Cells','Gal8 sum',...
    'Gal8/cell','# Foci','# Foci/cell','NP signal sum','# NPs',...
    '# NPs/cell', 'Gal8-NP correlation coefficient',...
    'Fraction Gal8 channel overlapping with NPs',...
    'Fraction NP channel overlapping with Gal8'};
outputArray = cell(num_images + 1, length(outputHeaders));
outputArray(1,:) = outputHeaders;

for j=1:num_images % all images
    
    currfile=strcat(workingdir,listing(j,1).name); % defines current file
    fn = listing(j,1).name;
    wellNum = extractBetween(fn, '(', ')');
    wellNum = wellNum{1};
  
    data = bfopen(currfile); % opens image using Bio-Formats
    series1 = data{1,1}; % stores the first image of the stack as series1

    gal8 = series1{1, 1}; % Gal-8 image is the first channel
    nuc = series1{2, 1}; % Nuclei image is the second channel
    uptake = series1{3, 1}; % Cy5 (or DiD) image is the third channel
    
    gal8_bg = mode(gal8,'all');
    nuc_bg = mode(nuc,'all');
    uptake_bg = mode(uptake,'all');
    
    fname = currfile; % sets fname to current file
    well = listing(j,1).name(1:end-4);
    
    % Notably, at this point, you should have your nuclear image living
    % within "nuc" and your gal8 image living within "gal8"
    
    % Identify gal8 foci.
    fociCells = identifyFoci(gal8, se_gal8th, se_gal8op);
    numFoci = fociCells{1};
    gal8ThreshClean = fociCells{2};
    gal8Labels = fociCells{3};
    focimapc=(label2rgb(gal8Labels,'jet','w','shuffle')); % rainbow map of foci
    gal8Circles = xor(imdilate(gal8ThreshClean, se10),...
        imdilate(gal8ThreshClean, se8));
    
    % Identify endocytosed NPs.
    NPCells = identifyFoci(uptake, se_uptaketh, se_uptakeop);
    numNPs = NPCells{1};
    NPThreshClean = NPCells{2};
    NPLabels = NPCells{3};
    npmapc=(label2rgb(NPLabels,'jet','w','shuffle')); % rainbow map of foci
    NPCircles = xor(imdilate(gal8ThreshClean, se10),...
        imdilate(gal8ThreshClean, se8));
    
    % Identify nuclei.
    nucCells = identifyNuclei(nuc, se_nucop);
    numNuclei = nucCells{1};
    nucThreshClean = nucCells{2};
    nucLabels = nucCells{3};
    
    % Calculate correlation between gal8 foci and NPs.
    mandersCells = Manders_ED(gal8ThreshClean, NPThreshClean);
    rP = mandersCells{1};
    rOverlap = mandersCells{2};
    rch1 = mandersCells{3};
    rch2 = mandersCells{4};
    ch1Overlap = mandersCells{5};
    ch2Overlap = mandersCells{6};
    
    % Generate output composite images
    comp=cat(3,NPThreshClean.*100, gal8ThreshClean.*100,...
        nucThreshClean.*50);
    circled=cat(3,gal8Circles.*2^16,gal8.*100-2e4,nuc.*50);
    
    % Integrates Cy5 channel for uptake
    cy5int = sum(uptake(:));

    % measurement
    gal8Sum=sum(gal8(gal8ThreshClean)); %Integrate Gal8 pixel intensities within gal8pos2 mask
    NPSum = sum(uptake(NPThreshClean));
    outputArray(j+1,:) = {run, well, wellNum, numNuclei,...
        gal8Sum, uint64(gal8Sum) / uint64(numNuclei), numFoci,...
        double(numFoci)/double(numNuclei), NPSum, numNPs,...
        double(numNPs)/double(numNuclei), rP, ch1Overlap,...
        ch2Overlap}
      
    %begin exports
    exportbase=strcat(exportdir,well,'_',run,'_');
    %imwrite(nucmap,strcat(exportbase,'nucmap','.png'),'png');
    %imwrite(circled,strcat(exportbase,'composite_circled','.png'),'png');
    imwrite(comp,strcat(exportbase,'composite','.png'),'png');
    
end

% Export an Excel sheet of your data
writetable(cell2table(outputArray),strcat(exportdir,'Output.xlsx')); 
