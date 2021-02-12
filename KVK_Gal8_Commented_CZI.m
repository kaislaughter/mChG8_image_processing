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

%% Parameters
debug=0;
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

% The Run number is used to track multiple runs of the software, and is used in
% export file names and in the DataCells array. Note: it is a character
% array / string!
run='1'; 

% These are values of Gal8 and Nuclear stain above background. 
galectin8_threshold=100;
nuclear_threshold=400;


%% Analysis
DataCells = [{'Run'},{'Well'},{'# Cells'},{'Gal8 sum'},{'Gal8/cell'},...
    {'# Foci'},{'# Foci/cell'}]; % Initializes cell variable for all data

    % The next few lines are specific to CZI images. Edit from here
    % if you have alternate arrangements. 

for j=1:length(listing) % all images
    
    currfile=strcat(workingdir,listing(j,1).name); % defines current file
  
    data = bfopen(currfile); % opens image using Bio-Formats
    series1 = data{1,1}; % stores the first image of the stack as series1

    gal8 = series1{1, 1}; % Gal-8 image is the first channel
    nuc = series1{2, 1}; % Nuclei image is the second channel

    fname = currfile; % sets fname to current file
    well=listing(j,1).name(1:end-4)
    
    % Notably, at this point, you should have your nuclear image living
    % within "nuc" and your gal8 image living within "gal8"
    
    gal8th=imtophat(gal8,se25); % Clean image with tophat filter for thresholding 
    gal8pos1=gal8th>galectin8_threshold; % threshold image
    gal8pos2=imopen(gal8pos1,se3); % open thresholded image
    circlelayer=xor(imdilate(gal8pos2,se10),imdilate(gal8pos2,se8));
    % Create circle layer for circled ouptut images
    
    % Generate output composite images
    comp=cat(3,circlelayer.*0,gal8.*100-2e4,nuc.*50);
    circled=cat(3,circlelayer.*2^16,gal8.*100-2e4,nuc.*50);
    
    % Count foci
    basinmap2=watershed(~gal8pos2); % watershed to count foci
    %fociarea=imdilate(gal8pos2,se1);
    fociarea=gal8pos2;
    focimap=basinmap2; 
    focimap(~fociarea)=0;
    focimapc=(label2rgb(focimap,'jet','w','shuffle')); % rainbow map of foci
    numfoci=max(focimap(:)); % stores count of foci
    
    % Code below counts cell number
    nthr=nuc>nuclear_threshold; % thresholding
    nuc1=(imopen(nthr,se25)); % remove small features
    disttrans=-bwdist(~nuc1); % distance transform
    mask=imextendedmin(disttrans,2); % removes noise from distance transform
    disttrans2=imimposemin(disttrans,mask);
    basinmap=watershed(disttrans2); % watershed to count nuclei
    cellarea=imdilate(nuc1, se3); % dilates image for output 
    cellmap=basinmap;
    cellmap(~cellarea)=0;
    nucmap=(label2rgb(cellmap,'jet','w','shuffle')); %Generates "sanity check" rainbow map
    
    % measurement
    numcell=max(cellmap(:)); % stores nuclei count as numcell
    galsum=sum(gal8(gal8pos2)); %Integrate Gal8 pixel intensities within gal8pos2 mask
    DataCells=[DataCells;{run},{well},{numcell},{galsum},{uint64(galsum)...
        /uint64(numcell)},{numfoci},{double(numfoci)/double(numcell)}]
      
    %begin exports
    exportbase=strcat(exportdir,well,'_',run,'_');
    imwrite(nucmap,strcat(exportbase,'nucmap','.png'),'png');
    imwrite(circled,strcat(exportbase,'composite_circled','.png'),'png');
    imwrite(comp,strcat(exportbase,'composite','.png'),'png');
    
end
writetable(cell2table(DataCells),strcat('Output.xls')); % Exports an csv sheet of your data
