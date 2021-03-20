
%% HTS Lysotracker quantification and export
% Written by Kameron V Kilchrist in the Biomedical Engineering Department
% at Vanderbilt University 2015 - 2018 in the course of his PhD studies
% advised by Craig Duvall (https://my.vanderbilt.edu/duvall/)

% University Email: kameron.v.kilchrist@vanderbilt.edu 
% Permanent Email: kameronkilchrist@gmail.com 
% Web: https://kameronkilchrist.com

% This code may be reused at will. 
% FigShare DOI: 10.6084/m9.figshare.7066490
% Web: http://doi.org/10.6084/m9.figshare.7066490
% Licensed under a Creative Commons Attribution 4.0 International License.

% Derivative code should be published at GitHub or FigShare. Derivative
% scientific works should cite Kilchrist, Dimobi, ..., Duvall; "Gal8
% Visualization of Endosome Disruption Predicts Carrier-mediated Biologic
% Drug Intracellular Bioavailability"

% This will calculate the Manders coefficients in each
% frame.

% This code assumes you have exported 2 frame images consisting of a
% drug stain in channel 1 (e.g., Alexa 488) and a Lysotracker
% fluorescent stain in channel 2 (e.g., Lysotracker Red). Notably, we used two
% channel Nikon ND2 files as our raw data, which are 12 bit data; however,
% MATLAB plays more nicely with 16 bit data so we scale this with ".*2^4"
% Our images were 2048 px squared with a 20x objective, so if you use
% larger or smaller images or higher mag objectives, you will likely need
% to edit this.

% Images were exported to 2 page TIF images using Nikon NIS Elements. If
% you use alternate file formats or use separate export files for each
% channel, please edit the code as appropriate.

% You *must* edit the workingdir and exportdir variables for this code to
% work. You will need to edit the following variables to fit your data:
% galectin8_threshold, nuclear_threshold, These and other parts of code
% that may need to be modified are tagged with "% EditHere" comments which
% can be easily found with Control-F text search

%% Initialize
clc, clear, clear all;
workingdir='C:\ [...] \ExportsND2\'; % With trailing slash! Where images to be analyzed live
exportdir='C:\ [...] \ExportsPNG\'; % With trailing slash! Where exported images go
filetype='png';
listing=dir(strcat(workingdir,'*.TIF'));


%% Parameters - set these based on controls. 
a488_threshold=20000; % EditHere
lysored_threshold=20000; % EditHere
%% Analysis
DataCells=[{'well'},{'sumred'},{'sumgreen'},{'Pearson Spatial Correlation'},{'m1 Red Over Green'},{'m2 Green Over Red'}]; %Clear the matrix for our data...

for j=1:length(listing) % all images
    
    currfile=strcat(workingdir,listing(j,1).name);
    
    fname = currfile;
    info = imfinfo(fname);
    num_images = numel(info);
    well=currfile(end-6:end-4)
    
    a488 = imread(fname, 1, 'Info', info);
    lysored = imread(fname, 2, 'Info', info);
    
    
    comp=cat(3,lysored,a488,lysored);
    
    
    red=lysored-lysored_threshold;
    green=a488-a488_threshold;
    
    rP=corr2(red,green);
    
    products=(double(red).*double(green));
    redsq=double(red).^2;
    greensq=double(green).^2;
    
    rOverlap=sum(products(:))/(sum(redsq(:))*sum(greensq(:)));
    
    rRed=sum(products(:))/sum(redsq(:)); %Overlap of product onto red
    rGreen=sum(products(:))/sum(greensq(:)); %Overlap of green onto red
    
    
    sumgreen=sum(green(:));
    sumred=sum(red(:));
    
    red_in_green_pos=red(green>0);
    green_in_red_pos=green(red>0);
    
    m1=sum(red_in_green_pos(:))/sumred;
    m2=sum(green_in_red_pos(:))/sumgreen;
    
    
    DataCells=[DataCells;{well},{sumred},{sumgreen},{rP},{m1},{m2}];
    
    
    % begin exports
    exportbase=strcat(exportdir,well,'_PurpGrn');
    imwrite(comp,strcat(exportbase,'.png'),'png');
    
end

xlswrite(strcat(exportdir,'DataCells.xls'),DataCells)

%% _________________________________________________________________
