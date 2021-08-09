%% Script for preparing composites for publication

%   Input: exported composite images from EE_disruption_uptake_CZI.m
%   (or similar).
%   IMPORTANT: the image format must conform to the following standard:
%    - Nuclei: blue
%    - Gal8: green
%    - Colloids/NPs: magenta (equal red and blue)

%   Output: individual channels plus annotated colocalization maps, in both
%   the original size and zoomed-in sections, with scale bars.

clc, clear, close all;

%% Config
%   Fill in this section before each run.

pixelSize = 0.1;  % uM PLACEHOLDER
zoomedWidth = 100;  % uM
zoomedHeight = 100;  % uM
scaleBar = 25;  % uM

%% Setup
disp('Choose input directory')
workingdir = [uigetdir(), filesep]; % Prompts user for input directory
disp('Choose output directory')
exportdir = [uigetdir(), filesep]; % Prompts user for output directory

listing = dir(strcat(workingdir, '*.png'));
numImages = length(listing);

%%  Load the images
for i = 1:numImages
    clc;
    disp(['Processing PNG image ', num2str(i), ' of ', num2str(numImages)]);
    title = listing(i,1).name(1:end-4);
    
    image = imread(strcat(workingdir, title, '.png'));

%%  Get the user to select ROIs for the zoomed sections
    % Display the image and a rectangular ROI that the user can move.
    h = imshow(image);
    sectionROI = images.roi.Rectangle(gca, 'Position',...
        [0, 0, floor(zoomedWidth / pixelSize)+1,...
        floor(zoomedHeight / pixelSize)+1]);
    while ishandle(h)
        % Continually store the ROI position until the user closes the
        % image.
        pos = round(sectionROI.Position);
        pause(0.02);
    end

%%  Export images and insets with added scale bars
    disp('Splitting and cropping the image...');
    % Isolate individual channels
    colloids = image(:, :, 1);
    gal8 = image(:, :, 2);
    nuclei = image(:, :, 3) - colloids;
    %   Crop the image to the inset.
    imageZoomed = image(...
        pos(2)+1:pos(2)+pos(4)-1, pos(1)+1:pos(1)+pos(3)-1, :);
    colloidsZoomed = imageZoomed(:, :, 1);
    gal8Zoomed = imageZoomed(:, :, 2);
    nucleiZoomed = imageZoomed(:, :, 3) - colloidsZoomed;
    
    largeHeight = height(nuclei);
    largeWidth = width(nuclei);
    smallHeight = floor(zoomedHeight/pixelSize);
    smallWidth = floor(zoomedWidth/pixelSize);
    
    %   Create images for the scale bar for both full and zoomed images.
    disp('Adding the scale bar...');
    zerosLarge = zeros(largeHeight, largeWidth, 'uint16');
    SBL = zerosLarge;
    zerosSmall = zeros(smallHeight, smallWidth, 'uint16');
    SBS = zerosSmall;
    SBL(largeHeight-50:largeHeight-25,...
        largeWidth-25-round(scaleBar/pixelSize):largeWidth-25) = 65535;
    SBS(smallHeight-50:smallHeight-25,...
        smallWidth-25-round(scaleBar/pixelSize):smallWidth-25) = 65535;
    
    
    %   Write the images.
    disp('Writing the images to disk...');
    exportBase = strcat(exportdir, title, '_');
    exportBaseZoomed = strcat(exportBase, 'zoomed_');
    imwrite(image + cat(3, SBL, SBL, SBL),...
        strcat(exportBase, 'composite.png'));
    imwrite(cat(3, colloids + SBL, SBL, colloids + SBL),...
        strcat(exportBase, 'colloids.png'));
    imwrite(cat(3, SBL, gal8 + SBL, SBL), strcat(exportBase, 'gal8.png'));
    imwrite(cat(3, SBL, SBL, nuclei + SBL),...
        strcat(exportBase, 'nuclei.png'));
    imwrite(imageZoomed + cat(3, SBS, SBS, SBS),...
        strcat(exportBaseZoomed, 'composite.png'));
    imwrite(cat(3, colloidsZoomed + SBS, SBS, colloidsZoomed + SBS),...
        strcat(exportBaseZoomed, 'colloids.png'));
    imwrite(cat(3, SBS, gal8Zoomed + SBS, SBS),...
        strcat(exportBaseZoomed, 'gal8.png'));
    imwrite(cat(3, SBS, SBS, nucleiZoomed + SBS),...
        strcat(exportBaseZoomed, 'nuclei.png'));
end
disp('Processing complete!');
