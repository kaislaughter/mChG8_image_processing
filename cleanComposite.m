%% Script for preparing composites for publication

%   Requirements:
%    - Image Processing Toolbox
%    - Machine Vision Toolbox

%   Input: exported composite images from EE_disruption_uptake_CZI.m
%   (or similar).
%   IMPORTANT: the image format must conform to the following standard:
%    - Nuclei: blue
%    - Gal8: green
%    - Colloids/NPs: magenta (equal red and blue)
%    - (Optional) phospholipidosis vesicles: cyan (equal blue and green)

%   Output: individual channels plus composite images, in both the original
%   size and a zoomed-in section, with scale bars.

clc, clear, close all;

%% Config
% Fill in this section before each run.

pixelSize = 0.114;  % um
zoomedWidth = 100;  % um
zoomedHeight = 100;  % um
scaleBar = 25;  % um
p = 25;  % pixels (padding and thickness of scale bar)
textSize = 48;  % point (font size for the scale bar label)

extraChannel = false;  % set to true if there's a 4th cyan channel

%% Setup
disp('Choose input directory')
workingdir = [uigetdir(), filesep]; % Prompts user for input directory
disp('Choose output directory')
exportdir = [uigetdir(), filesep]; % Prompts user for output directory

listing = dir(strcat(workingdir, '*.png'));
numImages = length(listing);

for i = 1:numImages
    clc;
    disp(['Processing image ', num2str(i), ' of ', num2str(numImages)]);
    title = listing(i,1).name(1:end-4);
    %% Load the image
    image = imread(strcat(workingdir, title, '.png'));

    %% Get the user to select ROIs for the zoomed sections
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

    %% Export images and insets with added scale bars
    disp('Splitting and cropping the image...');
    % Isolate individual channels
    colloids = image(:, :, 1);
    if extraChannel == true
        extra = image(:, :, 2) - colloids;
        gal8 = image(:, :, 2) - extra;
        nuclei = image(:, :, 3) - colloids - extra;
    else
        gal8 = image(:, :, 2);
        nuclei = image(:, :, 3) - colloids;
    end
    
    % Crop the image to the inset.
    imageZoomed = image(...
        pos(2)+1:pos(2)+pos(4)-1, pos(1)+1:pos(1)+pos(3)-1, :);
    colloidsZoomed = colloids(...
        pos(2)+1:pos(2)+pos(4)-1, pos(1)+1:pos(1)+pos(3)-1, :);
    gal8Zoomed = gal8(...
        pos(2)+1:pos(2)+pos(4)-1, pos(1)+1:pos(1)+pos(3)-1, :);
    nucleiZoomed = nuclei(...
        pos(2)+1:pos(2)+pos(4)-1, pos(1)+1:pos(1)+pos(3)-1, :);
    if extraChannel == true
        extraZoomed = extra(...
            pos(2)+1:pos(2)+pos(4)-1, pos(1)+1:pos(1)+pos(3)-1, :);
    end
    
    largeHeight = height(nuclei);
    largeWidth = width(nuclei);
    smallHeight = floor(zoomedHeight/pixelSize);
    smallWidth = floor(zoomedWidth/pixelSize);
    scaleBarWidth = floor(scaleBar/pixelSize);
    
    % Create images for the scale bar for both full and zoomed images.
    disp('Adding the scale bar...');
    zerosLarge = zeros(largeHeight, largeWidth, 'uint16');
    SBL = zeros(largeHeight, largeWidth, 3, 'uint16');
    zerosSmall = zeros(smallHeight, smallWidth, 'uint16');
    SBS = zeros(smallHeight, smallWidth, 3, 'uint16');
    % Create a white rectangle for the scale bar.
    SBL(largeHeight-2*p:largeHeight-p,...
        largeWidth-p-scaleBarWidth:largeWidth-p, :) = 65535;
    SBS(smallHeight-2*p:smallHeight-p,...
        smallWidth-p-scaleBarWidth:smallWidth-p, :) = 65535;
    % Add a text label to the scale bars.
    SBL = insertText(SBL,...
        [largeWidth-p-floor(scaleBarWidth/2), largeHeight-2*p],...
        [num2str(scaleBar) ' µm'], 'TextColor', 'white',...
        'BoxOpacity', 0, 'FontSize', textSize,...
        'AnchorPoint', 'CenterBottom');
    SBS = insertText(SBS,...
        [smallWidth-p-floor(scaleBarWidth/2), smallHeight-2*p],...
        [num2str(scaleBar) ' µm'], 'TextColor', 'white',...
        'BoxOpacity', 0, 'FontSize', textSize,...
        'AnchorPoint', 'CenterBottom');  
    
    % Write the images.
    disp('Writing the images to disk...');
    % Export full-size image.
    exportBase = strcat(exportdir, title, '_');
    imwrite(image + SBL, strcat(exportBase, 'composite.png'));
    imwrite(cat(3, colloids, zerosLarge, colloids) + SBL,...
        strcat(exportBase, 'colloids.png'));
    imwrite(cat(3, zerosLarge, gal8, zerosLarge) + SBL,...
        strcat(exportBase, 'gal8.png'));
    imwrite(cat(3, zerosLarge, zerosLarge, nuclei) + SBL,...
        strcat(exportBase, 'nuclei.png'));
    if extraChannel == true
        imwrite(cat(3, zerosLarge, extra, extra) + SBL,...
            strcat(exportBase, 'extra.png'));
    end
    % Export zoomed-in image.
    exportBaseZoomed = strcat(exportBase, 'zoomed_');
    imwrite(imageZoomed + SBS, strcat(exportBaseZoomed, 'composite.png'));
    imwrite(cat(3, colloidsZoomed, zerosSmall, colloidsZoomed) + SBS,...
        strcat(exportBaseZoomed, 'colloids.png'));
    imwrite(cat(3, zerosSmall, gal8Zoomed, zerosSmall) + SBS,...
        strcat(exportBaseZoomed, 'gal8.png'));
    imwrite(cat(3, zerosSmall, zerosSmall, nucleiZoomed) + SBS,...
        strcat(exportBaseZoomed, 'nuclei.png'));
    if extraChannel == true
        imwrite(cat(3, zerosSmall, extraZoomed, extraZoomed) + SBS,...
            strcat(exportBaseZoomed, 'extra.png'));
    end
end
disp('Processing complete!');
