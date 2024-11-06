%% Script for preparing composites fohr publication

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
scaleBarS = 25;  % um
scaleBarL = 100; % um
pS = 25;  % pixels (padding and thickness of scale bar) - small
pL = 100; % pixels (padding and thickness of scale bar) - large
scaleBarLabel = false;
textSize = 24;  % point (font size for the scale bar label)
fileType = '.png'; % type of file to process

gal8Channel = false; % set to true if there's a gal8 channel
gal8Remove = false; % set to true if there is a gal8 channel to remove
plChannel = false;  % set to true if there's a 4th cyan channel
nucTHSize = 10;  % pixels (used to separate nuclei from PL channel)
% Note: the nuclei channel is not perfectly separated from the PL channel

separateChannels = true;
numChannels = 2;
nucleiMultiplier=1;
NPMultiplier=2;
Gal8Multiplier=2;

%% Setup
disp('Choose input directory')
workingdir = [uigetdir(), filesep]; % Prompts user for input directory
disp('Choose output directory')
exportdir = [uigetdir(), filesep]; % Prompts user for output directory

listing = dir(strcat(workingdir,['*',fileType]));
if separateChannels == true
    numImages = length(listing)/(numChannels+1);
else
    numImages = length(listing);
end

for i = 1:numImages
    if separateChannels == true
        a = (numChannels+1)*i-(numChannels);
    else
        a=i;
    end
    clc;
    disp(['Processing image ', num2str(i), ' of ', num2str(numImages)]);
    title = listing(a,1).name(1:end-4);
    disp(title);
    %% Load the image
    if separateChannels == true
        if gal8Channel == true
            title = listing(a,1).name(1:end-4);
            imageGal8 = imread(strcat(workingdir, title, fileType));
            title = listing(a+1,1).name(1:end-4);
            imageNP = imread(strcat(workingdir, title, fileType));
            title = listing(a+2,1).name(1:end-4);
            imageNuc = imread(strcat(workingdir, title, fileType));
            title = listing(a+3,1).name(1:end-4);
            image = imread(strcat(workingdir, title, fileType));
        else
            title = listing(a,1).name(1:end-4);
            imageNP = imread(strcat(workingdir, title, fileType));
            title = listing(a+1,1).name(1:end-4);
            imageNuc = imread(strcat(workingdir, title, fileType));
            title = listing(a+2,1).name(1:end-4);
            image = imread(strcat(workingdir, title, fileType));
        end
    else
        image = imread(strcat(workingdir, title, fileType));
    end

    largeHeight = height(image);
    largeWidth = width(image);
    zerosLarge = zeros(largeHeight, largeWidth, 'uint16');

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
    
    if separateChannels == true
        colloids = imageNP(:,:,1)*NPMultiplier;
        nuclei = imageNuc(:,:,3)*nucleiMultiplier;
        if gal8Channel == true
            gal8 = imageGal8(:,:,2)*Gal8Multiplier;
            image = cat(3,colloids,gal8,nuclei+colloids);
        else
            image = cat(3,colloids,zerosLarge,nuclei+colloids);
        end
    else
        colloids = image(:, :, 1);
        if gal8Channel == true
            if plChannel == true
                BTH = imtophat(image(:, :, 3), strel('disk', nucTHSize));
                extra = BTH - colloids;
                gal8 = image(:, :, 2) - extra;
                nuclei = image(:, :, 3) - BTH;
            else
                gal8 = image(:, :, 2);
                nuclei = image(:, :, 3) - colloids;
            end
        else
            nuclei = image(:, :, 3) - colloids;
        end
    end
    
    % takes out gal8 channel if gal8Remove = true
    if gal8Remove == true
        image = image - gal8;
    end
    
    % Crop the image to the inset.
    imageZoomed = image(...
        pos(2)+1:pos(2)+pos(4)-1, pos(1)+1:pos(1)+pos(3)-1, :);
    colloidsZoomed = colloids(...
        pos(2)+1:pos(2)+pos(4)-1, pos(1)+1:pos(1)+pos(3)-1, :);
    nucleiZoomed = nuclei(...
        pos(2)+1:pos(2)+pos(4)-1, pos(1)+1:pos(1)+pos(3)-1, :);
    if gal8Channel == true
        gal8Zoomed = gal8(...
            pos(2)+1:pos(2)+pos(4)-1, pos(1)+1:pos(1)+pos(3)-1, :);
    end
    if plChannel == true
        extraZoomed = extra(...
            pos(2)+1:pos(2)+pos(4)-1, pos(1)+1:pos(1)+pos(3)-1, :);
    end
    
    largeHeight = height(nuclei);
    largeWidth = width(nuclei);
    smallHeight = floor(zoomedHeight/pixelSize);
    smallWidth = floor(zoomedWidth/pixelSize);
    scaleBarWidthS = floor(scaleBarS/pixelSize);
    scaleBarWidthL = floor(scaleBarL/pixelSize);
    
    % Create images for the scale bar for both full and zoomed images.
    disp('Adding the scale bar...');
    zerosLarge = zeros(largeHeight, largeWidth, 'uint16');
    SBL = zeros(largeHeight, largeWidth, 3, 'uint16');
    zerosSmall = zeros(smallHeight, smallWidth, 'uint16');
    SBS = zeros(smallHeight, smallWidth, 3, 'uint16');
    % Create a white rectangle for the scale bar.
    SBL(largeHeight-2*pL:largeHeight-pL,...
        largeWidth-pL-scaleBarWidthL:largeWidth-pL, :) = 65535;
    SBS(smallHeight-2*pS:smallHeight-pS,...
        smallWidth-pS-scaleBarWidthS:smallWidth-pS, :) = 65535;
    % Add a text label to the scale bars.
    if scaleBarLabel == true
        SBL = insertText(SBL,...
            [largeWidth-pL-floor(scaleBarWidthL/2), largeHeight-2*pL],...
            [num2str(scaleBarS) ' µm'], 'TextColor', 'white',...
            'BoxOpacity', 0, 'FontSize', textSize,...
            'AnchorPoint', 'CenterBottom');
        SBS = insertText(SBS,...
            [smallWidth-pS-floor(scaleBarWidthS/2), smallHeight-2*pS],...
            [num2str(scaleBarS) ' µm'], 'TextColor', 'white',...
            'BoxOpacity', 0, 'FontSize', textSize,...
            'AnchorPoint', 'CenterBottom');  
    end
    
    % Write the images.
    disp('Writing the images to disk...');
    % Export full-size image.
    exportBase = strcat(exportdir, title, '_');
    imwrite(image + SBL, strcat(exportBase, 'composite.png'));
    imwrite(cat(3, colloids, zerosLarge, colloids) + SBL,...
        strcat(exportBase, 'NPs.png'));
    
    imwrite(cat(3, zerosLarge, zerosLarge, nuclei) + SBL,...
        strcat(exportBase, 'nuclei.png'));
    
    if gal8Channel == true && gal8Remove == false
        imwrite(cat(3, zerosLarge, gal8, zerosLarge) + SBL,...
            strcat(exportBase, 'gal8.png'));
    end
    
    if plChannel == true
        imwrite(cat(3, zerosLarge, extra, extra) + SBL,...
            strcat(exportBase, 'extra.png'));
    end
    % Export zoomed-in image.
    exportBaseZoomed = strcat(exportBase, 'zoomed_');
    imwrite(imageZoomed + SBS, strcat(exportBaseZoomed, 'composite.png'));
    imwrite(cat(3, colloidsZoomed, zerosSmall, colloidsZoomed) + SBS,...
        strcat(exportBaseZoomed, 'NPs.png'));
    
    imwrite(cat(3, zerosSmall, zerosSmall, nucleiZoomed) + SBS,...
        strcat(exportBaseZoomed, 'nuclei.png'));
    
    if gal8Channel == true && gal8Remove == false
        imwrite(cat(3, zerosSmall, gal8Zoomed, zerosSmall) + SBS,...
            strcat(exportBaseZoomed, 'gal8.png'));
    end
    
    if plChannel == true
        imwrite(cat(3, zerosSmall, extraZoomed, extraZoomed) + SBS,...
            strcat(exportBaseZoomed, 'extra.png'));
    end
end
disp('Processing complete!');
