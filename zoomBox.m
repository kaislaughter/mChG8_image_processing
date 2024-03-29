%% Script for preparing images for publication with zoom box

%   Requirements:
%    - Image Processing Toolbox
%    - Machine Vision Toolbox
%    - Computer Vision Toolbox

%   Input: exported images from cleanComposite.m
%   (or similar).

%   Output: images with zoomed box in the top right corner

clc, clear, close all;

%% Config
% Fill in this section before each run.
% size of zoomed in area
zoomAreaW = 50;
zoomAreaH = 50;

% size of box in top right corner
LboxW = 300;
LboxH = 300;

% defining thickensses for annotations
boxBorder = 4;
lineThick = 2;

% Optionally select only a subset of the images in the input directory
%   This option is useful if there are multiple image categories in the
%   input folder.
%   Example: titleFilter = 'composite_zoomed_composite';
%   Use an empty string (i.e., '') to process all images.
titleFilter = '';

% Define an identifier for the output images that is added to the end of
% the filename.
fileEnding = 'inset';

%% Setup
disp('Choose input directory')
workingdir = [uigetdir(), filesep]; % Prompts user for input directory
disp('Choose output directory')
exportdir = [uigetdir(), filesep]; % Prompts user for output directory
%
listing = dir(strcat(workingdir, '*.png'));
numImages = length(listing);

% loop over all images
for i = 1:numImages
    clc;
    disp(['Processing image ', num2str(i), ' of ', num2str(numImages)]);
    title = listing(i,1).name(1:end-4);
    
    %% Only processes composite images
    if contains(title, titleFilter)  
        
        %% Load the image
        image = imread(strcat(workingdir, title, '.png'));

        % define image width, height, and number of layers
        [imW,imH,imD] = size(image);

        Lbox = [imW-LboxW-boxBorder/2,1+boxBorder/2,LboxW,LboxH];

        %% Get the user to select ROIs for the zoomed sections
        % Display the image and a rectangular ROI that the user can move.
        h = imshow(image);
        sectionROI = images.roi.Rectangle(gca, 'Position',...
            [50,imH-zoomAreaH-50, zoomAreaW, zoomAreaH], ...
            'InteractionsAllowed','translate');
        % shows area where large box will be located
        sectionROI2 = images.roi.Rectangle(gca, 'Position',...
            Lbox,'InteractionsAllowed','none');
        while ishandle(h)
            % Continually store the ROI position until the user closes the
            % image.
            pos = round(sectionROI.Position);
            pause(0.02);
        end
        
        % define zoomed in area of image
        imageZoomed = image(pos(2)+1:pos(2)+pos(4)-1,...
            pos(1)+1:pos(1)+pos(3)-1, :);
        
        % define export filename base
        exportBase = strcat(exportdir, title, '_');
        
        % define variable for output image
        im_out = image;
        
        % set top right corner as zoomed in image
        im_out(1:LboxH,...
            end-LboxW+1:end,:) = ...
            imresize(imageZoomed,[LboxH,LboxW],'nearest');
        
        % determines which quadrant the selected region is in to position
        % the zoom lines accordingly
        L1 = [pos(1)+boxBorder/2,pos(2)-boxBorder/2,...
            imW-LboxW-boxBorder/2,1];
        L2 = [pos(1)+pos(3)-boxBorder/2,pos(2)+pos(4)+boxBorder/2,...
            imW,LboxH+1+boxBorder/2];
        if pos(1) > imW-LboxW
            L1(2) = pos(2)+pos(4)+boxBorder/2;
            L1(4) = LboxH+1+boxBorder/2;
        elseif pos(2) + pos(4) < LboxH
            L2(1) = pos(1)+boxBorder/2;
            L2(3) = imW-LboxW-boxBorder/2;
        end

        % draw rectange and line annotations
        im_out = insertShape(im_out,'Rectangle',[Lbox;pos],...
            'Color','w','LineWidth',boxBorder);
        im_out = insertShape(im_out,'Line',...
            [L1;L2],'Color','w','LineWidth',lineThick);
        
        % write output image
        imwrite(im_out,strcat(exportBase,fileEnding,'.png'));

    end
end
disp('Processing complete!');
