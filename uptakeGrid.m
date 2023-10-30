%% Script for preparing image grids for publication

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
p = 25;  % pixels (border width)
textSize = 48;  % point (font size for labels)
leftMargin = 10;
brightnessMultiplier = 3;

zoomed = false;
manualEntry = true;
nonUniformImages = true;

if zoomed == true
    nucName = 'zoomed_nuclei';
    collName = 'zoomed_colloids';
    compName = 'zoomed_composite';
else
    nucName = 'composite_nuclei';
    collName = 'composite_colloids';
    compName = 'composite_composite';
end


imageList = [130,65,125,15,70,95];
%% Setup

disp('Choose directory')
workingdir = [uigetdir(), filesep]; % Prompts user for input directory
listing = dir(strcat(workingdir, '*.png'));

list_full = [];

if manualEntry == true
    numImages = length(imageList);
    for i = 1:length(imageList)
        list_full_temp = struct2table(listing);
        list_full = [list_full;list_full_temp(contains(list_full_temp.name,['(',num2str(imageList(i)),')']),:)];
    end
else
    numImages = length(listing);
    list_full = struct2table(listing);
end

    list_ZNuc = list_full(contains(list_full.name,nucName),:);
    list_ZColl = list_full(contains(list_full.name,collName),:);
    list_ZComp = list_full(contains(list_full.name,compName),:);

numGroups = height(list_ZNuc);
imOut = [];
p0 = p;

for i = 1:numImages
    clc;
    disp(['Processing group ', num2str(i), ' of ', num2str(numGroups)]);
    titleNuc = char(list_ZNuc(i,1).name);
    imNuc = imread(strcat(workingdir, titleNuc));
    titleColl = char(list_ZColl(i,1).name);
    imColl = imread(strcat(workingdir, titleColl));
    titleComp = char(list_ZComp(i,1).name);
    imComp = imread(strcat(workingdir, titleComp));
    
    if nonUniformImages == true
        p=p0 + 5300-size(imNuc,2);
    end
    
    imOut = [imOut;...
        ones(size(imNuc,1),leftMargin+p,3)*65535,...
        imNuc,...
        ones(size(imNuc,1),p,3)*65535,...
        imColl,...
        ones(size(imNuc,1),p,3)*65535,...
        imComp];
    imOut = insertText(imOut,...
        [leftMargin-200,(2*i-1)*(size(imNuc,1)+p)/2+p],...
        titleComp(1:end-33),'FontSize',100,'AnchorPoint','RightBottom',...
        'BoxOpacity',0,'TextColor','black');
    imOut = [imOut;...
        ones(p,size(imOut,2),3)*65535];

end

imshow(imOut*brightnessMultiplier)
imwrite(imOut*brightnessMultiplier,'grid_images.png');
