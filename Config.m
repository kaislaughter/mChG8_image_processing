%% Image processing configuration

% Copy this file into the same folder as your input image files and modify
% it depending on your setup.

% The Run number is used to track multiple runs of the software, and is used in
% export file names and in the DataCells array. Note: it is a character
% array / string!
RUN = '1';

% Channel selection.
% If one channel is not included, set the channel to 0 or false.
GAL8_CHANNEL = 1;
NUC_CHANNEL = 2;
NP_CHANNEL = 3;

% Choose file type for input microscope images
IMAGETYPE = '*.CZI';

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
    GROUP_TITLES = {'Test'};
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