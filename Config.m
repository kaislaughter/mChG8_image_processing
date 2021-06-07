%% Image processing configuration

% Copy this file into the same folder as your input image files and modify
% it depending on your setup.

% The Run number is used to track multiple runs of the software, and is used in
% export file names and in the DataCells array. Note: it is a character
% array / string!

config = containers.Map();
config('RUN') = '1';

% Channel selection.
% If one channel is not included, set the channel to 0 or false.
config('GAL8_CHANNEL') = 1;
config('NUC_CHANNEL') = 2;
config('NP_CHANNEL') = 3;

% Choose file type for input microscope images
config('IMAGETYPE') = '*.CZI';

% Choose which images will be exported and in what format.
config('FILETYPE') = 'png';
config('EXPORT_NUC_MAP') = false;
config('EXPORT_INDIVIDUAL_CHANNELS') = false;
config('EXPORT_GAL8_ANNOTATIONS') = true;
config('EXPORT_NP_ANNOTATIONS') = true;
config('EXPORT_OVERLAP_ANNOTATIONS') = true;
config('EXPORT_COMPOSITE') = true;
config('EXPORT_CORRELATION_PLOTS') = true;
config('GAL8_BRIGHTEN') = 200;  % Signal multiplier for display only.
config('NP_BRIGHTEN') = 200;  % Signal multiplier for display only.
config('NUC_BRIGHTEN') = 50;  % Signal multiplier for display only.

% Define thresholds and correction factors for image processing
config('OPTIMIZE_GAL8_THRESHOLD') = false;
config('GAL8_THRESHOLD') = 100;
config('NP_THRESHOLD') = 80;
config('GAL8_NP_CROSSTALK') = 0.27;
config('NP_GAL8_CROSSTALK') = 0.01;

% Suppress sorting to test individual images
config('TEST_MODE') = false;

% Define technical replicates and plate size. The script assumes technical
% replicates are placed in the same column.
% Note that these values are for annotation purposes only; not configuring
% them will not affect quantification in any way.

if config('TEST_MODE')
    config('TECH_REPLICATES') = 1;
    config('PLATE_COLUMNS') = 1;
    config('GROUP_TITLES') = {'Test'};
else
    config('TECH_REPLICATES') = 3;
    config('PLATE_COLUMNS') = 10;
    config('GROUP_TITLES') = {...
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

% Define structural elements to be used in processing images
config('SE_GAL8_TH') = strel('disk', 20); % Gal8 tophat
config('SE_GAL8_OP') = strel('disk', 2); % Gal8 open
config('SE_NUC_OP') = strel('disk', 25); % Nucleus open
config('SE_NP_TH') = strel('disk', 5); % Uptake tophat
config('SE_NP_OP') = strel('disk', 2); % Uptake open