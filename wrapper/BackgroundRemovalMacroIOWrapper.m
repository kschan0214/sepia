%% [localField,maskFinal] = BackgroundRemovalMacroIOWrapper(input,output,maskFullName,algorParam)
%
% Input
% --------------
% input         : input directory contains NIfTI (*totalfield* and
%                 *fieldmapsd*) files or structure containing filenames
% output        : output directory that stores the output (local field and final mask)
% maskFullName  : mask filename
% algorParam    : structure contains method and method specific parameters
%
% Output
% --------------
% localField    : local field (or tissue field) (in Hz)
% maskFinal     : final mask used for QSM
%
% Description:  This is a wrapper of BackgroundRemovalMacro.m for NIfTI
%               input/output
%
% Kwok-shing Chan @ DCCN
% k.chan@donders.ru.nl
% Date created: 17 April 2018
% Date modified: 26 August 2018
% Date modified: 29 March 2019
% Date modified: 8 March 2020 (v0.8.0)
% Date modified: 21 Jan 2020 (v0.8.1)
% Date modified: 13 August 2021 (v1.0)
%
%
function [localField,mask_QSM] = BackgroundRemovalMacroIOWrapper(input,output,maskFullName,algorParam)
%% add general Path
sepia_addpath;

sepia_universal_variables;

%% define variables
prefix = 'sepia_';

%% Check if output directory exists 
output_index    = strfind(output, filesep);
outputDir       = output(1:output_index(end));
% get prefix
if ~isempty(output(output_index(end)+1:end))
    prefix = [output(output_index(end)+1:end) '_'];
end
% if the output directory does not exist then create the directory
if exist(outputDir,'dir') ~= 7
    mkdir(outputDir);
end

% display output info
fprintf('Output directory       : %s\n',outputDir);
fprintf('Output filename prefix : %s\n',prefix);

outputFileList = construct_output_filename(outputDir, prefix, suffix);

%% Setting up Input
disp('---------');
disp('Load data');
disp('---------');

%%%%%% Step 1: get all required filenames
% input         : can be input directory or structure contains input filenames
% outputDir     : output directory (only for BIDS)
% prefix        : output basename (only for BIDS)
% inputDir      : intput directory of phase image
% inputFileList : structure contains all input filenames
[inputDir, inputFileList]	= io_01_get_input_file_list(input, outputDir, prefix);

%%%%% Step 2: validate input files
% 2.2 validate nifti files
% inputFileList         : structure contains all input filenames
% availableFileList     : data that is already available and validated
availableFileList           = io_02_validate_nifti_input(inputFileList);

%%%%% Step 3: get nifti template header for exporting output data
% availableFileList   	: structure contains all data filenames that are already available and validated
% outputNiftiTemplate   : nifti header with empty 'img' field
outputNiftiTemplate         = io_03_get_nifti_template(availableFileList);

% 3.2 load and validate SEPIA header 
if ~isempty(inputFileList(4).name)
    sepia_header = load([inputFileList(4).name]);
    disp('SEPIA header data is loaded.');
    % Validate header information
    sepia_header = validate_sepia_header_4wrapper(sepia_header, outputNiftiTemplate);

else
    error('Please specify a header required by SEPIA.');
end

% display some header info
display_sepia_header_info_4wrapper;

%%%%%% Step 4: get signal mask
disp('-----------');
disp('Signal mask');
disp('-----------');
% maskFullName          : mask filename
% inputDir              : intput directory of phase image
% sepia_header          : sepia header
% algorParam            : structure contains all pipeline parameters
% availableFileList     : structure contains all data filenames that are already available and validated
% outputFileList        : structure contains default output filenames
% outputNiftiTemplate   : nifti header with empty 'img' field
availableFileList           = io_04_get_signal_mask(maskFullName, inputDir, sepia_header, algorParam, availableFileList, outputFileList, outputNiftiTemplate);

%%%%%% Step 5: store some data to headerAndExtraData
% header
create_header_structure_4wrapper;

matrixSize  = double(sepia_header.matrixSize);
voxelSize   = double(sepia_header.voxelSize);

headerAndExtraData.availableFileList = availableFileList;
headerAndExtraData.outputDirectory   = outputDir; 

%% Background field removal
totalField   	= double(load_nii_img_only(availableFileList.totalField));
maskLocalfield	= double(load_nii_img_only(availableFileList.maskLocalField));

% core of background field removal
localField = BackgroundRemovalMacro(totalField,maskLocalfield,matrixSize,voxelSize,algorParam,headerAndExtraData);
clear totalField maskLocalfield % clear variables that no longer be needed

% generate new mask based on backgroudn field removal result
% mask_QSM = localField ~=0;
% 20230124 v1.2.2: make sure no holes inide ROIs
mask_QSM = imfill(localField ~= 0, 'holes');

% save results
fprintf('Saving local field map...');
save_nii_quick(outputNiftiTemplate,localField, outputFileList.localField);
fprintf('done!\n');
availableFileList.localField = outputFileList.localField;

fprintf('Saving mask for chi mapping...');
save_nii_quick(outputNiftiTemplate,mask_QSM, outputFileList.maskQSM);
fprintf('done!\n');
availableFileList.maskQSM = outputFileList.maskQSM;

disp('Processing pipeline is completed!');

end

%% Validate input data dimension
function check_input_dimension(img,matrixSize,TE)

matrixSize_img	= size(img);

% check matrix size between total field/fieldmap sd input and SEPIA header
if ~isequal(matrixSize_img(1:3),matrixSize)
    erro('Input NIfTI files and SEPIA header do not have the same matrix size. Please check these files and/or remove the ''matrixSize'' variable from the SEPIA header.')
end

% check echo dimension 
if ndims(img) == 4
    if matrixSize_img(4) ~= length(TE)
        error('Input NIfTI file and SEPIA header do not have the same number of echoes.  Please check these files.');
    end
end

end

%% I/O Step 1: get input file list
function [inputDir, inputNiftiList] = io_01_get_input_file_list(input,outputDir,prefix)

if isstruct(input)
    
    % Option 1: input are files
    inputNiftiList  = input;
    
    % take the phase data directory as reference input directory 
    [inputDir,~,~] = fileparts(inputNiftiList(1).name);
    
else
    
    % Option 2: input is a directory
    inputDir = input; 
    
    % First check with SEPIA default naming structure
    disp('Searching input directory based on SEPIA default naming structure...');
    filePattern = {'fieldmap','','noisesd','header'}; % don't change the order
    [isLoadSuccessful, inputNiftiList] = read_default_to_filelist(inputDir, filePattern);
    
    % If it doesn't work then check BIDS compatibility
    if ~isLoadSuccessful
        error('No file matches the SEPIA default naming structure. Please check the input directory again.')
    end
    
end

end

%% I/O Step 2: validate input data
function availableFileList          = io_02_validate_nifti_input(inputFileList)

availableFileList = struct();

% load NIFTI header for validating input images
% 2.1 Total field map 
if ~isempty(inputFileList(1).name)
    
    fprintf('Loading total field header files...')
    
    % load nifti structure for output template
    TotalFieldNIFTIHeader = load_untouch_header_only(inputFileList(1).name);
    % store it in the availableFileList
    availableFileList.totalField = inputFileList(1).name;
    
    fprintf('Done.\n');
    
else
    error('Fail! \nPlease specify a 3D total field map.');
end

% 2.2 Fieldmapsd data 
if ~isempty(inputFileList(3).name)
    
    fprintf('Loading fieldmap standard deviation header files...')
    
    % load true value from NIfTI
    fieldmapSDNIFTIHeader = load_untouch_header_only(inputFileList(3).name);
    % store it in the availableFileList
    availableFileList.fieldmapSD = inputFileList(3).name;

    fprintf('Done.\n');
    
else
    disp('No field map standard deviation data is loaded.');
end

disp('Input files are valid.')

end

%% I/O Step 3: get nifti template for nifti output
function outputNiftiTemplate        = io_03_get_nifti_template(availableFileList)

outputNiftiTemplate     = load_untouch_nii(availableFileList.totalField);
outputNiftiTemplate.img = [];

end

%% I/O Step 4: loading signal mask
function availableFileList          = io_04_get_signal_mask(maskFullName, inputDir, sepia_header, algorParam, availableFileList, outputFileList, outputNiftiTemplate)

matrixSize  = sepia_header.matrixSize;

% check if there is any mask specifically for BFR first
maskList = dir(fullfile(inputDir,'*mask_localfield*nii*'));
% if not then look for a general mask
if isempty(maskList)
    maskList = dir(fullfile(inputDir,'*mask*nii*'));
end

% Scenario: No specified mask file + there is a file called mask in the input directory
if isempty(maskFullName) && ~isempty(maskList)
    
    fprintf('No mask file is specified but a mask file is found in the input directory: %s\n',fullfile(inputDir, maskList(1).name));
    disp('Trying to load the file as signal mask');
    
    maskFullName = fullfile(inputDir, maskList(1).name);
end

% Scenario: User provided a mask file or above scenario was satified
if ~isempty(maskFullName)
    
    % load mask file
    mask = load_nii_img_only(maskFullName) > 0;
    
    % make sure the mask has the same dimension as other input data
    if ~isequal(size(mask),matrixSize)
        disp('The file does not have the same dimension as other images.')
        % display error message 
        error('No mask file is loaded. Please specify a valid mask file or put it in the input directory.');
    else
        availableFileList.maskLocalField = maskFullName;
        disp('Mask file is checked.');
    end
    
else
    % display error message if nothing is found
    error('No mask file is found. Please specify your mask file or put it in the input directory.');
end

end
