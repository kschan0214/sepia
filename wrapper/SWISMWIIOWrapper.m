%% chi = QSMMacroIOWrapper(input,output,maskFullName,algorParam)
%
% Input
% --------------
% input         :   input directory contains NIfTI (*localfield*, *magn* and
%                   *fieldmapsd*) files or structure containing filenames  
% output        :   output directory that stores the output (susceptibility map)
% maskFullName  :   mask filename
% algorParam    :   structure contains method and method specific parameters
%
% Output
% --------------
% chi           : magnetic susceptibility map (in ppm)
%
% Description: This is a wrapper of QSMMacro.m which for NIfTI input/output
%
% Kwok-shing Chan @ DCCN
% k.chan@donders.ru.nl
% Date created: 17 April 2018
% Date modified: 26 August 2018
% Date modified: 29 March 2019
% Date modified: 5 June 2019
% Date modified: 8 March 2020 (v0.8.0)
% Date modified: 21 Jan 2020 (v0.8.1)
% Date modified: 13 August 2021 (v1.0)
%
%
function SWISMWIIOWrapper(input,output,algorParam)
%% add general Path
sepia_addpath

sepia_universal_variables;

%% define variables
prefix = 'sepia_';

%% Check if output directory exists 
output_index = strfind(output, filesep);
outputDir = output(1:output_index(end));
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

% outputFileList = construct_output_filename(outputDir, prefix);
outputFileFullPrefix = fullfile(outputDir,prefix);

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
[~, inputFileList]	= io_01_get_input_file_list(input);

%%%%% Step 2: validate input files
% 2.2 validate nifti files
% inputFileList         : structure contains all input filenames
% algorParam            : structure contains all pipeline configuration
% availableFileList     : data that is already available and validated
availableFileList           = io_02_validate_nifti_input(inputFileList);

%%%%% Step 3: get nifti template header for exporting output data
% availableFileList   	: structure contains all data filenames that are already available and validated
% outputNiftiTemplate   : nifti header with empty 'img' field
outputNiftiTemplate         = io_03_get_nifti_template(availableFileList);

% 3.2 load and validate SEPIA header 
% if ~isempty(inputFileList(3).name)
if numel(inputFileList) > 2 && ~isempty(inputFileList(3).name)
    sepia_header = load([inputFileList(3).name]);
    disp('SEPIA header data is loaded.');
    % Validate header information
    sepia_header = validate_sepia_header_4wrapper(sepia_header, outputNiftiTemplate);

else
    
    warning('No SEPIA header is provided. Some method may require additional header info from the SEPIA header file. Initiate header with default values.');
    sepia_header.B0 = 3;
    sepia_header.TE = 15e-3;
    sepia_header = validate_sepia_header_4wrapper(sepia_header, outputNiftiTemplate);
end

% display some header info
display_sepia_header_info_4wrapper;

%%%%%% Step 5: store some data to headerAndExtraData
% header
create_header_structure_4wrapper;

% matrixSize  = double(sepia_header.matrixSize);
% voxelSize   = double(sepia_header.voxelSize);

headerAndExtraData.availableFileList    = availableFileList;
output_structure.outputFileFullPrefix   = outputFileFullPrefix;
output_structure.outputNiftiTemplate    = outputNiftiTemplate;

%% SWI or SMWI
% Unlike other QSM processing steps, loading of the input will be done
% inside the wrapper function because SWI and SMWI have different input
% data (i.e. phase and chi map). Export of the output will also be done
% inside the wrapper function

% General steps as follow in the wrapper function:
% 1. get the required input from 
% 2. main QSM algorithm
% 3. convert output unit to ppm
method = algorParam.swismwi.method;
for k = 1:length(wrapper_SWISMWI_function)
    if strcmpi(method,methodSWISMWIName{k})
        feval(wrapper_SWISMWI_function{k},headerAndExtraData,output_structure,algorParam);
    end
end


% localField   	= double(load_nii_img_only(availableFileList.localField));
% mask_QSM        = double(load_nii_img_only(availableFileList.maskQSM));

% % core of QSM
% [chi,mask_ref] = QSMMacro(localField,mask_QSM,matrixSize,voxelSize,algorParam,headerAndExtraData);
% clear localField mask_QSM
% 
% % save results
% fprintf('Saving susceptibility map...');
% save_nii_quick(outputNiftiTemplate, chi, outputFileList.QSM);
% clear chi
% 
% if ~isempty(mask_ref)
%     save_nii_quick(outputNiftiTemplate, mask_ref, outputFileList.maskRef);
% end
% 
% fprintf('Done!\n');

disp('Processing pipeline is completed!');

end

%% I/O Step 1: get input file list
function [inputDir, inputNiftiList] = io_01_get_input_file_list(input)

% if isstruct(input)
    
    % Option 1: input are files
    inputNiftiList  = input;
    
    % take the phase data directory as reference input directory 
    [inputDir,~,~] = fileparts(inputNiftiList(1).name);
    
% else
%     
%     % Option 2: input is a directory
%     inputDir = input; 
%     
%     % First check with SEPIA default naming structure
%     disp('Searching input directory based on SEPIA default naming structure...');
%     filePattern = {'localfield','mag','weights','header'}; % don't change the order
%     [isLoadSuccessful, inputNiftiList] = read_default_to_filelist(inputDir, filePattern);
%     
%     % If it doesn't work then check BIDS compatibility
%     if ~isLoadSuccessful
%         error('No file matches the SEPIA default naming structure. Please check the input directory again.')
%     end
%     
% end

end

%% I/O Step 2: validate input data
function availableFileList          = io_02_validate_nifti_input(inputFileList)

availableFileList = struct();

% load NIFTI header for validating input images
if ~isempty(inputFileList(1).name)
    
    fprintf('Loading Phase/QSM image header file...')
    
    % get header info from NIFTI for validation
    LocalFieldNIFTIHeader = load_untouch_header_only(inputFileList(1).name);
    
    % store the filename in the availableFileList structure
    availableFileList.phasechi = inputFileList(1).name;
    
    fprintf('Done.\n');
    
else
    error('Please specify a Phase image or a QSM map.');
end

% 2.2 magnitude data 
if ~isempty(inputFileList(2).name)
    
    fprintf('Loading magnitude image header file...')
    
    % get header info from NIFTI for validation
    magnitudeNIFTIHeader = load_untouch_header_only(inputFileList(2).name);
    
    % store the filename in the availableFileList structure
    availableFileList.magnitude = inputFileList(2).name;
    
    fprintf('Done.\n');
    
else
    disp('No magnitude data is loaded.');
end

disp('Input NIfTI files are valid.')

end

%% I/O Step 3: get nifti template for nifti output
function outputNiftiTemplate        = io_03_get_nifti_template(availableFileList)

outputNiftiTemplate     = load_untouch_nii(availableFileList.phasechi);
% just keep the header
outputNiftiTemplate.img = [];

end

% %% I/O Step 4: loading signal mask
% function availableFileList          = io_04_get_signal_mask(maskFullName, inputDir, sepia_header, algorParam, availableFileList, outputFileList, outputNiftiTemplate)
% 
% matrixSize  = sepia_header.matrixSize;
% 
% % check if there is any mask specifically for BFR first
% maskList = dir(fullfile(inputDir,'*mask_QSM*nii*'));
% % if not then look for a general mask
% if isempty(maskList)
%     maskList = dir(fullfile(inputDir,'*mask*nii*'));
% end
% 
% % Scenario: No specified mask file + there is a file called mask in the input directory
% if isempty(maskFullName) && ~isempty(maskList)
%     
%     fprintf('No mask file is specified but a mask file is found in the input directory: %s\n',fullfile(inputDir, maskList(1).name));
%     disp('Trying to load the file as signal mask');
%     
%     maskFullName = fullfile(inputDir, maskList(1).name);
% end
% 
% % Scenario: User provided a mask file or above scenario was satified
% if ~isempty(maskFullName)
%     
%     % load mask file
%     mask = load_nii_img_only(maskFullName) > 0;
%     
%     % make sure the mask has the same dimension as other input data
%     if ~isequal(size(mask),matrixSize)
%         disp('The file does not have the same dimension as other images.')
%         % display error message 
%         error('No mask file is loaded. Please specify a valid mask file or put it in the input directory.');
%     else
%         availableFileList.maskQSM = maskFullName;
%         disp('Mask file is checked.');
%     end
%     
% else
%     % display error message if nothing is found
%     error('No mask file is found. Please specify your mask file or put it in the input directory.');
% end
% 
% end