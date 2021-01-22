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
%
%
function [localField,maskFinal] = BackgroundRemovalMacroIOWrapper(input,output,maskFullName,algorParam)
%% add general Path
sepia_addpath;

sepia_universal_variables;

%% define variables
prefix = 'sepia_';

% make sure the input only load once (first one)
isTotalFieldLoad = false;
isFieldmapSDLoad = false;

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

%% Setting up Input
disp('---------');
disp('Load data');
disp('---------');

%%%%%% Step 1: get all required filenames
if isstruct(input)
    
    % Option 1: input are files
    inputNiftiList = input;
    
    % take the total field map directory as reference input directory 
    [inputDir,~,~] = fileparts(inputNiftiList(1).name);
    
else
    
    % Option 2: input is a directory
    disp('Searching input directory...');
    inputDir = input;
    
    % check and get filenames
    inputNiftiList = struct();
    filePattern = {'total-field','','noise-sd','header'}; % don't change the order
    for  k = 1:length(filePattern)
        if k ~= 2
            % get filename
            if k ~= 4   % NIFTI image input
                [file,numFiles] = get_filename_in_directory(inputDir,filePattern{k},'.nii');
            else        % SEPIA header
                [file,numFiles] = get_filename_in_directory(inputDir,filePattern{k},'.mat');
            end

            % actions given the number of files detected
            if numFiles == 1        % only one file -> get the name
                
                fprintf('One ''%s'' file is found: %s\n',filePattern{k},file.name);
                inputNiftiList(k).name = file.name;
                
            elseif numFiles == 0     % no file -> error
                
                if k ~= 3 % essential file 'total-field'
                    error(['No file with name containing string ''' filePattern{k} ''' is detected.']);
                else
                    disp(['No file with name containing string ''' filePattern{k} ''' is detected.']);
                end
                
            else % multiple files -> fatal error
                
                error(['Multiple files with name containing string ''' filePattern{k} ''' are detected. Make sure the input directory should contain only one file with string ''' filePattern{k} '''.']);
                
            end
        else
            inputNiftiList(k).name = '';
        end
    end
    
end

%%%%%% Step 2: load data
% 2.1 Total field map 
if ~isempty(inputNiftiList(1).name)
    
    % load nifti structure for output template
    inputTotalFieldNifti = load_untouch_nii([inputNiftiList(1).name]);
    % load true value from NIfTI
    totalField = load_nii_img_only(inputNiftiList(1).name);
    isTotalFieldLoad = true;
    disp('Total field map is loaded.')
    
else
    error('Please specify a 3D total field map.');
end

% 2.2 Fieldmapsd data 
if ~isempty(inputNiftiList(3).name)
    
    % load true value from NIfTI
    fieldmapSD = load_nii_img_only([inputNiftiList(3).name]);
    isFieldmapSDLoad = true;
    disp('Field map standard deviation (noise SD) data is loaded.')
    
else
    disp('No field map standard deviation data is loaded.');
end

%2.3 SEPIA header
if ~isempty(inputNiftiList(4).name)
    
    load([inputNiftiList(4).name]);
    disp('Header data is loaded.');
    
else
    error('Please specify a header required by SEPIA.');
end

% store the header of the NIfTI files, all following results will have
% the same header
outputNiftiTemplate = inputTotalFieldNifti;

%%%%%% Step 3: validate input
% 3.1 Validate header information
validate_sepia_header;

% 3.2 Validate NIfTI input
disp('Validating input NIfTI files...')
% make sure the dimension of input data is consistent
check_input_dimension(totalField,matrixSize,TE);
if exist('fieldmapSD','var'); check_input_dimension(fieldmapSD,matrixSize,TE); end
disp('Input NIfTI files are valid.')

% display some header info
display_sepia_header_info_4wrapper;

% ensure variables are double
totalField  = double(totalField);
voxelSize   = double(voxelSize);
matrixSize  = double(matrixSize);
if exist('fieldmapSD','var'); fieldmapSD = double(fieldmapSD); end

%%%%%% Step 5: store some data to headerAndExtraData
% header
create_header_structure_4wrapper;

if exist('fieldmapSD','var'); headerAndExtraData.N_std = fieldmapSD; end

clearvars inputTotalFieldNifti

%% get brain mask
disp('-----------');
disp('Signal mask');
disp('-----------');

% check if there is any mask specifically for BFR first
maskList = dir(fullfile(inputDir,'*mask-local_field*nii*'));
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
        disp('Mask file is loaded.');
    end
    
else
    % display error message if nothing is found
    error('No mask file is found. Please specify your mask file or put it in the input directory.');
end

mask = double(mask);

%% Background field removal
% core of background field removal
localField = BackgroundRemovalMacro(totalField,mask,matrixSize,voxelSize,algorParam,headerAndExtraData);
  
% generate new mask based on backgroudn field removal result
maskFinal = localField ~= 0;
  
% save results
fprintf('Saving local field map...');
save_nii_quick(outputNiftiTemplate,localField, [outputDir filesep prefix 'local-field.nii.gz']);
save_nii_quick(outputNiftiTemplate,maskFinal,  [outputDir filesep prefix 'mask-qsm.nii.gz']);
fprintf('Done!\n');

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
