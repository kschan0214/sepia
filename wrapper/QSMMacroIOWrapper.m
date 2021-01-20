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
%
%
function chi = QSMMacroIOWrapper(input,output,maskFullName,algorParam)
%% add general Path
sepia_addpath

sepia_universal_variables;

%% define variables
prefix = 'sepia_';
isInputDir = true;
% make sure the input only load once (first one)
isLocalFieldLoad    = false;
isWeightLoad        = false;
isMagnLoad          = false; 

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

%% Read input
disp('Reading data...');

% Step 1: check input for nifti files first
if isstruct(input)
    % Option 1: input are files
    inputNiftiList = input;
    isInputDir = false;
    
    % take the phase data directory as reference input directory 
    [inputDir,~,~] = fileparts(inputNiftiList(1).name);
else
    % Option 2: input is a directory
    inputDir = input;
%     inputNiftiList = dir([inputDir '/*.nii*']);
    % check and get filenames
    inputNiftiList = struct();
    filePattern = {'local-field','mag','weights','header'}; % don't change the order
    for  k = 1:length(filePattern)
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
        elseif numFiles == 0     % no file -> fatal error
            if k ~= 2 && k ~= 3 % essential file 'local-field' 
                error(['No file with name containing string ''' filePattern{k} ''' is detected.']);
            else
                disp(['No file with name containing string ''' filePattern{k} ''' is detected.']);
            end
        else % multiple files -> fatal error
            error(['Multiple files with name containing string ''' filePattern{k} ''' are detected. Make sure the input directory should contain only one file with string ''' filePattern{k} '''.']);
        end
    end
end

% Step 2: load data
%%%%%%%%%% Local field map %%%%%%%%%% 
if ~isempty(inputNiftiList(1).name)
    inputLocalFieldNifti = load_untouch_nii([inputNiftiList(1).name]);
    % load true value from NIfTI
    localField = load_nii_img_only([inputNiftiList(1).name]);
    
    isLocalFieldLoad = true;

    disp('Local field map is loaded.')
else
    error('Please specify a 3D local field map.');
end
%%%%%%%%%% magnitude data %%%%%%%%%%
if ~isempty(inputNiftiList(2).name)
    
    % load true value from NIfTI
    magn = load_nii_img_only([inputNiftiList(2).name]);
%             magn = double(inputMagnNifti.img);
    isMagnLoad = true;
    disp('Magnitude data is loaded.');
else
    disp('No magnitude data is loaded.');
    magn = ones(size(localField));
end
%%%%%%%%%% weights map %%%%%%%%%%
if ~isempty(inputNiftiList(3).name)
    
    % load true value from NIfTI
    weights = load_nii_img_only([inputNiftiList(3).name]);
    
    % check whether phase data contains DICOM values or wrapped
    if size(weights,4) > 1
        error('Input weighting map is 4D. Please specify a 3D weight data.');
    end
    isWeightLoad = true;
    disp('Weights data is loaded');
else
    disp('Default weighting method will be used for QSM.');
end
%%%%%%%%%% SEPIA header %%%%%%%%%%
if ~isempty(inputNiftiList(4).name)
    load([inputNiftiList(4).name]);
    disp('Header data is loaded.');
else
    error('Please specify a header required by SEPIA.');
end

% if no magnitude is loaded and MEDI is chosen -> fatal error
if ~isMagnLoad && strcmpi(algorParam.qsm.method,'MEDI')
    error('MEDI requires magnitude data. Please put the magnitude multi-echo data to input directory or use other algorithm');
end
if ~isWeightLoad && strcmpi(algorParam.qsm.method,'MEDI')
    error('MEDI requires a weighting map. Please put a (SNR) weighting map to input directory or use other algorithm');
end

% store the header the NIfTI files, all following results will have
% the same header
outputNiftiTemplate = inputLocalFieldNifti;
clearvars inputLocalFieldNifti

% In case some parameters are missing in the header file
if ~exist('matrixSize','var')
    matrixSize = size(magn);
    matrixSize = matrixSize(1:3);
end
if ~exist('voxelSize','var')
    voxelSize = outputNiftiTemplate.hdr.dime.pixdim(2:4);
end


% make sure the L2 norm of B0 direction = 1
B0_dir = B0_dir ./ norm(B0_dir);

% display some header info
disp('----------------------');
disp('Basic Data information');
disp('----------------------');
disp(['Voxel size(x,y,z mm^3)   =  ' num2str(voxelSize(1)) 'x' num2str(voxelSize(2)) 'x' num2str(voxelSize(3))]);
disp(['matrix size(x,y,z)       =  ' num2str(matrixSize(1)) 'x' num2str(matrixSize(2)) 'x' num2str(matrixSize(3))]);
disp(['B0 direction(x,y,z)      =  ' num2str(B0_dir(:)')]);
disp(['Field strength(T)        =  ' num2str(B0)]);

%% get brain mask
% look for qsm mask first
maskList = dir(fullfile(inputDir, '*mask-qsm*'));
% if no final mask then just look for normal mask
if isempty(maskList)
    maskList = dir(fullfile(inputDir, '*mask*'));
end

if ~isempty(maskFullName)
    % Option 1: mask file is provided
    maskFinal = load_nii_img_only(maskFullName) > 0;
    
elseif ~isempty(maskList) 
    % Option 2: input directory contains NIfTI file with name '*mask*'
    fprintf('A mask file is found in the input directory: %s\n',fullfile(inputDir, maskList(1).name));
    disp('Trying to load the file as signal mask');
    mask = load_nii_img_only(fullfile(inputDir, maskList(1).name)) > 0;
    
    % make sure the mask has the same dimension as other input data
    if ~isequal(size(mask),matrixSize)
        disp('The file does not have the same dimension as other images.')
        % display error message if nothing is found
        error('No mask file is loaded. Pleasee specific your mask file or put it in the input directory.');
    end
    
    disp('Mask file is loaded.');
    
else
    % display error message if nothing is found
    error('No mask file is loaded. Pleasee specific your mask file or put it in the input directory.');
    
end

%% make sure all variables are double
localField	= double(localField);
maskFinal   = double(maskFinal);
voxelSize   = double(voxelSize);
matrixSize  = double(matrixSize);
if exist('weights','var');  headerAndExtraData.weights  = double(weights).* maskFinal;	end
if exist('magn','var');     headerAndExtraData.magn     = double(magn);   	end

% create weighting map based on final mask
% for weighting map: higher SNR -> higher weights
% headerAndExtraData.weights = weights .* maskFinal;

headerAndExtraData.b0dir  	= B0_dir;
headerAndExtraData.b0     	= B0;
headerAndExtraData.te     	= TE;
headerAndExtraData.delta_TE	= delta_TE;
headerAndExtraData.CF    	= CF;

%% qsm
% core of QSM
[chi,mask_ref] = QSMMacro(localField,maskFinal,matrixSize,voxelSize,algorParam,headerAndExtraData);

% save results
fprintf('Saving susceptibility map...');

save_nii_quick(outputNiftiTemplate, chi, [outputDir filesep prefix 'QSM.nii.gz']);

if ~isempty(mask_ref)
    save_nii_quick(outputNiftiTemplate, mask_ref, [outputDir filesep prefix 'mask_reference_region.nii.gz']);
end

fprintf('Done!\n');

disp('Processing pipeline is completed!');

end

%% Validate nifti filenames with directory input
function CheckFileName(inputNiftiList)
% no. of files with particular name that has been read
numLocalFieldFile      = 0;

    % go through all files in the directory
    for klist = 1:length(inputNiftiList)
        if ContainName(lower(inputNiftiList(klist).name),'local-field')
            numLocalFieldFile = numLocalFieldFile + 1;
        end
    end
    
    % bring the error message if multiple files with the same string are
    % detected
    if numLocalFieldFile > 1
        error('Multiple files with name containing string ''local-field'' are detected. Please make sure only the magnitude data contains string ''local-field''');
    end
    
    % bring the error message if no file is detected
    if numLocalFieldFile == 0
        error('No file with name containing string ''local-field'' is detected. Please make sure the input directory contains a local field data');
    end

end
