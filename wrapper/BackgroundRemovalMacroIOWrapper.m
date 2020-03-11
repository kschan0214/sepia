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
%
%
function [localField,maskFinal] = BackgroundRemovalMacroIOWrapper(input,output,maskFullName,algorParam)
%% add general Path
sepia_addpath;

%% define variables
prefix = 'sepia_';

isInputDir = true;
% make sure the input only load once (first one)
isTotalFieldLoad = false;
isFieldmapSDLoad = false;

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
    
    % take the total field map directory as reference input directory 
    [inputDir,~,~] = fileparts(inputNiftiList(1).name);
else
    % Option 2: input is a directory
    inputDir = input;
    inputNiftiList = dir([inputDir '/*.nii*']);
end

% Step 2: load data
if ~isempty(inputNiftiList)
    
    if ~isInputDir
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %%%%%%%%%%%%%%%% Pathway 1: Input are NIfTI files %%%%%%%%%%%%%%%%%
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        disp('NIfTI input is being used.');
        
                        %%%%%%%%%% Total field map %%%%%%%%%%
        if ~isempty(inputNiftiList(1).name)
            inputTotalFieldNifti = load_untouch_nii([inputNiftiList(1).name]);
            % load true value from NIfTI
            totalField = load_nii_img_only([inputNiftiList(1).name]);
%             totalField = double(inputTotalFieldNifti.img);
            isTotalFieldLoad = true;
            
            disp('Total field map is loaded.')
        else
            error('Please specify a 3D total field map.');
        end
        
                         %%%%%%%%%% Fieldmapsd data %%%%%%%%%%
        if ~isempty(inputNiftiList(3).name)
%             inputFieldMapSDNifti = load_untouch_nii([inputNiftiList(3).name]);
            % load true value from NIfTI
            fieldmapSD = load_nii_img_only([inputNiftiList(3).name]);
%             fieldmapSD = double(inputFieldMapSDNifti.img);
            isFieldmapSDLoad = true;
            
            disp('Noise SD data is loaded.')
        else
            disp('No field map standard deviation data is loaded.');
        end
        
                        %%%%%%%%%% sepia header %%%%%%%%%%
        if ~isempty(inputNiftiList(4).name)
            load([inputNiftiList(4).name]);
            disp('Header data is loaded.');
        else
            error('Please specify a header required by SEPIA.');
        end
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    else
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %%%%%%%%%% Pathway 2: Input is a directory with NIfTI %%%%%%%%%%
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % validate indput directory files
        disp('Directory input is being used.');
        fprintf('Validating filenames in the directory....');
        CheckFileName(inputNiftiList);
        fprintf('Filenames are valid.\n');
        
        % loop all NIfTI files in the directory for total field map and fieldmap SD (optional)
        for klist = 1:length(inputNiftiList)

                        %%%%%%%%%% Total field map %%%%%%%%%%
            if ContainName(inputNiftiList(klist).name,'total-field') && ~isTotalFieldLoad
                inputTotalFieldNifti = load_untouch_nii([inputDir filesep inputNiftiList(klist).name]);
                % load true value from NIfTI
                totalField = load_nii_img_only([inputDir filesep inputNiftiList(klist).name]);
%                 totalField = double(inputTotalFieldNifti.img);
                isTotalFieldLoad = true;

                disp('Total field map is loaded.')
            end

                         %%%%%%%%%% Fieldmapsd data %%%%%%%%%%
            if ContainName(inputNiftiList(klist).name,'noise-sd') && ~isFieldmapSDLoad
%                 inputFieldMapSDNifti = load_untouch_nii([inputDir filesep inputNiftiList(klist).name]);
                % load true value from NIfTI
                fieldmapSD = load_nii_img_only([inputDir filesep inputNiftiList(klist).name]);
%                 fieldmapSD = double(inputFieldMapSDNifti.img);
                isFieldmapSDLoad = true;

                disp('Noise SD data is loaded.')
            end
        end

        % if no files matched the name format then displays error message
        if ~isTotalFieldLoad
            error('No total field map is loaded. Please make sure the input directory contains files with name *total-field*');
        end
        if ~isFieldmapSDLoad
            disp('No noise standard deviation data is loaded.');
        end

        %%%%%%%%%% SEPIA header file %%%%%%%%%%
        if ~isempty(dir([inputDir '/*header*']))
            % load header
            headerList = dir([inputDir '/*header*']);
            load([inputDir filesep headerList(1).name]);

            disp('Header data is loaded.');

        else
            error('Please specify a header required by SEPIA.');
        end
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    end
    
    % store the header the NIfTI files, all following results will have
    % the same header
    outputNiftiTemplate = inputTotalFieldNifti;
    
else
    error('This standalone only reads NIfTI format input data (*.nii or *.nii.gz).');
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
maskList = dir([inputDir '/*mask*']);

if ~isempty(maskFullName)
    % Option 1: mask file is provided
    mask = load_nii_img_only(maskFullName) > 0;
    
elseif ~isempty(maskList) 
    % Option 2: input directory contains NIfTI file with name '*mask*'
    inputMaskNii = load_untouch_nii([inputDir filesep maskList(1).name]);
	mask = inputMaskNii.img > 0;
    
else
    % display error message if nothing is found
    error('No mask file is found. Please specify your mask file or put it in the input directory.');
    
end

%% make sure all variables are double
totalField  = double(totalField);
mask       	= double(mask);
voxelSize   = double(voxelSize);
matrixSize  = double(matrixSize);

headerAndExtraData = [];
if exist('fieldmapSD','var'); headerAndExtraData.N_std = double(fieldmapSD); end

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

%% Validate nifti filenames with directory input
function CheckFileName(inputNiftiList)
% no. of files with particular name that has been read
numTotalFieldFile      = 0;

    % go through all files in the directory
    for klist = 1:length(inputNiftiList)
        if ContainName(lower(inputNiftiList(klist).name),'total-field')
            numTotalFieldFile = numTotalFieldFile + 1;
        end
    end
    
    % bring the error message if multiple files with the same string are
    % detected
    if numTotalFieldFile > 1
        error('Multiple files with name containing string ''total-field'' are detected. Please make sure only the magnitude data contains string ''total-field''');
    end
    
    % bring the error message if no file is detected
    if numTotalFieldFile == 0
        error('No file with name containing string ''total-field'' is detected. Please make sure the input directory contains a total field map');
    end

end