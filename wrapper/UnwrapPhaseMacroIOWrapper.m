%% [totalField,fieldmapSD,mask]=UnwrapPhaseMacroIOWrapper(input,output,maskFullName,algorParam)
% 
% Input
% --------------
% input         : input directory contains NIfTI (*mag* and *ph*) files or structure containing filenames
% output        : output directory that stores the output (Unwrapped total field and fieldMapSD)
% maskFullName  : mask filename
% algorParam    : structure contains method and method specific parameters
%
% Output
% --------------
% totalField    : total field (or background + tissue fields) (in Hz)
% fieldmapSD    : noise standard deviation in the phase
% mask          : signal (brain) mask
%
% Description: This is a wrapper of estimateTotalField.m which has the following objeectives:
%               (1) matches the input format of SEPIA.m
%               (2) save the results in NIfTI format
%
% Kwok-shing Chan @ DCCN
% k.chan@donders.ru.nl
% Date created: 17 April 2018
% Date modified: 16 September 2018
% Date modified: 29 March 2019
% Date modified: 9 March 2020
% Date modified: 20 Jan 2020 (v0.8.1)
%
%
function [totalField,fieldmapSD,mask]=UnwrapPhaseMacroIOWrapper(input,output,maskFullName,algorParam)
%% add general Path
sepia_addpath;

sepia_universal_variables;

%% define variables
prefix = 'sepia_';
% make sure the input only load once (first one)
isMagnLoad  = false;
isPhaseLoad = false;

%% Check output directory exist or not
output_index    = strfind(output, filesep);
outputDir       = output(1:output_index(end));
if ~isempty(output(output_index(end)+1:end))
    prefix = [output(output_index(end)+1:end) '_'];
end

if exist(outputDir,'dir') ~= 7
    % if not then create the directory
    mkdir(outputDir);
end

%% Check and set default algorithm parameters
algorParam          = check_and_set_SEPIA_algorithm_default(algorParam);
isInvert            = algorParam.general.isInvert;
isBET               = algorParam.general.isBET ;
fractional_threshold= algorParam.general.fractional_threshold;
gradient_threshold  = algorParam.general.gradient_threshold;
exclude_threshold	= algorParam.unwrap.excludeMaskThreshold;
exclude_method      = algorParam.unwrap.excludeMethod;
isEddyCorrect      	= algorParam.unwrap.isEddyCorrect;
isSaveUnwrappedEcho = algorParam.unwrap.isSaveUnwrappedEcho;

%% Read input
disp('Reading data...');

% Step 1: check input for nifti files first
if isstruct(input)
    % Option 1: input are files
    inputNiftiList = input;
    
    % take the phase data directory as reference input directory 
    [inputDir,~,~] = fileparts(inputNiftiList(1).name);
else
    % Option 2: input is a directory
    inputDir = input;
%     inputNiftiList = dir([inputDir '/*.nii*']);

    % check and get filenames
    inputNiftiList = struct();
    filePattern = {'ph','mag','','header'}; % don't change the order
    for  k = 1:length(filePattern)
        if k ~= 3
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
                error(['No file with name containing string ''' filePattern{k} ''' is detected.']);
            else % multiple files -> fatal error
                error(['Multiple files with name containing string ''' filePattern{k} ''' are detected. Make sure the input directory should contain only one file with string ''' filePattern{k} '''.']);
            end
        else
            inputNiftiList(k).name = '';    % third field is empty in this application
        end
    end
    
end

% Step 2: load data
%%%%%%%%%% 2.1 phase data %%%%%%%%%%
if ~isempty(inputNiftiList(1).name)
    
    % load true value from NIfTI
    fieldMap = load_nii_img_only([inputNiftiList(1).name]);
    
    disp('Phase data is loaded.')
    
    % check whether phase data contains DICOM values or wrapped
    % phase value
    if max(fieldMap(:))>4 || min(fieldMap(:))<-4
        disp('Values of input phase map exceed range of [-pi,pi]. DICOM value is assumed.')
        fprintf('Converting phase data from DICOM image value to radian unit...')
        inputPhaseNifti = load_untouch_nii([inputNiftiList(1).name]);
        fieldMap        = DICOM2Phase(inputPhaseNifti);
        fprintf('Done.\n')

        fprintf('Saving phase images in unit of radian...');
        save_nii_quick(inputPhaseNifti,fieldMap, [outputDir filesep prefix 'phase.nii.gz']);
        fprintf('Done.\n')
        
        clearvars inputPhaseNifti
    end
    
    isPhaseLoad = true;
    
else
    error('Please specify a single-echo/multi-echo phase data.');
end

%%%%%%%%%% magnitude data %%%%%%%%%%
if ~isempty(inputNiftiList(2).name)
    inputMagnNifti = load_untouch_nii([inputNiftiList(2).name]);
    % load true value from NIfTI
    magn = load_nii_img_only([inputNiftiList(2).name]);

    isMagnLoad = true;
    
    disp('Magnitude data is loaded.');
else
    error('Please specify a single-echo/multi-echo magnitude data.');
end

%%%%%%%%%% SEPIA header %%%%%%%%%%
if ~isempty(inputNiftiList(4).name)
    load([inputNiftiList(4).name]);
    disp('Header data is loaded.');
else
    error('Please specify a header required by SEPIA.');
end

% store the header of the NIfTI files, all following results will have
% the same header
outputNiftiTemplate = inputMagnNifti;
clearvars inputMagnNifti

% In case some parameters are missing in the header file
if ~exist('matrixSize','var')
    matrixSize = size(magn);
    matrixSize = matrixSize(1:3);
end
if ~exist('voxelSize','var')
    voxelSize = outputNiftiTemplate.hdr.dime.pixdim(2:4);
end

% validate loaded input
CheckInputDimension(magn,fieldMap,matrixSize,TE);

% in case user want to correct the frequency shift direction
if isInvert
    fieldMap = -fieldMap;
end

% make sure the L2 norm of B0 direction = 1
B0_dir = B0_dir ./ norm(B0_dir);

% display some header info
disp('----------------------');
disp('Basic Data information');
disp('----------------------');
disp(['Voxel size(x,y,z mm^3) =  ' num2str(voxelSize(1)) 'x' num2str(voxelSize(2)) 'x' num2str(voxelSize(3))]);
disp(['matrix size(x,y,z)     =  ' num2str(matrixSize(1)) 'x' num2str(matrixSize(2)) 'x' num2str(matrixSize(3))]);
disp(['B0 direction(x,y,z)    =  ' num2str(B0_dir(:)')]);
disp(['Field strength(T)      =  ' num2str(B0)]);
disp(['Number of echoes       = ' num2str(length(TE))]);

%% get brain mask
mask = [];
maskList = dir(fullfile(inputDir,'*mask*nii*'));
if ~isempty(maskFullName)
    % Option 1: mask file is provided
    mask = load_nii_img_only(maskFullName) > 0;
    
elseif ~isempty(maskList) 
    % Option 2: input directory contains NIfTI file with name '*mask*'
    fprintf('A mask file is found in the input directory: %s\n',fullfile(inputDir, maskList(1).name));
    disp('Trying to load the file as signal mask');
    mask = load_nii_img_only(fullfile(inputDir, maskList(1).name)) > 0;
    
    % make sure the mask has the same dimension as other input data
    if ~isequal(size(mask),matrixSize)
        disp('The file does not have the same dimension as other images.')
        mask = [];
    end
    
    disp('Mask file is loaded.');
    
end

% if no mask is found then display the following message
if isempty(mask) && ~isBET
    disp('No mask data is loaded, using FSL BET to obtain brain mask.');
end
    
% if BET is checked or no mask is found, run FSL's bet
if isempty(mask) || isBET
    sepia_addpath('MEDI');
    
    fprintf('Performing FSL BET...');
    % this is the BET function provided with MEDI toolbox
    mask = BET(magn(:,:,:,1),matrixSize,voxelSize,fractional_threshold,gradient_threshold);
    fprintf('done!\n');
    
    fprintf('Saving brain mask...')
    save_nii_quick(outputNiftiTemplate,mask, [outputDir filesep prefix 'mask.nii.gz']);
    fprintf('done!\n');

end

%% ensure all variable are double
% make sure the following variables are row vectors
matrixSize  = double(matrixSize(:).');
voxelSize   = double(voxelSize(:).');
fieldMap  	= double(fieldMap);
mask      	= double(mask);

headerAndExtraData.magn     = double(magn);
headerAndExtraData.b0       = B0;
headerAndExtraData.te       = TE;
headerAndExtraData.delta_TE = delta_TE;

%% Step 0: Eddy current correction for bipolar readout
if numel(TE) < 4 && isEddyCorrect
    warning('Bipolar readout correction requires data with at least 4 echoes.');
    disp('No bipolar readout correction is performed.');
    isEddyCorrect = false;
end
    
if isEddyCorrect
    
    % BipolarEddyCorrect requries complex-valued input
    imgCplx     = BipolarEddyCorrect(magn.*exp(1i*fieldMap),mask,algorParam);
    fieldMap    = angle(imgCplx);
    
    fprintf('Saving eddy current corrected phase data...');
    % save the eddy current corrected output
    save_nii_quick(outputNiftiTemplate,fieldMap,    [outputDir filesep prefix 'phase_eddy-correct.nii.gz']);
    fprintf('Done!\n');
    
    clear imgCplx
end

%% total field and phase unwrap
% Step 1: Phase unwrapping and echo phase combination
% core of temporo-spatial phase unwrapping
[totalField,fieldmapSD,fieldmapUnwrapAllEchoes] = estimateTotalField(fieldMap,mask,matrixSize,voxelSize,algorParam,headerAndExtraData);

% save unwrapped phase if chosen
if ~isempty(fieldmapUnwrapAllEchoes) && isSaveUnwrappedEcho
    % save the output                           
    fprintf('Saving unwrapped echo phase...');
    save_nii_quick(outputNiftiTemplate,fieldmapUnwrapAllEchoes,[outputDir filesep prefix 'unwrapped-phase.nii.gz']);
    clear fieldmapUnwrapAllEchoes
    fprintf('Done!\n');
end
            
%% Step 2: exclude unreliable voxel, based on monoexponential decay model with
% single freuqnecy shift
fprintf('Computing weighting map...');
% only work with multi-echo data
if length(TE) == 1 && ~isinf(exclude_threshold)
    fprintf('\n');
    warning('Excluding unreliable voxels can only work with multi-echo data.')
    disp('No voxels are excluded');
    exclude_threshold = inf;
end

if ~isinf(exclude_threshold)
    % multi-echo data
    r2s = R2star_trapezoidal(magn,TE);
    relativeResidual = ComputeResidualGivenR2sFieldmap(TE,r2s,totalField,magn.*exp(1i*fieldMap));
    maskReliable = relativeResidual < exclude_threshold;
else
    % single-echo & no threshold
    maskReliable = ones(size(totalField),'like',fieldmapSD);
end

switch exclude_method
    % threshold fieldmapSD with the reliable voxel mask
    case 'Weighting map'
        fieldmapSD = fieldmapSD .* maskReliable;
    % threshold brain mask with the reliable voxel mask
    case 'Brain mask'
        mask = mask .* maskReliable;
        
end

% 20180815: test with creating weights using relativeResidual
% weightResidual = 1-(relativeResidual./exclude_threshold);
% weightResidual(weightResidual>1) = 1;
% weightResidual(weightResidual<0) = 0;
% computing weights
weights                 = 1./fieldmapSD;
weights(isinf(weights))	= 0;
weights(isnan(weights))	= 0;
weights                 = weights./max(weights(and(mask>0,maskReliable>0)));
weights                 = weights .* and(mask>0,maskReliable);
fprintf('Done!\n');
             
% save the output                           
fprintf('Saving unwrapped field map...');

save_nii_quick(outputNiftiTemplate,totalField,  [outputDir filesep prefix 'total-field.nii.gz']);
save_nii_quick(outputNiftiTemplate,fieldmapSD,  [outputDir filesep prefix 'noise-sd.nii.gz']);
save_nii_quick(outputNiftiTemplate,weights,  	[outputDir filesep prefix 'weights.nii.gz']);

% additional output
if ~isinf(exclude_threshold)
    save_nii_quick(outputNiftiTemplate,maskReliable,   	[outputDir filesep prefix 'mask-reliable.nii.gz']);
    save_nii_quick(outputNiftiTemplate,relativeResidual,[outputDir filesep prefix 'relative-residual.nii.gz']);
    if strcmpi(exclude_method,'Brain mask')
       save_nii_quick(outputNiftiTemplate,mask,         [outputDir filesep prefix 'mask-local_field.nii.gz']); 
    end
end
fprintf('Done!\n');

disp('Processing pipeline is completed!');

end

%% Validate input data consistence
function CheckInputDimension(magn,phase,matrixSize,TE)
% check spatial dimension
for ndim = 1:3
    if size(magn,ndim) ~= matrixSize(ndim)
        error(['The ' num2str(ndim) 'st/nd/rd dimension of the magnitude data does not match with the header information']);
    end
    if size(phase,ndim) ~= matrixSize(ndim)
        error(['The ' num2str(ndim) 'st/nd/rd dimension of the phase data does not match with the header information']);
    end
end

% check echo time
if size(magn,4) ~= length(TE)
    error('The no. of echo in the magnitude data does not match with the header infomation');
end
if size(phase,4) ~= length(TE)
    error('The no. of echo in the phase data does not match with the header infomation');
end
if size(phase,4) ~= size(magn,4)
    error('The no. of echo in the magnitude data does not match with that of the phase data');
end

end