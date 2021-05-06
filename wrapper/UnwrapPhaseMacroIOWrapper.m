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
% Date modified: 21 Jan 2020 (v0.8.1)
% Date modified: 6 May 2021 (v0.8.1.1)
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

% display output info
fprintf('Output directory       : %s\n',outputDir);
fprintf('Output filename prefix : %s\n',prefix);

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

%% Setting up Input
disp('---------');
disp('Load data');
disp('---------');

%%%%%% Step 1: get all required filenames
if isstruct(input)
    
    % Option 1: input are files
    inputNiftiList = input;
    
    % take the phase data directory as reference input directory 
    [inputDir,~,~] = fileparts(inputNiftiList(1).name);
    
else
    % Option 2: input is a directory
    disp('Searching input directory...');
    inputDir = input;

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

%%%%%% Step 2: load data
% 2.1 phase data 
if ~isempty(inputNiftiList(1).name)
    
    % load true value from NIfTI
    fieldMap    = load_nii_img_only([inputNiftiList(1).name]); 
    isPhaseLoad = true;
    disp('Phase data is loaded.')
    
else
    error('Please specify a single-echo/multi-echo phase data.');
end

% 2.2 magnitude data 
if ~isempty(inputNiftiList(2).name)
    
    % load nifti structure for output template
    inputMagnNifti  = load_untouch_nii([inputNiftiList(2).name]);
    % load true value from NIfTI
    magn            = load_nii_img_only([inputNiftiList(2).name]);
    isMagnLoad      = true;
    disp('Magnitude data is loaded.');
    
else
    error('Please specify a single-echo/multi-echo magnitude data.');
end

% 2.3 SEPIA header 
if ~isempty(inputNiftiList(4).name)
    load([inputNiftiList(4).name]);
    disp('Header data is loaded.');
else
    error('Please specify a header required by SEPIA.');
end

% store the header of the NIfTI files, all following results will have
% the same header
outputNiftiTemplate = inputMagnNifti;

%%%%%% Step 3: validate input
% 3.1 Validate header information
validate_sepia_header_4wrapper;

% 3.2 Validate NIfTI input
disp('Validating input NIfTI files...')
% make sure the dimension of input data is consistent
check_input_dimension(magn,fieldMap,matrixSize,TE);
disp('Input NIfTI files are valid.')

%%%%%% Step 4: Basic correction
% 4.1: check whether phase data contains DICOM values or wrapped phase value
if abs(max(fieldMap(:))-pi)>1e-4 || abs(min(fieldMap(:))-(-pi))>1e-4 % allow small differences possibly due to data stype conversion

    disp('Values of input phase map exceed the range of [-pi,pi]. DICOM value is assumed.')
    fprintf('Rescaling phase data from DICOM image value to wrapped radian unit...')
    inputPhaseNifti = load_untouch_nii([inputNiftiList(1).name]);
    fieldMap        = DICOM2Phase(inputPhaseNifti);
    fprintf('Done.\n')

    fprintf('Saving phase images in unit of radian...');
    save_nii_quick(inputPhaseNifti,fieldMap, [outputDir filesep prefix 'phase.nii.gz']);
    fprintf('Done.\n')
    
    clearvars inputPhaseNifti
    
end

% 4.2: in case user want to reverse the frequency shift direction
if isInvert
    disp('Phase data is reversed.')
    fieldMap = -fieldMap;
end

% display some header info
display_sepia_header_info_4wrapper;

% ensure variables are double
magn        = double(magn);
fieldMap    = double(fieldMap);
TE          = double(TE);
matrixSize  = double(matrixSize);  
voxelSize   = double(voxelSize);  

%%%%%% Step 5: store some data to headerAndExtraData
% header
create_header_structure_4wrapper;

headerAndExtraData.magn     = magn;
headerAndExtraData.phase    = fieldMap;

clearvars inputMagnNifti

%% get brain mask
disp('-----------');
disp('Signal mask');
disp('-----------');

mask        = [];
maskList    = dir(fullfile(inputDir,'*mask*nii*'));

% Scenario: No specified mask file + No check BET + there is a file called mask in the input directory
if isempty(maskFullName) && ~isempty(maskList) && ~isBET
    
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
        mask = [];
    else
        disp('Mask file is loaded.');
    end
end

% if no mask is found then display the following message
if isempty(mask) && ~isBET
    disp('No mask data is loaded. Using FSL BET to obtain brain mask.');
end
    
% if BET is checked or no mask is found, run FSL's bet
if isempty(mask) || isBET
    
    sepia_addpath('MEDI');
    
    disp('Performing FSL BET...');
    % Here uses MEDI toolboxes MEX implementation
    mask = BET(magn(:,:,:,1),matrixSize,voxelSize,fractional_threshold,gradient_threshold);
    disp('Signal mask is obtained.');
    
    fprintf('Saving signal mask...')
    save_nii_quick(outputNiftiTemplate,mask, [outputDir filesep prefix 'mask.nii.gz']);
    fprintf('done!\n');

end

% put the mask to headerAndExtraData
headerAndExtraData.mask = mask;

%% Step 0: Eddy current correction for bipolar readout
if numel(TE) < 4 && isEddyCorrect
    
    warning('Bipolar readout correction requires data with at least 4 echoes.');
    disp('Bipolar readout correction is not performed.');
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
    
    % update phase in headerAndExtraData
    headerAndExtraData.phase    = fieldMap;
    
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

%% Validate input data dimension
function check_input_dimension(magn,phase,matrixSize,TE)

matrixSize_magn     = size(magn);
matrixSize_phase    = size(phase);

if ndims(magn) == 3
    matrixSize_magn(4) = 1;
end
if ndims(phase) == 3
    matrixSize_phase(4) = 1;
end

% check matrix size between magnitude data and phase data
if ~isequal(matrixSize_magn,matrixSize_phase)
    error('Input phase and magnitude data do not have the same (3D/4D) matrix size. Please check the NIfTi files.');
end

% check matrix size between NIfTi input and SEPIA header
if ~isequal(matrixSize_magn(1:3),matrixSize)
    erro('Input NIfTI data and SEPIA header do not have the same matrix size. Please check these files and/or remove the ''matrixSize'' variable from the SEPIA header.')
end

% check echo dimension 
if matrixSize_magn(4) ~= length(TE)
    error('Input NIfTI data and SEPIA header do not have the same number of echoes.  Please check these files.');
end

end