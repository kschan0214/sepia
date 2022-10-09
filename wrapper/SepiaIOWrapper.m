%% [chi,localField,totalField,fieldmapSD]=SepiaIOWrapper(inputDir,outputDir,varargin)
%
% Input
% --------------
% input         :   input directory contains NIfTI files or structure containing filenames  
% output        :   output directory that stores the output (susceptibility map)
% maskFullName  :   mask filename
% algorParam    :   structure contains method and method specific parameters
%
% Output
% --------------
% totalField            : unwrapped field map (in Hz)
% fieldmapSD            : relative standard deviation of field map,
%                         esimated using Eq. 11 of Robinson et al. NMR Biomed 2017 (doi:10.1002/nbm.3601)
% localField            : local field (or tissue field) (in Hz)
% chi                   : quantitative susceptibility map (in ppm)
%
% Description: This is a wrapper of estimateTotalField.m which has the following objeectives:
%               (1) matches the input format of sepia.m
%               (2) save the results in NIfTI format
%
% Kwok-shing Chan @ DCCN
% k.chan@donders.ru.nl
% Date created: 14 September 2017
% Date modified: 26 August 2018
% Date modified: 29 March 2019
% Date modified: 27 Feb 2020 (v0.8.0)
% Date modified: 21 Jan 2020 (v0.8.1)
% Date modified: 6 May 2021 (v0.8.1.1)
% Date modified: 13 August 2021 (v1.0)
% Date modified: 12 September 2022 (v1.1)
%
function [chi,localField,totalField,fieldmapSD]=SepiaIOWrapper(input,output,maskFullName,algorParam)
%% add general Path and universal variables
sepia_addpath

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

outputFileList = construct_output_filename(outputDir, prefix);

%% Check and set default algorithm parameters
algorParam          = check_and_set_SEPIA_algorithm_default(algorParam);
% generl algorithm parameters
exclude_threshold	= algorParam.unwrap.excludeMaskThreshold;
exclude_method      = algorParam.unwrap.excludeMethod;
isSaveUnwrappedEcho = algorParam.unwrap.isSaveUnwrappedEcho;

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

%%%%%% Step 4: Check whether phase data contains DICOM values or wrapped phase value
% availableFileList	: structure contains all data filenames that are already available and validated
% outputFileList  	: structure contains default output filenames
availableFileList           = io_04_true_phase_value(availableFileList, outputFileList);

%%%%%% Step 5: in case user want to reverse the frequency shift direction
% availableFileList	: structure contains all data filenames that are already available and validated
% outputFileList  	: structure contains default output filenames
% algorParam        : structure contains all pipeline parameters
availableFileList           = io_05_reverse_phase(availableFileList, outputFileList, algorParam);

% display some header info
display_sepia_header_info_4wrapper;

%%%%%% Step 6: get signal mask
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
availableFileList           = io_06_get_signal_mask(maskFullName, inputDir, sepia_header, algorParam, availableFileList, outputFileList, outputNiftiTemplate);

%%%%%% Step 7: refine signal mask
% sepia_header          : sepia header
% algorParam            : structure contains all pipeline parameters
% availableFileList     : structure contains all data filenames that are already available and validated
% outputFileList        : structure contains default output filenames
% outputNiftiTemplate   : nifti header with empty 'img' field
availableFileList           = io_07_refine_signal_mask(sepia_header, algorParam, availableFileList, outputFileList, outputNiftiTemplate);

%%%%%% store some data to headerAndExtraData
% header
create_header_structure_4wrapper;

matrixSize  = double(sepia_header.matrixSize);
voxelSize   = double(sepia_header.voxelSize);
TE          = double(sepia_header.TE);

headerAndExtraData.availableFileList = availableFileList;
headerAndExtraData.outputDirectory   = outputDir; 

%% Main QSM processing - Step 1: total field and phase unwrap

%%%%%%%%%% Step 0: Eddy current correction for bipolar readout %%%%%%%%%%
% sepia_header          : sepia header
% algorParam            : structure contains all pipeline parameters
% availableFileList     : structure contains all data filenames that are already available and validated
% outputFileList        : structure contains default output filenames
% outputNiftiTemplate   : nifti header with empty 'img' field
availableFileList          = tf_00_bipolar_correction(sepia_header, algorParam, availableFileList, outputFileList, outputNiftiTemplate);

%%%%%%%%%% Step 1: Phase unwrapping and echo phase combination %%%%%%%%%%
fieldMap    = load_nii_img_only(availableFileList.phase);
mask        = load_nii_img_only(availableFileList.mask);

headerAndExtraData.availableFileList = availableFileList;

% core of temporo-spatial phase unwrapping
[totalField,fieldmapSD,fieldmapUnwrapAllEchoes,mask] = estimateTotalField(fieldMap,mask,matrixSize,voxelSize,algorParam,headerAndExtraData);

% save unwrapped phase if chosen
if ~isempty(fieldmapUnwrapAllEchoes) && isSaveUnwrappedEcho
    % save the output                           
    fprintf('Saving unwrapped echo phase...');
    save_nii_quick(outputNiftiTemplate,fieldmapUnwrapAllEchoes, outputFileList.unwrappedPhase);
    fprintf('Done!\n');
    
    availableFileList.unwrappedPhase = outputFileList.unwrappedPhase;
end
clear fieldmapUnwrapAllEchoes

% save the total fieldmap                       
fprintf('Saving unwrapped fieldmap...');
save_nii_quick(outputNiftiTemplate,totalField,  outputFileList.totalField);
fprintf('Done.\n');
availableFileList.totalField = outputFileList.totalField;

%%%%%%%%%% Step 2: exclude unreliable voxel, based on monoexponential decay model %%%%%%%%%%
% only work with multi-echo data
if length(TE) == 1 && ~isinf(exclude_threshold)
    warning('\nExcluding unreliable voxels can only work with multi-echo data.')
    disp('No voxels are excluded');
    exclude_threshold = inf;
end
    
    
if ~isinf(exclude_threshold)
    
    magn = double(load_nii_img_only(availableFileList.magnitude));
    
    % multi-echo data
    r2s                 = R2star_trapezoidal(magn,TE);
    relativeResidual    = ComputeResidualGivenR2sFieldmap(TE,r2s,totalField,magn.*exp(1i*fieldMap));
    maskReliable        = relativeResidual < exclude_threshold;
    % v1.1: 20220919
    relativeResidualWeights = relativeResidual;
    % clipping
    relativeResidualWeights(relativeResidualWeights>exclude_threshold) = exclude_threshold;
    % weightsRelativeResidual should be between [0,1]
    relativeResidualWeights = (exclude_threshold - relativeResidualWeights) ./ exclude_threshold;
    
    clear r2s magn 
    
    fprintf('Saving other output...');
    save_nii_quick(outputNiftiTemplate,maskReliable,   	outputFileList.maskReliable);
    save_nii_quick(outputNiftiTemplate,relativeResidual,outputFileList.relativeResidual);
    save_nii_quick(outputNiftiTemplate,relativeResidualWeights,outputFileList.relativeResidualWeights);
    fprintf('Done.\n');
    
    clear relativeResidual
    
    availableFileList.maskReliable              = outputFileList.maskReliable;
    availableFileList.relativeResidual          = outputFileList.relativeResidual;
    availableFileList.relativeResidualWeights   = outputFileList.relativeResidualWeights;
    
else
    % single-echo & no threshold
    maskReliable = ones(size(totalField),'like',totalField);
end


switch exclude_method
    % threshold fieldmapSD with the reliable voxel mask
    case 'Weighting map'
        fieldmapSD = fieldmapSD .* maskReliable;
    % threshold brain mask with the reliable voxel mask
    case 'Brain mask'
        mask = mask .* maskReliable;
        
end
save_nii_quick(outputNiftiTemplate,fieldmapSD,  outputFileList.fieldmapSD);
save_nii_quick(outputNiftiTemplate,mask,        outputFileList.maskLocalField); 

availableFileList.fieldmapSD        = outputFileList.fieldmapSD;
availableFileList.maskLocalField    = outputFileList.maskLocalField;

% create weighting map 
% for weighting map: higher SNR -> higher weighting
if ~isfield(availableFileList, 'weights')
    
    fprintf('Computing weighting map...');
    % weights = sepia_utils_compute_weights_v0p8(fieldmapSD,and(mask>0,maskReliable>0));
    weights = sepia_utils_compute_weights_v1(fieldmapSD,mask);
    weights = weights .* mask;

    % modulate weighting map by relativa residual
    if ~isinf(exclude_threshold) && strcmp(exclude_method,'Weighting map')
        weights = weights .* relativeResidualWeights;
    end

    save_nii_quick(outputNiftiTemplate,weights,	outputFileList.weights);
    availableFileList.weights = outputFileList.weights;
    
    fprintf('Done!\n');
else
    if ~isinf(exclude_threshold) % if user select thresholding their own weight
        % load user weights
        weights = double(load_nii_img_only(availableFileList.weights));
        
        % mask out unreliable voxel
        weights = weights .* mask;
%         weights = weights .* and(mask>0,maskReliable);

        % modulate weighting map by relativa residual
        if ~isinf(exclude_threshold) && strcmp(exclude_method,'Weighting map')
            weights = weights .* relativeResidualWeights;
        end
        
        % export modified weights and update filelist
        save_nii_quick(outputNiftiTemplate,weights,	outputFileList.weights);
        availableFileList.weights = outputFileList.weights;
    end
end

% clear variable that no longer be needed
clear fieldMap fieldmapSD weights maskReliable mask relativeResidualWeights

% update availableFileList
headerAndExtraData.availableFileList = availableFileList;

%% Background field removal
totalField   	= double(load_nii_img_only(availableFileList.totalField));
maskLocalfield	= double(load_nii_img_only(availableFileList.maskLocalField));

localField = BackgroundRemovalMacro(totalField,maskLocalfield,matrixSize,voxelSize,algorParam,headerAndExtraData);
clear totalField maskLocalfield % clear variables that no longer be needed

% generate new mask based on backgroudn field removal result
mask_QSM = localField ~=0;

fprintf('Saving local field map...');
save_nii_quick(outputNiftiTemplate,localField, outputFileList.localField);
fprintf('done!\n');
availableFileList.localField = outputFileList.localField;
clear localField

% save results
fprintf('Saving mask for chi mapping...');
save_nii_quick(outputNiftiTemplate,mask_QSM, outputFileList.maskQSM);
fprintf('done!\n');
availableFileList.maskQSM = outputFileList.maskQSM;
clear mask_QSM

% update availableFileList
headerAndExtraData.availableFileList = availableFileList;

%% QSM
% make sure all variables are double
localField   	= double(load_nii_img_only(availableFileList.localField));
mask_QSM        = double(load_nii_img_only(availableFileList.maskQSM));

% Apply final mask to weights
% headerAndExtraData.weights = headerAndExtraData.weights .* mask_QSM;

% core of QSM
[chi,mask_ref] = QSMMacro(localField,mask_QSM,matrixSize,voxelSize,algorParam,headerAndExtraData);
clear localField mask_QSM

% save results
fprintf('Saving susceptibility map...');
save_nii_quick(outputNiftiTemplate, chi, outputFileList.QSM);
clear chi

if ~isempty(mask_ref)
    save_nii_quick(outputNiftiTemplate, mask_ref, outputFileList.maskRef);
end
fprintf('done!\n');

disp('Processing pipeline is completed!');
          
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
    filePattern = {'ph','mag','weights','header'}; % don't change the order
    [isLoadSuccessful, inputNiftiList] = read_default_to_filelist(inputDir, filePattern);
    
    % If it doesn't work then check again for 'phase' instead of 'ph'
    if ~isLoadSuccessful
        filePattern = {'phase','mag','weights','header'}; % don't change the order
        [isLoadSuccessful, inputNiftiList] = read_default_to_filelist(inputDir, filePattern);
    end
    
    % If it doesn't work then check BIDS compatibility
    if ~isLoadSuccessful
        disp('Searching input directory based on BIDS...');
        inputNiftiList = read_bids_to_filelist(inputDir,fullfile(outputDir,prefix));
    end
    
end

end

%% I/O Step 2: validate input data
function availableFileList          = io_02_validate_nifti_input(inputFileList)

availableFileList = struct();

% load NIFTI header for validating input images
% 2.2 phase data 
if ~isempty(inputFileList(1).name)
    
    fprintf('Loading phase header files...')
    
    % get header info from NIFTI for validation
    phaseNIFTIHeader = load_untouch_header_only(inputFileList(1).name);
    
    % store the filename in the availableFileList structure
    availableFileList.phase = inputFileList(1).name;
    
    fprintf('Done.\n');
    
else
    error('Fail! \nPlease specify a single-echo(3D0/multi-echo(4D) phase data.');
end

% 2.2 magnitude data 
if ~isempty(inputFileList(2).name)
    
    fprintf('Loading magnitude header files...')
    
    % get header info from NIFTI for validation
    magnitudeNIFTIHeader = load_untouch_header_only(inputFileList(2).name);
    
    % store the filename in the availableFileList structure
    availableFileList.magnitude = inputFileList(2).name;
    
    fprintf('Done.\n');
    
else
    error('Fail! \nPlease specify a single-echo(3D-/multi-echo(4D) magnitude data.');
end

fprintf('Validating input phase and magnitude images...')
% make sure input phase and magnitude have the same dimension
matrixSize_magn     = magnitudeNIFTIHeader.dime.dim(2:5);
matrixSize_phase    = phaseNIFTIHeader.dime.dim(2:5);
% check matrix size between magnitude data and phase data
if ~isequal(matrixSize_magn,matrixSize_phase)
    error('Fail! \nInput phase and magnitude data do not have the same (3D/4D) matrix size. Please check the NIfTi files.');
else
    fprintf('Passed.\n');
end


% 2.3 Weights data 
if ~isempty(inputFileList(3).name)
    
    % get header info from NIFTI for validation
    weightsNIFTIHeader = load_untouch_header_only(inputFileList(3).name);
    
    availableFileList.weights = inputFileList(3).name;
   
else
    disp('No weighting map is loaded. Default QSM weighting method will be used for QSM.');
end

% check dimension of weights
if exist('weightsNIFTIHeader','var')
    if weightsNIFTIHeader.dime.dim(1) > 3
        error('Input weighting map is 4D. SEPIA accepts weighting map to be 3D only.');
    end
end

disp('Input files are valid.')

end

%% I/O Step 3: get nifti template for nifti output
function outputNiftiTemplate        = io_03_get_nifti_template(availableFileList)

outputNiftiTemplate     = load_untouch_nii(availableFileList.magnitude);
outputNiftiTemplate.img = [];

end

%% I/O Step 4: convert phase to radian unit if required
function availableFileList          = io_04_true_phase_value(availableFileList, outputFileList)

% load phase image to check if the 
phaseNIFTI = load_untouch_nii(availableFileList.phase);
phaseIMG    = load_nii_img_only(availableFileList.phase);

if abs(max(phaseIMG(:))-pi)>0.1 || abs(min(phaseIMG(:))-(-pi))>0.1 % allow small differences possibly due to data stype conversion or DICOM digitisation
% if abs(max(phaseNIFTI.img(:))-pi)>0.1 || abs(min(phaseNIFTI.img(:))-(-pi))>0.1 % allow small differences possibly due to data stype conversion or DICOM digitisation

    disp('Values of input phase map exceed the range of [-pi,pi]. DICOM value is assumed.')
    fprintf('Rescaling phase data from DICOM image value to wrapped radian unit...')
    phase = DICOM2Phase(phaseNIFTI);
    fprintf('Done.\n')

    fprintf('Saving phase images in unit of radian...');
    save_nii_quick(phaseNIFTI, phase, outputFileList.phaseRadian);
    fprintf('Done.\n')
    
    % update the phase data for QSM processing
    availableFileList.phase = outputFileList.phaseRadian;
    
end

end

%% I/O Step 5: reverse phase rotation if required
function availableFileList          = io_05_reverse_phase(availableFileList, outputFileList, algorParam)

if algorParam.general.isInvert
    
    phaseNIFTI = load_untouch_nii(availableFileList.phase);
    phaseNIFTI.img = -phaseNIFTI.img;
    
    disp('Phase data is reversed.')
    
    fprintf('Saving reversed phase images...');
    save_nii_quick(phaseNIFTI, phaseNIFTI.img, outputFileList.phaseReversed);
    fprintf('Done.\n')
    
    % update the phase data for QSM processing
    availableFileList.phase = outputFileList.phaseReversed;
    
end
    
end

%% I/O Step 6: loading signal mask
function availableFileList          = io_06_get_signal_mask(maskFullName, inputDir, sepia_header, algorParam, availableFileList, outputFileList, outputNiftiTemplate)

isBET               = algorParam.general.isBET;
fractional_threshold= algorParam.general.fractional_threshold;
gradient_threshold  = algorParam.general.gradient_threshold;

matrixSize  = sepia_header.matrixSize;
voxelSize   = sepia_header.voxelSize;

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
        availableFileList.mask = maskFullName;
        disp('Mask file is checked.');
    end
end

% if no mask is found then display the following message
if isempty(mask) && ~isBET
    disp('No mask data is loaded. Using FSL BET to obtain brain mask.');
end
    
% if BET is checked or no mask is found, run FSL's bet
if isempty(mask) || isBET
    
    magn = load_nii_img_only(availableFileList.magnitude);
    
    sepia_addpath('MEDI');
    
    disp('Performing FSL BET...');
    % Here uses MEDI toolboxes MEX implementation
    mask = BET(magn(:,:,:,1),matrixSize,voxelSize,fractional_threshold,gradient_threshold);
    disp('Signal mask is obtained.');
    
    fprintf('Saving signal mask...')
    save_nii_quick(outputNiftiTemplate,mask, outputFileList.maskBrain);
    fprintf('Done!\n');
    
    availableFileList.mask = outputFileList.maskBrain;
end

end

%% I/O Step 7: refine brain mask
function availableFileList          = io_07_refine_signal_mask(sepia_header, algorParam, availableFileList, outputFileList, outputNiftiTemplate)

TE          = sepia_header.TE;
voxelSize   = sepia_header.voxelSize;
isMultiEcho         = numel(TE)>1;
isRefineBrainMask   = algorParam.general.isRefineBrainMask;

if ~isMultiEcho
    isRefineBrainMask = 0;
    disp('Refine brain mask only works with multi-echo data');
end

if isRefineBrainMask
    disp('Refine brain using R2* info');
    magn        = double(load_nii_img_only(availableFileList.magnitude));
    mask        = double(load_nii_img_only(availableFileList.mask));
    r2s         = R2star_trapezoidal(magn, TE);
    mask_refine = refine_brain_mask_using_r2s(r2s,mask,voxelSize);

    % save the eddy current corrected output
    fprintf('Saving refined brain mask...');
    save_nii_quick(outputNiftiTemplate, mask_refine, outputFileList.maskRefine);
    fprintf('Done!\n');

    % update availableFileList
    availableFileList.mask = outputFileList.maskRefine;
end

end

%% TF Step 0: bipolar readout phase correction
function availableFileList          = tf_00_bipolar_correction(sepia_header, algorParam, availableFileList, outputFileList, outputNiftiTemplate)

isEddyCorrect   = algorParam.unwrap.isEddyCorrect;

TE              = sepia_header.TE;

if numel(TE) < 4 && isEddyCorrect
    
    warning('Bipolar readout correction requires data with at least 4 echoes.');
    disp('Bipolar readout correction is  not performed.');
    isEddyCorrect = false;
    
end

if isEddyCorrect
    
    % load data
    magn        = double(load_nii_img_only(availableFileList.magnitude));
    fieldMap    = double(load_nii_img_only(availableFileList.phase));
    mask        = double(load_nii_img_only(availableFileList.mask));

    % BipolarEddyCorrect requries complex-valued input
%     [imgCplx,bipolar_phase]	= BipolarEddyCorrect(magn.*exp(1i*fieldMap),mask,algorParam);
    [imgCplx,bipolar_phase]	= FastBipolarCorrect(magn.*exp(1i*fieldMap),mask);
    fieldMap            	= double(angle(imgCplx));
    
    % save the eddy current corrected output
    fprintf('Saving eddy current corrected phase data...');
    save_nii_quick(outputNiftiTemplate, fieldMap, outputFileList.phaseEddyCorr);
    save_nii_quick(outputNiftiTemplate, bipolar_phase, outputFileList.phase_bipolar);
    fprintf('Done!\n');
    
    % update availableFileList
    availableFileList.phase = outputFileList.phaseEddyCorr;
    
end

end