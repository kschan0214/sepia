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
% Date modified: 16 November 2025 (v1.3): Add support to bids fMRI/fQSM format
%
function [chi,localField,totalField,fieldmapSD,chi_para,chi_dia] = SepiaIOWrapper(input,output,maskFullName,algorParam)
%% add general Path and universal variables
sepia_addpath

sepia_universal_variables;
suffix = get_nifti_extension_from_input(input);

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
fprintf('Output filename suffix : %s\n',suffix);

write_bids_dataset_description(outputDir);

%% Check and set default algorithm parameters
algorParam          = check_and_set_SEPIA_algorithm_default(algorParam);
% general algorithm parameters
exclude_threshold	= algorParam.unwrap.excludeMaskThreshold;
exclude_method      = algorParam.unwrap.excludeMethod;
isSaveUnwrappedEcho = algorParam.unwrap.isSaveUnwrappedEcho;
isSaveR2s           = algorParam.unwrap.isSaveR2s;
isMagnitudeCombine  = algorParam.unwrap.isMagnitudeCombine;

outputFileList = construct_output_filename(outputDir, prefix, algorParam, suffix);

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
[inputDir, inputNiftiCell]	= io_01_get_input_file_list(input, outputDir, prefix);

nVol = numel(inputNiftiCell);
for v = 1:nVol

if nVol > 1; fprintf('Processing #%i/%i volume\n',v,nVol); end

% load current inputFileList
inputFileList = inputNiftiCell{v}.inputNIFTIList;
if nVol > 1; outputFileList = construct_output_filename(outputDir, strcat(prefix,'vol-',num2str(v),'_'), algorParam); end

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
if numel(inputFileList) < 4 || isempty(inputFileList(4).name)
    error('Please specify a header required by SEPIA.');
else
    sepia_header = load([inputFileList(4).name]);
    disp('SEPIA header data is loaded.');
    % Validate header information
    sepia_header = validate_sepia_header_4wrapper(sepia_header, outputNiftiTemplate);
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
% for multi-volume, make sure the same mask is used throughout the process
if nVol > 1; maskFullName = availableFileList.mask; end

%%%%%% Step 8: Tensor-MPPCA denoising
% sepia_header          : sepia header
% algorParam            : structure contains all pipeline parameters
% availableFileList     : structure contains all data filenames that are already available and validated
% outputFileList        : structure contains default output filenames
% outputNiftiTemplate   : nifti header with empty 'img' field
availableFileList          = io_08_denoising(sepia_header, algorParam, availableFileList, outputFileList, outputNiftiTemplate);

%%%%%% Step 9: Upsampling
% sepia_header          : sepia header
% algorParam            : structure contains all pipeline parameters
% availableFileList     : structure contains all data filenames that are already available and validated
% outputFileList        : structure contains default output filenames
% outputNiftiTemplate   : nifti header with empty 'img' field
[availableFileList,sepia_header,outputNiftiTemplate] = io_09_upsampling(sepia_header, algorParam, availableFileList, outputFileList, outputNiftiTemplate);

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

% 20230124 v1.2.2: apply mask on derived map
totalField = totalField .* double(mask);

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
save_json_sidecar(outputFileList.totalField, struct( ...
    'Description', 'Unwrapped total field map estimated by temporo-spatial phase unwrapping.', ...
    'Units',       'Hz', ...
    'Sources',     {{get_relative_source_path(outputDir, availableFileList.phase)}}, ...
    'Parameters',  algorParam.unwrap));
fprintf('Done.\n');
availableFileList.totalField = outputFileList.totalField;

%%%%%%%%%% Step 2: exclude unreliable voxel, based on monoexponential decay model %%%%%%%%%%
% only work with multi-echo data
if isscalar(TE) && ~isinf(exclude_threshold)
    warning('\nExcluding unreliable voxels can only work with multi-echo data.')
    disp('No voxels are excluded');
    exclude_threshold = inf;
end
    
    
if ~isinf(exclude_threshold)

    magn = double(load_nii_img_only(availableFileList.magnitude));

    % multi-echo data
    % R2* map only needs to be computed once; reuse it if it is already
    % available, otherwise compute it and make it available for later use
    if isfield(availableFileList,'r2s') && exist(availableFileList.r2s,'file')
        disp('R2* map is already available. Loading it from disk...');
        r2s = double(load_nii_img_only(availableFileList.r2s));
    else
        r2s = R2star_trapezoidal(magn,TE);

        fprintf('Saving R2* map...');
        save_nii_quick(outputNiftiTemplate, r2s, outputFileList.r2s);
        save_json_sidecar(outputFileList.r2s, struct( ...
            'Description', 'R2* map estimated from multi-echo magnitude data (trapezoidal method), used to exclude unreliable voxels.', ...
            'Units',       '1/s', ...
            'Sources',     {{get_relative_source_path(outputDir, availableFileList.magnitude)}}));
        fprintf('Done!\n');

        availableFileList.r2s = outputFileList.r2s;
    end
    relativeResidual    = ComputeResidualGivenR2sFieldmap(TE,r2s,totalField,magn.*exp(1i*fieldMap));
    maskReliable        = relativeResidual < exclude_threshold;
    % v1.1: 20220919
    relativeResidualWeights = relativeResidual;
    % clipping
    relativeResidualWeights(relativeResidualWeights>exclude_threshold) = exclude_threshold;
    % weightsRelativeResidual should be between [0,1]
    relativeResidualWeights = (exclude_threshold - relativeResidualWeights) ./ exclude_threshold;

    % Save r2s and optimal combined magnitude
    if isSaveR2s
        fprintf('Saving R2star map...');
        save_nii_quick(outputNiftiTemplate,r2s,   	                outputFileList.r2s);
        save_json_sidecar(outputFileList.r2s, struct( ...
            'Description', 'R2* map estimated from multi-echo magnitude data (trapezoidal method), used to exclude unreliable voxels.', ...
            'Units',       '1/s', ...
            'Sources',     {{get_relative_source_path(outputDir, availableFileList.magnitude)}}));
    end
    if isMagnitudeCombine
        fprintf('Combining multi-echo data optimally...');
        optimalCombinedMagnitude = ComputeOptimalCombinedMagnitude(TE,r2s,magn);
        save_nii_quick(outputNiftiTemplate,optimalCombinedMagnitude,outputFileList.optimalCombinedMagnitude);
        save_json_sidecar(outputFileList.optimalCombinedMagnitude, struct( ...
            'Description', 'Optimally combined multi-echo magnitude image.', ...
            'Units',       'arbitrary', ...
            'Sources',     {{get_relative_source_path(outputDir, availableFileList.magnitude)}}));

        clear optimalCombinedMagnitude
    end

    clear r2s magn

    fprintf('Saving other output...');
    save_nii_quick(outputNiftiTemplate,maskReliable,   	outputFileList.maskReliable);
    save_json_sidecar(outputFileList.maskReliable, struct( ...
        'Description', 'Reliable voxel mask based on a monoexponential decay model of the multi-echo magnitude data.', ...
        'Sources',     {{get_relative_source_path(outputDir, availableFileList.magnitude)}}));
    save_nii_quick(outputNiftiTemplate,relativeResidual,outputFileList.relativeResidual);
    save_json_sidecar(outputFileList.relativeResidual, struct( ...
        'Description', 'Relative residual of the monoexponential decay model fit used to identify unreliable voxels.', ...
        'Units',       'arbitrary', ...
        'Sources',     {{get_relative_source_path(outputDir, availableFileList.magnitude)}}));
    save_nii_quick(outputNiftiTemplate,relativeResidualWeights,outputFileList.relativeResidualWeights);
    save_json_sidecar(outputFileList.relativeResidualWeights, struct( ...
        'Description', 'Weighting map derived from the relative residual of the monoexponential decay model fit.', ...
        'Units',       'arbitrary', ...
        'Sources',     {{get_relative_source_path(outputDir, availableFileList.magnitude)}}));
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
save_json_sidecar(outputFileList.fieldmapSD, struct( ...
    'Description', 'Noise standard deviation of the total field map.', ...
    'Units',       'arbitrary', ...
    'Sources',     {{get_relative_source_path(outputDir, availableFileList.phase)}}, ...
    'Parameters',  algorParam.unwrap));
save_nii_quick(outputNiftiTemplate,mask,        outputFileList.maskLocalField);
save_json_sidecar(outputFileList.maskLocalField, struct( ...
    'Description', 'Signal mask used for background field removal, after excluding unreliable voxels.', ...
    'Sources',     {{get_relative_source_path(outputDir, availableFileList.phase)}}));

availableFileList.fieldmapSD        = outputFileList.fieldmapSD;
availableFileList.maskLocalField    = outputFileList.maskLocalField;

% create weighting map 
% for weighting map: higher SNR -> higher weighting
if ~isfield(availableFileList, 'weights')
    
    fprintf('Computing weighting map...');
    % weights = sepia_utils_compute_weights_v0p8(fieldmapSD,and(mask>0,maskReliable>0));
    weights = sepia_utils_compute_weights_v1(fieldmapSD,mask);
    weights = weights .* mask;

    % modulate weighting map by relative residual
    if ~isinf(exclude_threshold) && strcmp(exclude_method,'Weighting map')
        weights = weights .* relativeResidualWeights;
    end

    save_nii_quick(outputNiftiTemplate,weights,	outputFileList.weights);
    save_json_sidecar(outputFileList.weights, struct( ...
        'Description', 'Weighting map based on field map noise standard deviation, for use in QSM dipole inversion.', ...
        'Units',       'arbitrary', ...
        'Sources',     {{get_relative_source_path(outputDir, availableFileList.phase)}}));
    availableFileList.weights = outputFileList.weights;

    fprintf('Done!\n');
else
    if ~isinf(exclude_threshold) % if user select thresholding their own weight
        % load user weights
        weights = double(load_nii_img_only(availableFileList.weights));
        
        % mask out unreliable voxel
        weights = weights .* mask;
%         weights = weights .* and(mask>0,maskReliable);

        % modulate weighting map by relative residual
        if ~isinf(exclude_threshold) && strcmp(exclude_method,'Weighting map')
            weights = weights .* relativeResidualWeights;
        end
        
        % export modified weights and update filelist
        save_nii_quick(outputNiftiTemplate,weights,	outputFileList.weights);
        save_json_sidecar(outputFileList.weights, struct( ...
            'Description', 'User-supplied weighting map, after masking out unreliable voxels.', ...
            'Units',       'arbitrary', ...
            'Sources',     {{get_relative_source_path(outputDir, availableFileList.phase)}}));
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

% generate new mask based on background field removal result
% mask_QSM = localField ~=0;
% 20230124 v1.2.2: make sure no holes inide ROIs
mask_QSM = imfill(localField ~= 0, 'holes');

fprintf('Saving local field map...');
save_nii_quick(outputNiftiTemplate,localField, outputFileList.localField);
save_json_sidecar(outputFileList.localField, struct( ...
    'Description', 'Local (tissue) field map after background field removal.', ...
    'Units',       'Hz', ...
    'Sources',     {{get_relative_source_path(outputDir, availableFileList.totalField)}}, ...
    'Parameters',  algorParam.bfr));
fprintf('done!\n');
availableFileList.localField = outputFileList.localField;
clear localField

% save results
fprintf('Saving mask for chi mapping...');
save_nii_quick(outputNiftiTemplate,mask_QSM, outputFileList.maskQSM);
save_json_sidecar(outputFileList.maskQSM, struct( ...
    'Description', 'Signal mask for QSM dipole inversion, derived from the background field removal result.', ...
    'Sources',     {{get_relative_source_path(outputDir, availableFileList.totalField)}}));
fprintf('done!\n');
availableFileList.maskQSM = outputFileList.maskQSM;
clear mask_QSM

% update availableFileList
headerAndExtraData.availableFileList = availableFileList;

%% QSM
% make sure all variables are double
localField   	= double(load_nii_img_only(availableFileList.localField));

% Apply final mask to weights
% headerAndExtraData.weights = headerAndExtraData.weights .* mask_QSM;

% Two-pass masking
if not(strcmpi(algorParam.qsm.isTwoPass, 'None'))
    % backup algorParam, availableFileList, and outputFileList
    algorParamTwoPass = algorParam;
    algorParamTwoPass.msk.refineMethod = algorParam.qsm.isTwoPass;
    algorParamTwoPass.msk.threshold    = algorParam.qsm.twopass_lambda;
    availableFileListTwoPass = availableFileList;
    availableFileListTwoPass.mask      = availableFileList.maskQSM;
    outputFileListTwoPass = outputFileList;
    outputFileListTwoPass.maskReliable = outputFileList.maskQSM2pass;
    availableFileList = MaskRefinementIOWrapper(sepia_header, ...
                                                algorParamTwoPass, ...
                                                availableFileListTwoPass, ...
                                                outputFileListTwoPass, ...
                                                outputNiftiTemplate);
    mask_QSM_pass_2 = double(load_nii_img_only(availableFileList.maskReliable));
    mask_QSM_pass_1 = double(load_nii_img_only(availableFileList.maskQSM));
    mask_QSM{1} = mask_QSM_pass_1;
    mask_QSM{2} = mask_QSM_pass_2;
else
    mask_QSM        = double(load_nii_img_only(availableFileList.maskQSM));
end

% core of QSM
% 20260822 KC: expanded for chi-sep type output
[chi,mask_ref,chi_para,chi_dia] = QSMMacro(localField,mask_QSM,matrixSize,voxelSize,algorParam,headerAndExtraData);
clear localField mask_QSM

% save results
%20250903 PSF: expanded for two-pass masking output
fprintf('Saving susceptibility map...');
if iscell(chi)
    save_nii_quick(outputNiftiTemplate, chi{1}, outputFileList.QSM);
    save_json_sidecar(outputFileList.QSM, struct( ...
        'Description', 'Quantitative susceptibility map from dipole field inversion.', ...
        'Units',       'ppm', ...
        'Sources',     {{get_relative_source_path(outputDir, availableFileList.localField)}}, ...
        'Parameters',  algorParam.qsm));
    save_nii_quick(outputNiftiTemplate, chi{2}, outputFileList.QSMpass1);
    save_json_sidecar(outputFileList.QSMpass1, struct( ...
        'Description', 'Quantitative susceptibility map from dipole field inversion with first mask.', ...
        'Units',       'ppm', ...
        'Sources',     {{get_relative_source_path(outputDir, availableFileList.localField)}}, ...
        'Parameters',  algorParam.qsm));
    save_nii_quick(outputNiftiTemplate, chi{3}, outputFileList.QSMpass2);
    save_json_sidecar(outputFileList.QSMpass2, struct( ...
        'Description', 'Quantitative susceptibility map from dipole field inversion with second mask.', ...
        'Units',       'ppm', ...
        'Sources',     {{get_relative_source_path(outputDir, availableFileList.localField)}}, ...
        'Parameters',  algorParam.qsm));
else
    save_nii_quick(outputNiftiTemplate, chi, outputFileList.QSM);
    save_json_sidecar(outputFileList.QSM, struct( ...
        'Description', 'Quantitative susceptibility map from dipole field inversion.', ...
        'Units',       'ppm', ...
        'Sources',     {{get_relative_source_path(outputDir, availableFileList.localField)}}, ...
        'Parameters',  algorParam.qsm));
end
% 20260822 KC: expanded for chi-sep type output
if ~isempty(chi_para)
    if iscell(chi_para)
        chi_para = chi_para{1};
    end
    save_nii_quick(outputNiftiTemplate, chi_para, outputFileList.QSMpara);
    save_json_sidecar(outputFileList.QSMpara, struct( ...
        'Description', 'Paramagnetic susceptibility component map from chi-separation.', ...
        'Units',       'ppm', ...
        'Sources',     {{get_relative_source_path(outputDir, availableFileList.localField)}}, ...
        'Parameters',  algorParam.qsm));
end
if ~isempty(chi_dia)
    if iscell(chi_dia)
        chi_dia = chi_dia{1};
    end
    save_nii_quick(outputNiftiTemplate, chi_dia, outputFileList.QSMdia);
    save_json_sidecar(outputFileList.QSMdia, struct( ...
        'Description', 'Diamagnetic susceptibility component map from chi-separation.', ...
        'Units',       'ppm', ...
        'Sources',     {{get_relative_source_path(outputDir, availableFileList.localField)}}, ...
        'Parameters',  algorParam.qsm));
end
% clear chi chi_para chi_dia

if ~isempty(mask_ref)
    save_nii_quick(outputNiftiTemplate, mask_ref, outputFileList.maskRef);
    save_json_sidecar(outputFileList.maskRef, struct( ...
        'Description', 'Reference region mask used for susceptibility referencing.', ...
        'Sources',     {{get_relative_source_path(outputDir, availableFileList.localField)}}));
end
fprintf('done!\n');

disp('Processing pipeline is completed!');

end

% concatenate multiple volumes if needed
if nVol > 1
    tmp = [];
    for v = 1:nVol  
        tmp = cat(4,tmp,load_nii_img_only(fullfile(outputDir, strcat(prefix,'vol-',num2str(v),'_Chimap.nii.gz'))));
    end
    save_nii_img_only(fullfile(outputDir, strcat(prefix,'vol-',num2str(v),'_Chimap.nii.gz')), fullfile(outputDir, strcat(prefix,'Chimap.nii.gz')),tmp);
end
end

%% I/O Step 1: get input file list
function [inputDir, inputNiftiCell] = io_01_get_input_file_list(input,outputDir,prefix)

if isstruct(input)
    
    % Option 1: input are files
    inputNiftiCell{1}.inputNIFTIList  = input;
    
    % take the phase data directory as reference input directory 
    [inputDir,~,~] = fileparts(inputNiftiCell{1}.inputNIFTIList(1).name);
    
else
    
    % Option 2: input is a directory
    inputDir = input; 
    
    % First check with SEPIA default naming structure
    disp('Searching input directory based on SEPIA default naming structure...');
    filePattern = {'ph','mag','weights','header'}; % don't change the order
    [isLoadSuccessful, inputNIFTIList] = read_default_to_filelist(inputDir, filePattern);
    
    % If it doesn't work then check again for 'phase' instead of 'ph'
    if ~isLoadSuccessful
        filePattern = {'phase','mag','weights','header'}; % don't change the order
        [isLoadSuccessful, inputNIFTIList] = read_default_to_filelist(inputDir, filePattern);
    end
    inputNiftiCell{1}.inputNIFTIList = inputNIFTIList;
    
    % If it doesn't work then check BIDS compatibility
    if ~isLoadSuccessful
        disp('Searching input directory based on BIDS...');
        inputNiftiCell = read_bids_to_filelist(inputDir,fullfile(outputDir,prefix));
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
if numel(inputFileList) < 3 || isempty(inputFileList(3).name)
    disp('No weighting map is loaded. Default QSM weighting method will be used for QSM.');
else
    % get header info from NIFTI for validation
    weightsNIFTIHeader = load_untouch_header_only(inputFileList(3).name);
    
    availableFileList.weights = inputFileList(3).name;    
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

% PSF20251110: Separate wrapper to ensure SepiaIOWrapper and 
% UnwrapPhaseMacroIOWrapper use the same masking structure and backend
availableFileList = MaskWrapper(maskFullName, inputDir, sepia_header, algorParam, availableFileList, outputFileList, outputNiftiTemplate);

end

%% I/O Step 7: refine brain mask
function availableFileList          = io_07_refine_signal_mask(sepia_header, algorParam, availableFileList, outputFileList, outputNiftiTemplate)

% PSF20251110: Separate wrapper for mask refinement, also used for two-pass
% masking, and to unify the sub-modules of UnwrapPhaseMacroIOWrapper and
% SepiaIOWrapper
if algorParam.general.isRefineBrainMask
    algorParam.msk.refineMethod = 'r2s-refine'; % hard-coded for the time being, need to confim with PF
    availableFileList = MaskRefinementIOWrapper(sepia_header, algorParam, availableFileList, outputFileList, outputNiftiTemplate);
end

end

%% I/O Step 8: image denoising
function availableFileList          = io_08_denoising(sepia_header, algorParam, availableFileList, outputFileList, outputNiftiTemplate)

sepia_universal_variables;

if algorParam.general.isDenoise

    % check if tensor MPPCA code exist, if not then download form GitHub
    setup_tMPPCA_toolbox();

    kernel = ceil(algorParam.general.denoiseKernel ./ sepia_header.voxelSize);
    if any(kernel<3)
        warning('Denoising kernel size too small. [3x3x3] voxels kernel will be used');
    end
    kernel = max(kernel,[3,3,3]); % minimum window are 3 voxels

    disp('Tensor-MP-PCA denoising in progress (can take some time)...')
    magn        = double(load_nii_img_only(availableFileList.magnitude));
    phase       = double(load_nii_img_only(availableFileList.phase));
    mask        = double(load_nii_img_only(availableFileList.mask)) >0;

    % create complex-valued image
    img         = magn .* exp(1i*phase);
    
    tic
    [img_denoise,sigma,P,snrgain] = denoise_recursive_tensor(img,kernel,'mask',mask);
    toc
    
    save_nii_quick(outputNiftiTemplate, abs(img_denoise),   outputFileList.magDenoise);
    save_nii_quick(outputNiftiTemplate, angle(img_denoise), outputFileList.phaseDenoise);
    save_nii_quick(outputNiftiTemplate, sigma,              outputFileList.sigma);
    save_nii_quick(outputNiftiTemplate, P,                  outputFileList.P);
    save_nii_quick(outputNiftiTemplate, snrgain,            outputFileList.snrgain);

    % update availableFileList
    availableFileList.magnitude = outputFileList.magDenoise;
    availableFileList.phase     = outputFileList.phaseDenoise;

    disp('Done!');
end
end

%% I/O Step 9: image upsampling
function [availableFileList,sepia_header,outputNiftiTemplate] = io_09_upsampling(sepia_header, algorParam, availableFileList, outputFileList, outputNiftiTemplate)

sepia_universal_variables;


if algorParam.general.isUpsample

    disp('Upsampling in progress...')
    magn        = double(load_nii_img_only(availableFileList.magnitude));
    phase       = double(load_nii_img_only(availableFileList.phase));
    mask        = double(load_nii_img_only(availableFileList.mask)) >0;

    % create complex-valued image
    img         = magn .* exp(1i*phase);

    scaleFactor = sepia_header.voxelSize./algorParam.general.target_resolution ;
    if any(scaleFactor < 1)
        warning('You are downsampling the data. This step will be skipped.');
        return
    end
    matrixSize_upSample = round(scaleFactor .* sepia_header.matrixSize);

    % Assuming `data` is [X Y Z Echo] complex GRE data
    numEchoes   = size(img, 4);
    img_us = zeros([matrixSize_upSample numEchoes], 'like', img);
    
    for e = 1:numEchoes
        img_us(:,:,:,e) = fft_upsample_complex(img(:,:,:,e), matrixSize_upSample);
    end

    upsampledMask = imresize3(double(mask), matrixSize_upSample, 'cubic')> 0.2;

    % update sepia header
    sepia_header.matrixSize = matrixSize_upSample;
    sepia_header.voxelSize  = ones(size(sepia_header.voxelSize))*algorParam.general.target_resolution;

    % save output
    outputNiftiTemplate.hdr.dime.pixdim(2:4) = sepia_header.voxelSize;

    save_nii_quick(outputNiftiTemplate, abs(img_us),    outputFileList.magUpsample);
    save_nii_quick(outputNiftiTemplate, angle(img_us),  outputFileList.phaseUpsample);
    save_nii_quick(outputNiftiTemplate, upsampledMask,  outputFileList.maskUpsample);
    TE = sepia_header.TE; B0 = sepia_header.B0;
    save(outputFileList.sepiaHeaderUpsample,'TE','B0')

    % update availableFileList
    availableFileList.magnitude     = outputFileList.magUpsample;
    availableFileList.phase         = outputFileList.phaseUpsample;
    availableFileList.mask          = outputFileList.maskUpsample;
    availableFileList.sepiaheader   = outputFileList.sepiaHeaderUpsample;

    disp('Done!');
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


