%% [chi,localField,totalField,fieldmapSD]=SepiaIOWrapper(inputDir,outputDir,varargin)
%
% Input
% --------------
% inputDir              : input directory contains DICOM or NIfTI (*magn* and *phase*) files 
% outputDir             : output directory that stores the output (Unwrapped total field and fieldMapSD)
% varargin ('Name','Value' pair)
% ---------
% 'FSLBet'              : boolen brain extraction using FSL's bet (defaut: false)
% 'mask'                : mask file (in NIfTI format) full name 
% 'unwrap'              : phase unwrapping method (default 'Laplacian')
% 'subsampling'         : subsamling factor for graph-cut unwrapping (default:1) (unsupport yet)
% 'exclude_threshold'   : threshold to exclude high SD voxel in analysis, [0,1] (default: 1 (unthrehold))
% 'eddy'                : boolean eddy current correction for bipolar readout data (default: false)
% 'BFR'                 : background field removal method (default: 'LBV')
% 'refine'              : remove the remainin B1 field inhomogeneity using 4th order polynomial fitting (defualt: true) 
% 'BFR_tol'             : tolerance of LBV or PDF (overloaded) (default: 0.01)
% 'depth'               : depth of LBV (default: 5)
% 'peel'                : layers to be peeled of LBV (default: 2)
% 'BFR_iteration'       : no. of iterations of PDF or iHARPERELLA (overloaded) (default: 50)
% 'BFR_padsize'         : zeropad size of PDF (default: 40)
% 'BFR_radius'          : radius of spherical mean kernel of SHARP, RESHARP, VSHARP and VSHARPSTI
%                         (overloaded) (default: 4)
% 'BFR_alpha'           : regularisation parameter of RESHARP (default: 0.01)
% 'BFR_threshold'       : threshold of SHARP (default: 0.03)
% 'QSM'                 : QSM method (default: 'TKD')
% 'QSM_threshold'       : threshold of TKD (defualt: 0.15) 
% 'QSM_lambda'          : regularisation parameter of TKD, CFL2, iLSQR, FANSI and MEDI (overloaded) (default: 0.13)
% 'QSM_optimise'        : boolean automatically estimate regularisation parameter based on L-curve approach of CFL2 and iLSQR (overloaded) (default: false)
% 'QSM_tol'             : tolerance of iLSQR and FANSI (overloaded)(default: 1e-3)
% 'QSM_iteration'       : no. of iterations of iLSQR, STISuiteiLQR and FANSI (overloaded) (default: 50)
% 'QSM_tol1'            : step 1 tolerance of STISuiteiLQR (default: 0.01)
% 'QSM_tol2'            : step 2 tolerance of STISuiteiLQR (default: 0.001)
% 'QSM_padsize'         : pad size of STISuiteiLQR (default: [4,4,4])
% 'QSM_mu'              : regularisation parameter of data consistency of FANSI (default: 5e-5)
% 'QSM_zeropad'         : size of zero-padding of MEDI (default: 0)
% 'QSM_wData'           : weighting of data of MEDI (default: 1)
% 'QSM_wGradient'       : weighting of gradient regularisation of MEDI (default: 1)
% 'QSM_radius'          : radius for the spherical mean value operator of MEDI (default: 5)
% 'QSM_isSMV'           : boolean using spherical mean value operator of MEDI (default: false)
% 'QSM_merit'           : boolean model error reduction through iterative tuning of MEDI (default: false)
% 'QSM_isLambdaCSF'     : boolen automatic zero reference (MEDI+0) (required CSF mask) (default: false)
% 'QSM_lambdaCSF'       : regularisation parameter of CSF reference of MEDI (default: 100)
% varargin ('flag')
% 'linear'              : linear solver for FANSI (default)
% 'non-linear'          : non-linear solver for FANSI
% 'TV'                  : Total variation constraint for FANSI (default)
% 'TGV'                 : total generalisaed variation for FANSI 
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
%
%
function [chi,localField,totalField,fieldmapSD]=SepiaIOWrapper(input,output,maskFullName,algorParam)
%% add general Path
sepia_addpath

sepia_universal_variables;

%% define variables
prefix = 'sepia_';   
isInputDir = true;
% make sure the input only load once (first one)
isMagnLoad      = false;
isPhaseLoad     = false;
isWeightLoad    = false;

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

%% Check and set default algorithm parameters
algorParam          = check_and_set_SEPIA_algorithm_default(algorParam);
% generl algorithm parameters
isInvert            = algorParam.general.isInvert;
isBET               = algorParam.general.isBET;
fractional_threshold = algorParam.general.fractional_threshold;
gradient_threshold   = algorParam.general.gradient_threshold;
% phase unwrap algorithm parameters
isEddyCorrect      	= algorParam.unwrap.isEddyCorrect;
exclude_threshold	= algorParam.unwrap.excludeMaskThreshold;
exclude_method      = algorParam.unwrap.excludeMethod;
isSaveUnwrappedEcho = algorParam.unwrap.isSaveUnwrappedEcho;


%% Read input
disp('-----');
disp('Input');
disp('-----');
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
    inputNiftiList = dir([inputDir '/*.nii*']);
end

% Step 2: load data
if ~isempty(inputNiftiList)
    
    if ~isInputDir
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %%%%%%%%%%%%%%%% Pathway 1: Input are NIfTI files %%%%%%%%%%%%%%%%%
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        disp('NIfTI input is being used.');
        
                        %%%%%%%%%% magnitude data %%%%%%%%%%
        if ~isempty(inputNiftiList(2).name)
            inputMagnNifti = load_untouch_nii([inputNiftiList(2).name]);
            % load true value from NIfTI
            magn = load_nii_img_only([inputNiftiList(2).name]);
%             magn = double(inputMagnNifti.img);
            isMagnLoad = true;
            disp('Magnitude data is loaded.');
        else
            error('Please specify a single-echo/multi-echo magnitude data.');
        end
        
                         %%%%%%%%%% phase data %%%%%%%%%%
        if ~isempty(inputNiftiList(1).name)
            inputPhaseNifti = load_untouch_nii([inputNiftiList(1).name]);
            % load true value from NIfTI
            fieldMap = load_nii_img_only([inputNiftiList(1).name]);
%             fieldMap = double(inputPhaseNifti.img);
            % check whether phase data contains DICOM values or wrapped
            % phase value
            if max(fieldMap(:))>4 || min(fieldMap(:))<-4
                disp('Converting phase data from DICOM image value to radian unit...')
                fieldMap = DICOM2Phase(inputPhaseNifti);

                disp('Saving phase images in unit of radian...');
                save_nii_quick(inputPhaseNifti,fieldMap, [outputDir filesep prefix 'phase.nii.gz']);

            end
            isPhaseLoad = true;
            disp('Phase data is loaded.')
        else
            error('Please specify a single-echo/multi-echo phase data.');
        end

                         %%%%%%%%%% Weights data %%%%%%%%%%
        if ~isempty(inputNiftiList(3).name)
%             inputWeightNifti = load_untouch_nii([inputNiftiList(3).name]);
            % load true value from NIfTI
            weights = load_nii_img_only([inputNiftiList(3).name]);
%             weights = double(inputWeightNifti.img);
            % check whether phase data contains DICOM values or wrapped
            if size(weights,4) > 1
                error('QSM weighting image has to be 3D.');
            end
            isWeightLoad = true;
        else
            disp('Default QSM weighting method will be used for QSM.');
        end
        
                        %%%%%%%%%% SEPIA header %%%%%%%%%%
        if ~isempty(inputNiftiList(4).name)
            load([inputNiftiList(4).name]);
            disp('Header data is loaded.');
        else
            error('Please specify a header required by SEPIA.');
        end
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    else
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %%%%%%%%%%% Pathway 2: Input is a directory with NIfTI %%%%%%%%%%%%
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % validate indput directory files
        disp('Directory input is being used.');
        fprintf('Validating filenames in the directory....');
        CheckFileName(inputNiftiList);
        fprintf('Filenames are valid.\n');
        
        % loop all NIfTI files in the directory for magnitude and phase files
        for klist = 1:length(inputNiftiList)
            
                            %%%%%%%%%% magnitude data %%%%%%%%%%
            if ContainName(lower(inputNiftiList(klist).name),'mag') && ~isMagnLoad
                inputMagnNifti = load_untouch_nii([inputDir filesep inputNiftiList(klist).name]);
                % load true value from NIfTI
                magn = load_nii_img_only([inputDir filesep inputNiftiList(klist).name]);
%                 magn = double(inputMagnNifti.img);
                isMagnLoad = true;
                disp('Magnitude data is loaded.');
            end

                            %%%%%%%%%% phase data %%%%%%%%%%
            if ContainName(lower(inputNiftiList(klist).name),'ph') && ~isPhaseLoad
                inputPhaseNifti = load_untouch_nii([inputDir filesep inputNiftiList(klist).name]);
                % load true value from NIfTI
                fieldMap = load_nii_img_only([inputDir filesep inputNiftiList(klist).name]);
%                 fieldMap = double(inputPhaseNifti.img);

                % if input fieldmap is directly converted from nifti converter
                % then converts the fieldmap to unit of radian and save to output dir
                if max(fieldMap(:))>4 || min(fieldMap(:))<-4
                    disp('Converting phase data from DICOM image value to radian unit ...')
                    fieldMap = DICOM2Phase(inputPhaseNifti);

                    fprintf('Saving phase images in unit of radian...');
                    save_nii_quick(inputPhaseNifti,fieldMap, [outputDir filesep prefix 'phase.nii.gz']);
                    fprintf('Done!\n');

                end
                isPhaseLoad = true;
                disp('Phase data is loaded.')
            end
        end

                        %%%%%%%%%% sepia header file %%%%%%%%%%
        if ~isempty(dir([inputDir '/*header*']))
            % load header
            headerList = dir([inputDir '/*header*']);
            load([inputDir filesep headerList(1).name]);

            disp('Header data is loaded.');

        else
            error('Please specify a header file for Sepia.');

        end
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    end
    
    % store the header of the NIfTI files, all following results will have
    % the same header
    outputNiftiTemplate = inputMagnNifti;
    
end

% validate loaded input
CheckInputDimension(magn,fieldMap,matrixSize,TE);

% in case user want to reverse the frequency shift direction
if isInvert
    fieldMap = -fieldMap;
end

% make sure the L2 norm of B0 direction = 1
B0_dir = B0_dir ./ norm(B0_dir);

% display some header info
disp('----------------------');
disp('Basic data information');
disp('----------------------');
disp(['Voxel size(x,y,z mm^3)   =  ' num2str(voxelSize(1)) 'x' num2str(voxelSize(2)) 'x' num2str(voxelSize(3))]);
disp(['matrix size(x,y,z)       =  ' num2str(matrixSize(1)) 'x' num2str(matrixSize(2)) 'x' num2str(matrixSize(3))]);
disp(['B0 direction(x,y,z)      =  ' num2str(B0_dir(:)')]);
disp(['Field strength(T)        =  ' num2str(B0)]);
disp(['Number of echoes         =  ' num2str(length(TE))]);


% ensure variables are double
magn        = double(magn);
fieldMap    = double(fieldMap);
TE          = double(TE);
matrixSize  = double(matrixSize(:).');  % row vectors
voxelSize   = double(voxelSize(:).');   % row vectors

%% get brain mask
mask = [];
maskList = dir([inputDir '/*mask*nii*']);
if ~isempty(maskFullName)
    % Option 1: mask file is provided
    mask = load_nii_img_only(maskFullName) > 0;
    
elseif ~isempty(maskList) 
    % Option 2: input directory contains NIfTI file with name '*mask*'
    inputMaskNii = load_untouch_nii([inputDir filesep maskList(1).name]);
	mask = inputMaskNii.img > 0;
    
end

% if no mask is found then display the following message
if isempty(mask) && ~isBET
    disp('No mask data is loaded, using FSL BET to obtain brain mask.');
end
    
% if BET is checked or no mask is found, run FSL's bet
if isempty(mask) || isBET
    sepia_addpath('MEDI');
    
    fprintf('Performing FSL BET...');
    % Here uses MEDI toolboxes MEX implementation
    mask = BET(magn(:,:,:,1),matrixSize,voxelSize,fractional_threshold,gradient_threshold);
    fprintf('done!\n');
    
    fprintf('Saving brain mask...')
    save_nii_quick(outputNiftiTemplate,mask, [outputDir filesep prefix 'mask.nii.gz']);
    fprintf('done!\n');

end

%% store some data to headerAndExtraData

% header
headerAndExtraData.b0       = B0;
headerAndExtraData.b0dir  	= B0_dir;
headerAndExtraData.te       = TE;
headerAndExtraData.delta_TE = delta_TE;
headerAndExtraData.CF    	= CF;

headerAndExtraData.magn     = double(magn);
headerAndExtraData.phase    = fieldMap;
%% total field and phase unwrap

%%%%%%%%%% Step 0: Eddy current correction for bipolar readout %%%%%%%%%%
if isEddyCorrect

    % BipolarEddyCorrect requries complex-valued input
    imgCplx     = BipolarEddyCorrect(magn.*exp(1i*fieldMap),mask,algorParam);
    fieldMap    = double(angle(imgCplx));
    
    % save the eddy current corrected output
    fprintf('Saving eddy current corrected phase data...');
    save_nii_quick(outputNiftiTemplate,fieldMap,    [outputDir filesep prefix 'phase_eddy-correct.nii.gz']);
    fprintf('Done!\n');
    
    % update phase in headerAndExtraData
    headerAndExtraData.phase    = fieldMap;
    
    clear imgCplx
    
end

%%%%%%%%%% Step 1: Phase unwrapping and echo phase combination %%%%%%%%%%
% core of temporo-spatial phase unwrapping
[totalField,fieldmapSD,fieldmapUnwrapAllEchoes] = estimateTotalField(fieldMap,mask,matrixSize,voxelSize,algorParam,headerAndExtraData);

% save unwrapped phase if chosen
if ~isempty(fieldmapUnwrapAllEchoes) && isSaveUnwrappedEcho
    % save the output                           
    fprintf('Saving unwrapped echo phase...');
    save_nii_quick(outputNiftiTemplate,fieldmapUnwrapAllEchoes,[outputDir filesep prefix 'unwrapped-phase.nii.gz']);
    fprintf('Done!\n');
    
    clear fieldmapUnwrapAllEchoes
end

%%%%%%%%%% Step 2: exclude unreliable voxel, based on monoexponential decay model %%%%%%%%%%
fprintf('Computing weighting map...');
% only work with multi-echo data
if length(TE) == 1
    warning('\nExcluding unreliable voxels can only work with multi-echo data.')
    disp('No voxels are excluded');
    exclude_threshold = inf;
end
    
    
if ~isinf(exclude_threshold)
    % multi-echo data
    r2s                 = R2star_trapezoidal(magn,TE);
    relativeResidual    = ComputeResidualGivenR2sFieldmap(TE,r2s,totalField,magn.*exp(1i*fieldMap));
    maskReliable        = relativeResidual < exclude_threshold;
    
    clear r2s
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
fprintf('Done!\n');

% create weighting map 
% for weighting map: higher SNR -> higher weighting
if ~isWeightLoad
    weights                 = 1./fieldmapSD;
    weights(isinf(weights)) = 0;
    weights(isnan(weights)) = 0;
    weights                 = weights./max(weights(and(mask>0,maskReliable>0)));
end
weights = weights .* and(mask>0,maskReliable);
             
% save the output                           
fprintf('Saving unwrapped field map...');

save_nii_quick(outputNiftiTemplate,totalField,  [outputDir filesep prefix 'total-field.nii.gz']);
save_nii_quick(outputNiftiTemplate,fieldmapSD,  [outputDir filesep prefix 'noise-sd.nii.gz']);
if ~isWeightLoad
    save_nii_quick(outputNiftiTemplate,weights,	[outputDir filesep prefix 'weights.nii.gz']);
end

if ~isinf(exclude_threshold)
    
    save_nii_quick(outputNiftiTemplate,maskReliable,   	[outputDir filesep prefix 'mask-reliable.nii.gz']);
    save_nii_quick(outputNiftiTemplate,relativeResidual,[outputDir filesep prefix 'relative-residual.nii.gz']);
    
    if strcmpi(exclude_method,'Brain mask')
       save_nii_quick(outputNiftiTemplate,mask,         [outputDir filesep prefix 'mask-local_field.nii.gz']); 
    end
    
    clear relativeResidual
end
fprintf('Done!\n');

% store some output to extradata
headerAndExtraData.N_std    = double(fieldmapSD);
headerAndExtraData.weights  = double(weights);

% clear variable that no longer be needed
clear fieldMap magn fieldmapSD weights maskReliable

%% Background field removal
% make sure all variables are double
totalField  = double(totalField);
mask       	= double(mask);

localField = BackgroundRemovalMacro(totalField,mask,matrixSize,voxelSize,algorParam,headerAndExtraData);
  
% generate new mask based on backgroudn field removal result
maskFinal = double(localField ~=0);
  
% save results
fprintf('Saving local field map...');

save_nii_quick(outputNiftiTemplate,localField,  [outputDir filesep prefix 'local-field.nii.gz']);
save_nii_quick(outputNiftiTemplate,maskFinal,	[outputDir filesep prefix 'mask-qsm.nii.gz']);
fprintf('done!\n');

% clear variables that no longer be needed
clear totalField mask

%% QSM
% make sure all variables are double
localField	= double(localField);
maskFinal   = double(maskFinal);

% Apply final mask to weights
headerAndExtraData.weights = headerAndExtraData.weights .* maskFinal;

% core of QSM
[chi,mask_ref] = QSMMacro(localField,maskFinal,matrixSize,voxelSize,algorParam,headerAndExtraData);
  
% save results
fprintf('Saving susceptibility map...');
save_nii_quick(outputNiftiTemplate, chi, [outputDir filesep prefix 'QSM.nii.gz']);
if ~isempty(mask_ref)
    save_nii_quick(outputNiftiTemplate, mask_ref, [outputDir filesep prefix 'mask_reference_region.nii.gz']);
end
fprintf('done!\n');

disp('Processing pipeline is completed!');
          
end

%% Validate nifti filenames with directory input
function CheckFileName(inputNiftiList)
% no. of files with particular name that has been read
numMagFile      = 0;
numPhaseFile    = 0;

    % go through all files in the directory
    for klist = 1:length(inputNiftiList)
        if ContainName(lower(inputNiftiList(klist).name),'mag')
            numMagFile = numMagFile + 1;
        end
        if ContainName(lower(inputNiftiList(klist).name),'ph')
            numPhaseFile = numPhaseFile + 1;
        end
    end
    
    % bring the error message if multiple files with the same string are
    % detected
    if numMagFile > 1
        error('Multiple files with name containing string ''mag'' are detected. Please make sure only the magnitude data contains string ''mag''');
    end
    if numPhaseFile > 1
        error('Multiple files with name containing string ''ph'' are detected. Please make sure only the phase data contains string ''ph''');
    end
    
    % bring the error message if no file is detected
    if numMagFile == 0
        error('No file with name containing string ''mag'' is detected. Please make sure the input directory contains a magnitude data');
    end
    if numPhaseFile == 0
        error('No file with name containing string ''ph'' is detected. Please make sure the input directory contains a phase data');
    end

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