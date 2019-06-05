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
%
%
function [chi,localField,totalField,fieldmapSD]=SepiaIOWrapper(input,output,maskFullName,algorParam)
%% add general Path
sepia_addpath

%% define variables
prefix = 'sepia_';   
gyro = 42.57747892;
isInputDir = true;
% make sure the input only load once (first one)
isMagnLoad      = false;
isPhaseLoad     = false;
isWeightLoad    = false;
maskCSF = [];

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
algorParam = CheckAndSetDefault(algorParam);
% generl algorithm parameters
isInvert            = algorParam.general.isInvert;
isBET               = algorParam.general.isBET ;
isGPU               = algorParam.general.isGPU;
% phase unwrap algorithm parameters
isEddyCorrect      	= algorParam.unwrap.isEddyCorrect;
phaseCombMethod    	= algorParam.unwrap.echoCombMethod;
unwrap              = algorParam.unwrap.unwrapMethod;
subsampling         = algorParam.unwrap.subsampling;
exclude_threshold	= algorParam.unwrap.excludeMaskThreshold;
% background field removal algorithm parameters
BFR                 = algorParam.bfr.method;
refine              = algorParam.bfr.refine;
BFR_tol             = algorParam.bfr.tol;
BFR_depth           = algorParam.bfr.depth;
BFR_peel            = algorParam.bfr.peel;
BFR_iteration       = algorParam.bfr.iteration;
BFR_padSize         = algorParam.bfr.padSize;
BFR_radius          = algorParam.bfr.radius;
BFR_alpha           = algorParam.bfr.alpha;
BFR_threshold       = algorParam.bfr.threshold;
BFR_erode_radius    = algorParam.bfr.erode_radius;
% QSM algorithm parameters
QSM_method          = algorParam.qsm.method; 
QSM_threshold       = algorParam.qsm.threshold; 
QSM_lambda          = algorParam.qsm.lambda; 
QSM_optimise        = algorParam.qsm.optimise; 
QSM_tol             = algorParam.qsm.tol;   
QSM_maxiter         = algorParam.qsm.maxiter;
QSM_tol1            = algorParam.qsm.tol1;  
QSM_tol2            = algorParam.qsm.tol2; 
QSM_padsize         = algorParam.qsm.padsize;
QSM_mu1             = algorParam.qsm.mu1; 
QSM_mu2             = algorParam.qsm.mu2;  
QSM_solver          = algorParam.qsm.solver;  
QSM_constraint      = algorParam.qsm.constraint; 
QSM_radius          = algorParam.qsm.radius;
QSM_zeropad         = algorParam.qsm.zeropad;   
QSM_wData           = algorParam.qsm.wData; 
QSM_wGradient       = algorParam.qsm.wGradient;
QSM_isLambdaCSF     = algorParam.qsm.isLambdaCSF;
QSM_lambdaCSF       = algorParam.qsm.lambdaCSF; 
QSM_isSMV           = algorParam.qsm.isSMV;
QSM_merit           = algorParam.qsm.merit;  

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
        
                        %%%%%%%%%% qsm hub header %%%%%%%%%%
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
disp('Basic data information');
disp(['Voxel size(x,y,z mm^3) =  ' num2str(voxelSize(1)) 'x' num2str(voxelSize(2)) 'x' num2str(voxelSize(3))]);
disp(['matrix size(x,y,z) =  ' num2str(matrixSize(1)) 'x' num2str(matrixSize(2)) 'x' num2str(matrixSize(3))]);
disp(['B0 direction(x,y,z) =  ' num2str(B0_dir(:)')]);
disp(['Field strength(T) =  ' num2str(B0)]);
disp(['Number of echoes = ' num2str(length(TE))]);

% make sure the following variables are row vectors
matrixSize = matrixSize(:).';
voxelSize = voxelSize(:).';

% ensure variables are double
magn        = double(magn);
fieldMap    = double(fieldMap);
TE          = double(TE);
matrixSize  = double(matrixSize);
voxelSize   = double(voxelSize);

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
    sepia_addpath('bet');
    
    fprintf('Performing FSL BET...');
    % Here uses MEDI toolboxes MEX implementation
    mask = BET(magn(:,:,:,1),matrixSize,voxelSize);
    fprintf('done!\n');
    
    fprintf('Saving brain mask...')
    save_nii_quick(outputNiftiTemplate,mask, [outputDir filesep prefix 'mask.nii.gz']);
    fprintf('done!\n');

end

% double data type
mask = double(mask);

%% total field and phase unwrap

% Step 0: Eddy current correction for bipolar readout
if isEddyCorrect
    disp('Step 0: Correcting eddy current effect on bipolar readout data...');
    
    % BipolarEddyCorrect requries complex-valued input
    imgCplx = BipolarEddyCorrect(magn.*exp(1i*fieldMap),mask,unwrap);
    fieldMap = angle(imgCplx);
    
    % save the eddy current corrected output
    fprintf('Saving eddy current corrected phase data...');
    save_nii_quick(outputNiftiTemplate,fieldMap,    [outputDir filesep prefix 'phase_eddy-correct.nii.gz']);
    fprintf('Done!\n');
    clear imgCplx
    
    % convert data to single type to reduce memory usage
    fieldMap = double(fieldMap);
end

% Step 1: Phase unwrapping and echo phase combination
disp('Step 1: Calculating field map...');

% fix the output field map unit in Hz
unit = 'Hz';

% core of phase unwrapping
try 
    [totalField,fieldmapSD,fieldmapUnwrapAllEchoes] = estimateTotalField(fieldMap,magn,...
                        matrixSize,voxelSize,...
                        'method',phaseCombMethod,'Unwrap',unwrap,...
                        'TE',TE,'B0',B0,'unit',unit,...
                        'Subsampling',subsampling,'mask',mask);
                    
    if ~isempty(fieldmapUnwrapAllEchoes)
        fprintf('Saving unwrapped echo phase...');
        save_nii_quick(outputNiftiTemplate,fieldmapUnwrapAllEchoes,[outputDir filesep prefix 'unwrapped-phase.nii.gz']);
        clear fieldmapUnwrapAllEchoes
        fprintf('Done!\n');
    end
    
catch
    % if the above selected method is not working then do Laplacian with
    % optimum weights
    warning('The selected method is not supported in this system. Using Laplacian algorithm for phase unwrapping.')
    if ~isinf(exclude_threshold)
        warning('Unreliable voxels will not be excluded due to switching of the algorithm.')
        exclude_threshold = Inf;
    end
    
    unwrap = 'laplacian'; 
    sepia_addpath(unwrap);
    
    [totalField,fieldmapSD] = estimateTotalField(fieldMap,magn,matrixSize,voxelSize,...
                            'Unwrap',unwrap,'TE',TE,'B0',B0,'unit',unit,...
                            'Subsampling',subsampling,'mask',mask);
end

% Step 2: exclude unreliable voxel, based on monoexponential decay model with
% single freuqnecy shift
fprintf('Computing weighting map...');
if length(TE) > 1 && ~isinf(exclude_threshold)
    % multi-echo data
    r2s = R2star_trapezoidal(magn,TE);
    relativeResidual = ComputeResidualGivenR2sFieldmap(TE,r2s,totalField,magn.*exp(1i*fieldMap));
    maskReliable = relativeResidual < exclude_threshold;
    clear r2s
else
    % single-echo & no threshold
    maskReliable = ones(size(totalField),'like',totalField);
end

% threshold fieldmapSD with the reliable voxel mask
fieldmapSD = fieldmapSD .* maskReliable;
fprintf('Done!\n');

% 20180815: test with creating weights using relativeResidual
% weightResidual = 1-(relativeResidual./exclude_threshold);
% weightResidual(weightResidual>1) = 1;
% weightResidual(weightResidual<0) = 0;

% deprecated
% mask = and(mask,maskReliable);
             
% save the output                           
fprintf('Saving unwrapped field map...');

save_nii_quick(outputNiftiTemplate,totalField,  [outputDir filesep prefix 'total-field.nii.gz']);
save_nii_quick(outputNiftiTemplate,fieldmapSD,  [outputDir filesep prefix 'noise-sd.nii.gz']);

if ~isinf(exclude_threshold)
    save_nii_quick(outputNiftiTemplate,maskReliable,   	[outputDir filesep prefix 'mask-reliable.nii.gz']);
    save_nii_quick(outputNiftiTemplate,relativeResidual,[outputDir filesep prefix 'relative-residual.nii.gz']);
    clear relativeResidual
end
fprintf('done!\n');

% clear variable that no longer be needed
clear fieldMap

%% Background field removal
disp('Step 2: Recovering local field...');

% core of background field removal
if isGPU
    localField = cuBackgroundRemovalMacro(...
                        totalField,mask,matrixSize,voxelSize,...
                        'method',BFR,'refine',refine,'erode',BFR_erode_radius,'tol',BFR_tol,'depth',BFR_depth,...
                        'peel',BFR_peel,'b0dir',B0_dir,'iteration',BFR_iteration,'padsize',...
                        BFR_padSize,'noisestd',fieldmapSD,'radius',BFR_radius,'alpha',BFR_alpha,...
                        'threshold',BFR_threshold); 
else
    localField = BackgroundRemovalMacro(...
                        totalField,mask,matrixSize,voxelSize,...
                        'method',BFR,'refine',refine,'erode',BFR_erode_radius,'tol',BFR_tol,'depth',BFR_depth,...
                        'peel',BFR_peel,'b0dir',B0_dir,'iteration',BFR_iteration,'padsize',...
                        BFR_padSize,'noisestd',fieldmapSD,'radius',BFR_radius,'alpha',BFR_alpha,...
                        'threshold',BFR_threshold); 
end
  
% generate new mask based on backgroudn field removal result
maskFinal = double(localField ~=0);
  
% save results
fprintf('Saving local field map...');

save_nii_quick(outputNiftiTemplate,localField,  [outputDir filesep prefix 'local-field.nii.gz']);
save_nii_quick(outputNiftiTemplate,maskFinal,  [outputDir filesep prefix 'mask-qsm.nii.gz']);
fprintf('done!\n');

% create weighting map based on final mask
% for weighting map: higher SNR -> higher weighting
% wmap = fieldmapSD./norm(fieldmapSD(maskFinal==1));   
if isWeightLoad
    wmap = weights;
else
    wmap = 1./fieldmapSD;
    wmap(isinf(wmap)) = 0;
    wmap(isnan(wmap)) = 0;
    wmap = wmap./max(wmap(and(maskFinal>0,maskReliable>0)));
    wmap = wmap .* maskFinal;
end
% wmap = wmap .* weightResidual;
wmap = wmap .* double(maskReliable);
if ~isWeightLoad || ~isinf(exclude_threshold)
    save_nii_quick(outputNiftiTemplate,wmap,  [outputDir filesep prefix 'weights.nii.gz']);
end

% clear variables that no longer be needed
clear maskReliable totalField mask

%% qsm
disp('Step 3: Computing QSM...');

% some QSM algorithms work better with certain unit of the local field map
switch lower(QSM_method)
    case 'closedforml2'
        
    case 'ilsqr'
        
    case 'stisuiteilsqr'
        % The order of local field values doesn't affect the result of chi  
        % in STI suite v3 implementation, i.e. 
        % chi = method(localField,...) = method(localField*C,...)/C, where
        % C is a constant.
        % Therefore, because of the scaling factor in their implementation,
        % the local field map is converted to rad
        localField = localField * 2*pi * delta_TE; 
        
    case 'fansi'
        % if both data are loaded
        if isWeightLoad && isMagnLoad
            disp('Both weighting map and magnitude images are loaded.');
            disp('Only the weighing map will be used.');
        end
        % if only magnitude images are loaded
        if ~isWeightLoad && isMagnLoad
            disp('The normalised RMS magnitude image will be used as the weighting map.');
            magn = sqrt(mean(magn.^2,4));
            wmap = magn/max(magn(:)) .* (maskFinal); 
        end
        % if nothing is loaded
        if ~isWeightLoad && ~isMagnLoad
            warning('Providing a weighing map or magnitude images can potentially improve the QSM map quality.');
        end
        
        % FANSI default parameters are optimised for ppm
        localField = localField/(B0*gyro);
        
    case 'ssvsharp'
        % not support yet
        
    case 'star'
        % Unlike the iLSQR implementation, the order of local field map
        % values will affect the Star-QSM result, i.e. 
        % chi = method(localField,...) ~= method(localField*C,...)/C, where
        % C is a constant. Lower order of local field magnitude will more 
        % likely produce chi map with streaking artefact. 
        % In the STI_Templates.m example, Star-QSM expecting local field in 
        % the unit of rad. However, value of the field map in rad will 
        % vary with echo time. Therefore, data acquired with short
        % delta_TE will be prone to streaking artefact. To mitigate this
        % potential problem, local field map is converted from Hz to radHz
        % here and the resulting chi will be normalised by the same factor 
        % 
        localField = localField * 2*pi;
        
    case 'medi_l1'
        
        if length(TE) == 1 && QSM_isLambdaCSF
            QSM_isLambdaCSF = false;
            warning('Single echo data cannot use CSF regularisation.');
            warning('No CSF regularisation will be used');
        end
            
        % zero reference using CSF requires CSF mask
        if QSM_isLambdaCSF && isMagnLoad
            sepia_addpath('medi_l1');
            
            fprintf('Extracting CSF mask....');
            % R2* mapping & computing CSF mask
            r2s = arlo(TE,magn);
            maskCSF = extract_CSF(r2s,maskFinal,voxelSize)>0;
            
            fprintf('done!');
            magn = sqrt(sum(magn.^2,4));
        end
        
        % MEDI input expects local field in rad
        localField = localField*2*pi*delta_TE;
end

% core of QSM
if isGPU
    chi = cuQSMMacro(localField,maskFinal,matrixSize,voxelSize,...
                     'method',QSM_method,'threshold',QSM_threshold,'lambda',QSM_lambda,...
                     'optimise',QSM_optimise,'tol',QSM_tol,'iteration',QSM_maxiter,'weight',wmap,...
                     'b0dir',B0_dir,'tol_step1',QSM_tol1,'tol_step2',QSM_tol2,'TE',delta_TE,'B0',B0,...
                     'padsize',QSM_padsize,'mu',QSM_mu1,'mu2',QSM_mu2,QSM_solver,QSM_constraint,...
                     'noisestd',fieldmapSD,'magnitude',magn,'data_weighting',QSM_wData,...
                     'gradient_weighting',QSM_wGradient,'merit',QSM_merit,'smv',QSM_isSMV,'zeropad',QSM_zeropad,...
                     'lambda_CSF',QSM_lambdaCSF,'CF',CF,'radius',QSM_radius,'Mask_CSF',maskCSF);
else
    chi = QSMMacro(localField,maskFinal,matrixSize,voxelSize,...
                   'method',QSM_method,'threshold',QSM_threshold,'lambda',QSM_lambda,...
                   'optimise',QSM_optimise,'tol',QSM_tol,'iteration',QSM_maxiter,'weight',wmap,...
                   'b0dir',B0_dir,'tol_step1',QSM_tol1,'tol_step2',QSM_tol2,'TE',delta_TE,'B0',B0,...
                   'padsize',QSM_padsize,'mu',QSM_mu1,'mu2',QSM_mu2,QSM_solver,QSM_constraint,...
                   'noisestd',fieldmapSD,'magnitude',magn,'data_weighting',QSM_wData,...
                   'gradient_weighting',QSM_wGradient,'merit',QSM_merit,'smv',QSM_isSMV,'zeropad',QSM_zeropad,...
                   'lambda_CSF',QSM_lambdaCSF,'CF',CF,'radius',QSM_radius,'Mask_CSF',maskCSF);
end

% convert the susceptibility map into ppm
switch lower(QSM_method)
    case 'tkd'
        chi = chi/(B0*gyro);
    case 'closedforml2'
        chi = chi/(B0*gyro);
    case 'ilsqr'
        chi = chi/(B0*gyro);
    case 'stisuiteilsqr'
        % STI suite v3 implementation already converted the chi map to ppm
    case 'fansi'
        % FANSI default parameters are optimised for ppm
    case 'ssvsharp'
        chi = chi/(B0*gyro);
    case 'star'
        % STI suite v3 implementation already nomalised the output by B0
        % and delta_TE, since the input is radHz here, we have to
        % multiply the reuslt by delta_TE here
        chi = chi * delta_TE;
    case 'medi_l1'
        % MEDI implementation already normalised the output to ppm
end
  
% save results
fprintf('Saving susceptibility map...');
save_nii_quick(outputNiftiTemplate, chi, [outputDir filesep prefix 'QSM.nii.gz']);
fprintf('done!\n');

disp('Processing pipeline is completed!');
          
end

%% check and set all algorithm parameters
function algorParam2 = CheckAndSetDefault(algorParam)
algorParam2 = algorParam;

try algorParam2.general.isInvert   	= algorParam.general.isInvert; 	catch; algorParam2.general.isInvert	= false;	end
try algorParam2.general.isBET   	= algorParam.general.isBET; 	catch; algorParam2.general.isBET   	= false;	end
try algorParam2.general.isGPU       = algorParam.general.isGPU;     catch; algorParam2.general.isGPU    = false;    end

% default method is MEDI nonlinear fitting + Laplacian + no eddy correct + no voxel exclusion
try algorParam2.unwrap.echoCombMethod       = algorParam.unwrap.echoCombMethod;         catch; algorParam2.unwrap.echoCombMethod        = 'MEDI nonlinear fit';	end
% default phase unwrapping method is Laplacian
try algorParam2.unwrap.unwrapMethod         = algorParam.unwrap.unwrapMethod;           catch; algorParam2.unwrap.unwrapMethod          = 'Laplacian';          end
try algorParam2.unwrap.isEddyCorrect        = algorParam.unwrap.isEddyCorrect;          catch; algorParam2.unwrap.isEddyCorrect         = 0;                    end
try algorParam2.unwrap.excludeMaskThreshold	= algorParam.unwrap.excludeMaskThreshold;	catch; algorParam2.unwrap.excludeMaskThreshold	= Inf;                  end
% for the rest, if the parameter does not exist then initiates it with an empty array
try algorParam2.unwrap.subsampling          = algorParam.unwrap.subsampling;            catch; algorParam2.unwrap.subsampling           = [];                   end

% default background field removal method is VSHARP
try algorParam2.bfr.method      = algorParam.bfr.method;        catch; algorParam2.bfr.method = 'vsharpsti';end
try algorParam2.bfr.radius      = algorParam.bfr.radius;        catch; algorParam2.bfr.radius = 10;         end
try algorParam2.bfr.refine      = algorParam.bfr.refine;        catch; algorParam2.bfr.refine = false;      end
try algorParam2.bfr.erode_radius= algorParam.bfr.erode_radius;	catch; algorParam2.bfr.erode_radius = 1;    end
% for the rest, if the parameter does not exist then initiates it with an empty array
try algorParam2.bfr.tol         = algorParam.bfr.tol;     	catch; algorParam2.bfr.tol = [];            end
try algorParam2.bfr.depth       = algorParam.bfr.depth;  	catch; algorParam2.bfr.depth = [];          end
try algorParam2.bfr.peel        = algorParam.bfr.peel;   	catch; algorParam2.bfr.peel = [];           end
try algorParam2.bfr.iteration   = algorParam.bfr.iteration;	catch; algorParam2.bfr.iteration = [];      end
try algorParam2.bfr.padSize     = algorParam.bfr.padSize; 	catch; algorParam2.bfr.padSize = [];        end
try algorParam2.bfr.alpha       = algorParam.bfr.alpha;    	catch; algorParam2.bfr.alpha = [];          end
try algorParam2.bfr.threshold   = algorParam.bfr.threshold;	catch; algorParam2.bfr.threshold = [];      end

% default background field removal method is TKD
try algorParam2.qsm.method      = algorParam.qsm.method;        catch; algorParam2.qsm.method       = 'tkd';	end
try algorParam2.qsm.threshold   = algorParam.qsm.threshold;     catch; algorParam2.qsm.threshold    = 0.15;     end
% for the rest, if the parameter does not exist then initiates it with an empty array
try algorParam2.qsm.radius      = algorParam.qsm.radius;        catch; algorParam2.qsm.radius       = [];       end
try algorParam2.qsm.lambda      = algorParam.qsm.lambda;        catch; algorParam2.qsm.lambda       = [];       end
try algorParam2.qsm.optimise   	= algorParam.qsm.optimise;     	catch; algorParam2.qsm.optimise     = [];       end
try algorParam2.qsm.maxiter   	= algorParam.qsm.maxiter;       catch; algorParam2.qsm.maxiter      = [];     	end
try algorParam2.qsm.tol1        = algorParam.qsm.tol1;          catch; algorParam2.qsm.tol1         = [];   	end
try algorParam2.qsm.tol2        = algorParam.qsm.tol2;          catch; algorParam2.qsm.tol2         = [];      	end
try algorParam2.qsm.padsize     = algorParam.qsm.padsize;       catch; algorParam2.qsm.padsize      = [];   	end
try algorParam2.qsm.tol         = algorParam.qsm.tol;           catch; algorParam2.qsm.tol          = [];      	end
try algorParam2.qsm.mu1         = algorParam.qsm.mu1;           catch; algorParam2.qsm.mu1          = [];      	end
try algorParam2.qsm.mu2         = algorParam.qsm.mu2;           catch; algorParam2.qsm.mu2          = [];    	end
try algorParam2.qsm.solver    	= algorParam.qsm.solver;    	catch; algorParam2.qsm.solver       = [];    	end
try algorParam2.qsm.constraint 	= algorParam.qsm.constraint;   	catch; algorParam2.qsm.constraint   = [];       end
try algorParam2.qsm.wData       = algorParam.qsm.wData;         catch; algorParam2.qsm.wData        = [];       end
try algorParam2.qsm.wGradient  	= algorParam.qsm.wGradient;    	catch; algorParam2.qsm.wGradient    = [];       end
try algorParam2.qsm.zeropad  	= algorParam.qsm.zeropad;    	catch; algorParam2.qsm.zeropad      = [];      	end
try algorParam2.qsm.isSMV    	= algorParam.qsm.isSMV;         catch; algorParam2.qsm.isSMV        = [];    	end
try algorParam2.qsm.isLambdaCSF	= algorParam.qsm.isLambdaCSF;	catch; algorParam2.qsm.isLambdaCSF  = [];   	end
try algorParam2.qsm.lambdaCSF 	= algorParam.qsm.lambdaCSF;    	catch; algorParam2.qsm.lambdaCSF    = [];     	end
try algorParam2.qsm.merit       = algorParam.qsm.merit;         catch; algorParam2.qsm.merit        = [];      	end

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