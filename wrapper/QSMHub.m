%% [chi,localField,totalField,fieldmapSD]=QSMHub(inputDir,outputDir,varargin)
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
%               (1) matches the input format of qsm_hub.m
%               (2) save the results in NIfTI format
%
% Kwok-shing Chan @ DCCN
% k.chan@donders.ru.nl
% Date created: 14 September 2017
% Date last modified: 26 August 2018
%
%
function [chi,localField,totalField,fieldmapSD]=QSMHub(input,output,maskFullName,varargin)
%% add general Path
qsm_hub_AddMethodPath

%% define variables
prefix = 'squirrel_';   
gyro = 42.57747892;
isInputDir = true;
% make sure the input only load once (first one)
isMagnLoad = false;
isPhaseLoad = false;
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

%% Parse input argument
[isInvert,isGPU,isBET,isEddyCorrect,...
    phaseCombMethod,unwrap,subsampling,exclude_threshold,...
    BFR,refine,BFR_tol,BFR_depth,BFR_peel,BFR_iteration,BFR_padSize,BFR_radius,BFR_alpha,BFR_threshold,...
    QSM_method,QSM_threshold,QSM_lambda,QSM_optimise,QSM_tol,QSM_maxiter,...
    QSM_tol1,QSM_tol2,QSM_padsize,QSM_mu1,QSM_mu2,QSM_solver,QSM_constraint,...
    QSM_radius,QSM_zeropad,QSM_wData,QSM_wGradient,QSM_isLambdaCSF,QSM_lambdaCSF,...
    QSM_isSMV,QSM_merit] = parse_varargin_QSMHub(varargin);

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
        
                        %%%%%%%%%% magnitude data %%%%%%%%%%
        if ~isempty(inputNiftiList(2).name)
            inputMagnNifti = load_untouch_nii([inputNiftiList(2).name]);
            % make sure the data is multi-echo magnitude data
            if size(inputMagnNifti.img,4) > 1
                magn = double(inputMagnNifti.img);
                isMagnLoad = true;
                disp('Magnitude data is loaded.');
            else
                error('QSM Hub only works with 4D data.');
            end
        else
            error('Please specify a 4D magnitude data.');
        end
        
                         %%%%%%%%%% phase data %%%%%%%%%%
        if ~isempty(inputNiftiList(1).name)
            inputPhaseNifti = load_untouch_nii([inputNiftiList(1).name]);
            fieldMap = double(inputPhaseNifti.img);
            % check whether phase data contains DICOM values or wrapped
            % phase value
            if max(fieldMap(:))>1000
                disp('Converting phase data from DICOM image value to radian unit...')
                fieldMap = DICOM2Phase(inputPhaseNifti);

                disp('Saving phase images in unit of radian...');
                save_nii_quick(inputPhaseNifti,fieldMap, [outputDir filesep prefix 'phase.nii.gz']);

            end
            isPhaseLoad = true;
            disp('Phase data is loaded.')
        else
            error('Please specify a 4D Phase data.');
        end
        
                        %%%%%%%%%% qsm hub header %%%%%%%%%%
        if ~isempty(inputNiftiList(4).name)
            load([inputNiftiList(4).name]);
            disp('Header data is loaded.');
        else
            error('Please specify a qsm_hub header.');
        end
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    else
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %%%%%%%%%%% Pathway 2: Input is a directory with NIfTI %%%%%%%%%%%%
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        
        % loop all NIfTI files in the directory for magnitude and phase files
        for klist = 1:length(inputNiftiList)
                            %%%%%%%%%% magnitude data %%%%%%%%%%
            if ContainName(lower(inputNiftiList(klist).name),'magn') && ~isMagnLoad
                inputMagnNifti = load_untouch_nii([inputDir filesep inputNiftiList(klist).name]);
                % only load multi-echo magnitude data
                if size(inputMagnNifti.img,4) > 1
                    magn = double(inputMagnNifti.img);
                    isMagnLoad = true;
                    disp('Magnitude data is loaded.')
                end
            end

                            %%%%%%%%%% phase data %%%%%%%%%%
            if ContainName(lower(inputNiftiList(klist).name),'phase') && ~isPhaseLoad
                inputPhaseNifti = load_untouch_nii([inputDir filesep inputNiftiList(klist).name]);
                fieldMap = double(inputPhaseNifti.img);

                % if input fieldmap is directly converted from nifti converter
                % then converts the fieldmap to unit of radian and save to output dir
                if max(fieldMap(:))>1000
                    disp('Converting phase data from DICOM image value to radian unit ...')
                    fieldMap = DICOM2Phase(inputPhaseNifti);

                    disp('Saving phase images in unit of radian...');
                    save_nii_quick(inputPhaseNifti,fieldMap, [outputDir filesep prefix 'phase.nii.gz']);

                end
                isPhaseLoad = true;
                disp('Phase data is loaded.')
            end
        end

        % if no files matched the name format then displays error message
        if ~isMagnLoad
            error('No magnitude data is loaded. Please make sure the input directory contains NIfTI files with name *magn*');
        end
        if ~isPhaseLoad
            error('No phase data is loaded. Please make sure the input directory contains files with name *phase*');
        end

                        %%%%%%%%%% qsm hub header file %%%%%%%%%%
        if ~isempty(dir([inputDir '/*header*']))
            % load header
            headerList = dir([inputDir '/*header*']);
            load([inputDir filesep headerList(1).name]);

            disp('Header data is loaded.');

        else
            disp('No header for qsm_hub is found. Creating synthetic header based on NIfTI header...');

            % create synthetic header in case no qsm_hub's header is found
            [B0,B0_dir,voxelSize,matrixSize,TE,delta_TE,CF]=SyntheticQSMHubHeader(inputMagnNifti);

            % look for text file for TEs information
            teTextFullName = dir([inputDir filesep '*txt']);
            if ~isempty(teTextFullName)
                te_ = readTEfromText([inputDir filesep teTextFullName(1).name]);
                te_ = te_(:);
                if ~isempty(te_)
                    TE = te_;
                    if length(TE) > 1
                        delta_TE = TE(2)-TE(1);
                    end
                end
            end

            % if no header file then save the synthetic header in output dir
            save([outputDir filesep 'SyntheticQSMhub_header'],'voxelSize','matrixSize','CF','delta_TE',...
            'TE','B0_dir','B0');

            disp('The synthetic header is saved in output directory.');

        end
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    end
    
    % store the header of the NIfTI files, all following results will have
    % the same header
    outputNiftiTemplate = inputMagnNifti;
    
else
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%%%%%%%%%% Pathway 3: Input is a directory with DICOM %%%%%%%%%%%%%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    qsm_hub_AddMethodPath('dicom');
    
    % if no nifti file then check for DICOM files
    [iField,voxelSize,matrixSize,CF,delta_TE,TE,B0_dir]=Read_DICOM(inputDir);
    
    % deprecated
    %     [iField,voxelSize,matrixSize,CF,delta_TE,TE,B0_dir]=Read_Siemens_DICOM_old(inputDir);
    
    B0 = CF/(gyro*1e6);
    
    % after testing with a couple of dataset it seems to me that the field
    % is inverted with DICOM input, so apply conjugate here
    fieldMap = angle(conj(iField));
    magn = abs(iField);
    
    isMagnLoad = true;
    isPhaseLoad = true;

    % save magnitude and phase images as nifti files
    disp('Saving DICOM data into NIfTI format...');
    
    % save magnitude and phase data as NIfTI_GZ format
    nii_fieldMap    = make_nii(single(fieldMap),    voxelSize,[],16);
    nii_magn        = make_nii(single(magn),        voxelSize,[],16);
    save_nii(nii_fieldMap,  [outputDir filesep prefix 'phase.nii.gz']);
    save_nii(nii_magn,      [outputDir filesep prefix 'magn.nii.gz']);
    % save important header in .mat format
    save([outputDir filesep prefix 'header.mat'],'voxelSize','matrixSize','CF','delta_TE',...
        'TE','B0_dir','B0');
    
    % reload the NIfTI template so that later can use save_untouch_nii for
    % all results
    outputNiftiTemplate = load_untouch_nii([outputDir filesep prefix 'magn.nii.gz']);
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
end

% in case user want to reverse the frequency shift direction
if isInvert
    fieldMap = -fieldMap;
end

% display some header info
disp('Basic DICOM information');
disp(['Voxel size(x,y,z mm^3) =  ' num2str(voxelSize(1)) 'x' num2str(voxelSize(2)) 'x' num2str(voxelSize(3))]);
disp(['matrix size(x,y,z) =  ' num2str(matrixSize(1)) 'x' num2str(matrixSize(2)) 'x' num2str(matrixSize(3))]);
disp(['B0 direction(x,y,z) =  ' num2str(B0_dir(:)')]);
disp(['Field strength(T) =  ' num2str(B0)]);

% make sure the following variables are row vectors
matrixSize = matrixSize(:).';
voxelSize = voxelSize(:).';

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
    qsm_hub_AddMethodPath('bet');
    disp('Performing FSL BET...');
    
    % Here uses MEDI toolboxes MEX implementation
    mask = BET(magn(:,:,:,1),matrixSize,voxelSize);
    
    disp('Saving brain mask...')
    save_nii_quick(outputNiftiTemplate,mask, [outputDir filesep prefix 'mask.nii.gz']);

end

%% total field and phase unwrap

% Step 0: Eddy current correction for bipolar readout
if isEddyCorrect
    disp('Correcting eddy current effect on bipolar readout data...');
    
    % BipolarEddyCorrect requries complex-valued input
    imgCplx = BipolarEddyCorrect(magn.*exp(1i*fieldMap),mask,unwrap);
    fieldMap = angle(imgCplx);
    magn = abs(imgCplx);
    
    % save the eddy current corrected output
    save_nii_quick(outputNiftiTemplate,fieldMap,    [outputDir filesep prefix 'phase_eddy-correct.nii.gz']);
    save_nii_quick(outputNiftiTemplate,magn,        [outputDir filesep prefix 'magn_eddy-correct.nii.gz']);
    
end

% Step 1: Phase unwrapping and echo phase combination
disp('Calculating field map...');

% fix the output field map unit in Hz
unit = 'Hz';

% core of phase unwrapping
try 
    switch phaseCombMethod
        % optimum weight method from Robinson et al. NMR Biomed 2017
        case 'Optimum weights'
            [totalField,fieldmapSD] = estimateTotalField(fieldMap,magn,...
                                    matrixSize,voxelSize,...
                                    'Unwrap',unwrap,'TE',TE,'B0',B0,'unit',unit,...
                                    'Subsampling',subsampling,'mask',mask);
                                
        % nonlinear ffitting method from MEDI toolbox
        case 'MEDI nonlinear fit'
            [totalField,fieldmapSD] = estimateTotalField_MEDI(fieldMap,magn,...
                                    matrixSize,voxelSize,...
                                    'Unwrap',unwrap,'TE',TE,'B0',B0,'unit',unit,...
                                    'Subsampling',subsampling,'mask',mask);
    end
    
catch
    % if the above selected method is not working then do Laplacian with
    % optimum weights
    disp('The selected method is not supported in this system. Using Laplacian algorithm for phase unwrapping...')
    
    unwrap = 'laplacian'; exclude_threshold = Inf;
    qsm_hub_AddMethodPath(unwrap);
    
    [totalField,fieldmapSD] = estimateTotalField(fieldMap,magn,matrixSize,voxelSize,...
                            'Unwrap',unwrap,'TE',TE,'B0',B0,'unit',unit,...
                            'Subsampling',subsampling,'mask',mask);
end

% Step 2: exclude unreliable voxel, based on monoexponential decay model with
% single freuqnecy shift
r2s = R2star_trapezoidal(magn,TE);
relativeResidual = ComputeResidualGivenR2sFieldmap(TE,r2s,totalField,magn.*exp(1i*fieldMap));

maskReliable = relativeResidual < exclude_threshold;

mask = and(mask,maskReliable);
             
% save the output                           
disp('Saving unwrapped field map...');

save_nii_quick(outputNiftiTemplate,totalField,  [outputDir filesep prefix 'total-field.nii.gz']);
save_nii_quick(outputNiftiTemplate,fieldmapSD,  [outputDir filesep prefix 'fieldmap-sd.nii.gz']);

if relativeResidual ~= Inf
    save_nii_quick(outputNiftiTemplate,mask,        [outputDir filesep prefix 'mask-reliable.nii.gz']);
    save_nii_quick(outputNiftiTemplate,relativeResidual,[outputDir filesep prefix 'relative-residual.nii.gz']);
end

%% Background field removal
disp('Recovering local field...');

% core of background field removal
if isGPU
    localField = cuBackgroundRemovalMacro(...
                        totalField,mask,matrixSize,voxelSize,...
                        'method',BFR,'refine',refine,'tol',BFR_tol,'depth',BFR_depth,...
                        'peel',BFR_peel,'b0dir',B0_dir,'iteration',BFR_iteration,'padsize',...
                        BFR_padSize,'noisestd',fieldmapSD,'radius',BFR_radius,'alpha',BFR_alpha,...
                        'threshold',BFR_threshold); 
else
    localField = BackgroundRemovalMacro(...
                        totalField,mask,matrixSize,voxelSize,...
                        'method',BFR,'refine',refine,'tol',BFR_tol,'depth',BFR_depth,...
                        'peel',BFR_peel,'b0dir',B0_dir,'iteration',BFR_iteration,'padsize',...
                        BFR_padSize,'noisestd',fieldmapSD,'radius',BFR_radius,'alpha',BFR_alpha,...
                        'threshold',BFR_threshold); 
end
  
% generate new mask based on backgroudn field removal result
maskFinal = localField ~=0;
  
% save results
disp('Saving local field map...');

save_nii_quick(outputNiftiTemplate,localField,  [outputDir filesep prefix 'local-field.nii.gz']);
save_nii_quick(outputNiftiTemplate,maskFinal,  [outputDir filesep prefix 'mask-qsm.nii.gz']);
            
% create weighting map based on final mask
% for weighting map: higher SNR -> higher weighting
% wmap = fieldmapSD./norm(fieldmapSD(maskFinal==1));    
wmap = 1./fieldmapSD;
wmap(isinf(wmap)) = 0;
wmap(isnan(wmap)) = 0;
wmap = wmap./max(wmap(maskFinal>0));


%% qsm
disp('Computing QSM...');

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
        % zero reference using CSF requires CSF mask
        if QSM_isLambdaCSF && isMagnLoad
            disp('Extracting CSF mask....');
            
            % R2* mapping
            r2s = arlo(TE,magn);
            maskCSF = extract_CSF(r2s,maskFinal,voxelSize)>0;
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
                     'noisestd',fieldmapSD,'magnitude',sqrt(sum(magn.^2,4)),'data_weighting',QSM_wData,...
                     'gradient_weighting',QSM_wGradient,'merit',QSM_merit,'smv',QSM_isSMV,'zeropad',QSM_zeropad,...
                     'lambda_CSF',QSM_lambdaCSF,'CF',CF,'radius',QSM_radius,'Mask_CSF',maskCSF);
else
    chi = QSMMacro(localField,maskFinal,matrixSize,voxelSize,...
                   'method',QSM_method,'threshold',QSM_threshold,'lambda',QSM_lambda,...
                   'optimise',QSM_optimise,'tol',QSM_tol,'iteration',QSM_maxiter,'weight',wmap,...
                   'b0dir',B0_dir,'tol_step1',QSM_tol1,'tol_step2',QSM_tol2,'TE',delta_TE,'B0',B0,...
                   'padsize',QSM_padsize,'mu',QSM_mu1,'mu2',QSM_mu2,QSM_solver,QSM_constraint,...
                   'noisestd',fieldmapSD,'magnitude',sqrt(sum(magn.^2,4)),'data_weighting',QSM_wData,...
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
disp('Saving susceptibility map...');

save_nii_quick(outputNiftiTemplate, chi, [outputDir filesep prefix 'QSM.nii.gz']);

disp('Done!');
          
end
