%% chi = qsmMacroIOWrapper(inputDir,outputDir,varargin)
%
% Input
% --------------
% inputDir              : input directory contains NIfTI (*localfield*, *magn* and *fieldmapsd*) files 
% outputDir             : output directory that stores the output (susceptibility map)
% varargin ('Name','Value' pair)
% ---------
% 'mask'                : mask file (in NIfTI format) full name 
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
% chi                   : quantitative susceptibility map (in ppm)
%
% Description: This is a wrapper of BackgroundRemovalMac1ro.m which has the following objeectives:
%               (1) matches the input format of qsm_hub.m
%               (2) save the results in NIfTI format
%
% Kwok-shing Chan @ DCCN
% k.chan@donders.ru.nl
% Date created: 17 April 2018
% Date modified: 26 August 2018
% Date modified: 29 March 2019
% Date modified: 5 June 2019
% Date modified: 27 Feb 2020 (v0.8.0)
%
%
function chi = QSMMacroIOWrapper(input,output,maskFullName,algorParam)
%% add general Path
sepia_addpath

%% define variables
prefix = 'sepia_';
gyro = 42.57747892;
isInputDir = true;
% make sure the input only load once (first one)
isLotalFieldLoad    = false;
isWeightLoad        = false;
isMagnLoad          = false; 
maskCSF             = [];

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
algorParam      = check_and_set_SEPIA_algorithm_default(algorParam);
isGPU           = algorParam.general.isGPU;
QSM_method      = algorParam.qsm.method; 
QSM_threshold   = algorParam.qsm.threshold; 
QSM_lambda      = algorParam.qsm.lambda; 
QSM_optimise 	= algorParam.qsm.optimise; 
QSM_tol         = algorParam.qsm.tol;   
QSM_maxiter     = algorParam.qsm.maxiter;
QSM_tol1        = algorParam.qsm.tol1;  
QSM_tol2        = algorParam.qsm.tol2; 
QSM_padsize     = algorParam.qsm.padsize;
QSM_mu1         = algorParam.qsm.mu1; 
QSM_mu2         = algorParam.qsm.mu2;  
QSM_solver      = algorParam.qsm.solver;  
QSM_constraint	= algorParam.qsm.constraint; 
QSM_radius      = algorParam.qsm.radius;
QSM_zeropad     = algorParam.qsm.zeropad;   
QSM_wData       = algorParam.qsm.wData; 
QSM_wGradient 	= algorParam.qsm.wGradient;
QSM_isLambdaCSF	= algorParam.qsm.isLambdaCSF;
QSM_lambdaCSF	= algorParam.qsm.lambdaCSF; 
QSM_isSMV       = algorParam.qsm.isSMV;
QSM_merit       = algorParam.qsm.merit;  
QSM_stepSize   	= algorParam.qsm.stepSize;  
QSM_percentage  = algorParam.qsm.percentage;  

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
        
                        %%%%%%%%%% Local field map %%%%%%%%%% 
        if ~isempty(inputNiftiList(1).name)
            inputLocalFieldNifti = load_untouch_nii([inputNiftiList(1).name]);
            % load true value from NIfTI
            localField = load_nii_img_only([inputNiftiList(1).name]);
%             localField = double(inputLocalFieldNifti.img);
            isLotalFieldLoad = true;
            
            disp('Local field map is loaded.')
        else
            error('Please specify a 3D local field map.');
        end
        
                        %%%%%%%%%% magnitude data %%%%%%%%%%
        if ~isempty(inputNiftiList(2).name)
%             inputMagnNifti = load_untouch_nii([inputNiftiList(2).name]);
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
%             inputWeightNifti = load_untouch_nii([inputNiftiList(3).name]);
            % load true value from NIfTI
            weights = load_nii_img_only([inputNiftiList(3).name]);
%             weights = double(inputWeightNifti.img);
            % check whether phase data contains DICOM values or wrapped
            if size(weights,4) > 1
                error('Please specify a 3D weight data.');
            end
            isWeightLoad = true;
            disp('Weights data is loaded');
        else
            disp('Default weighting method will be used for QSM.');
        end
        
                        %%%%%%%%%% qsm hub header %%%%%%%%%%
        if ~isempty(inputNiftiList(4).name)
            load([inputNiftiList(4).name]);
            disp('Header data is loaded.');
        else
            error('Please specify a header required by SEPIA.');
        end
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    else
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %%%%%%%%%%% Pathway 2: Input is a directory with NIfTI %%%%%%%%%%%%
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % validate indput directory files
        disp('Directory input is being used.');
        fprintf('Validating filenames in the directory....');
        CheckFileName(inputNiftiList);
        fprintf('Filenames are valid.\n');
        
        % loop all NIfTI files in the directory for magnitude,localField and fieldmapSD
        for klist = 1:length(inputNiftiList)

                        %%%%%%%%%% Local field map %%%%%%%%%%
            if ContainName(lower(inputNiftiList(klist).name),'local-field') && ~isLotalFieldLoad
                inputLocalFieldNifti = load_untouch_nii([inputDir filesep inputNiftiList(klist).name]);
                % load true value from NIfTI
                localField = load_nii_img_only([inputDir filesep inputNiftiList(klist).name]);
%                 localField = double(inputLocalFieldNifti.img);
                isLotalFieldLoad = true;
                disp('Local field map is loaded.')
            end

                        %%%%%%%%%% magnitude data %%%%%%%%%%
            if ContainName(lower(inputNiftiList(klist).name),'mag') && ~ContainName(lower(inputNiftiList(klist).name),'brain') && ~isMagnLoad
%                 inputMagnNifti = load_untouch_nii([inputDir filesep inputNiftiList(klist).name]);
                % load true value from NIfTI
                magn = load_nii_img_only([inputDir filesep inputNiftiList(klist).name]);
                isMagnLoad = true;
%                 magn = double(inputMagnNifti.img);
                disp('Magnitude data is loaded.')
            end

                        %%%%%%%%%% weights map %%%%%%%%%%
            if ContainName(lower(inputNiftiList(klist).name),'weights') && ~isWeightLoad
%                 inputWeightNifti = load_untouch_nii([inputDir filesep inputNiftiList(klist).name]);
                % load true value from NIfTI
                weights = load_nii_img_only([inputDir filesep inputNiftiList(klist).name]);
%                 weights = double(inputWeightNifti.img);
                isWeightLoad = true;
                disp('Weights data is loaded.')
            end

        end

        % if no files matched the name format then displays error message
        if ~isLotalFieldLoad
            error('No local field map is loaded. Please make sure the input directory contains files with name *localfield*');
        end

                    %%%%%%%%%% qsm hub header file %%%%%%%%%%
        if ~isempty(dir([inputDir '/*header*']))
            % load header
            headerList = dir([inputDir '/*header*']);
            load([inputDir filesep headerList(1).name]);

            disp('Header data is loaded.');

        else
            error('Please specify a header required by Sepia.');
        end

        % if no magnitude found then creates one with all voxels have the same value
        if ~isMagnLoad && strcmpi(QSM_method,'medi_l1')
            error('MEDI requires magnitude data. Please put the magnitude multi-echo data to input directory or use other algorithm');
        elseif ~isMagnLoad
            disp('No magnitude data is loaded.');
            magn = ones(matrixSize,'like',matrixSize);
        end
    end
    
    % if no fieldmapSD found then creates one with all voxels have the same value
    if ~isWeightLoad
        disp('No weights file is loaded.');
%         fieldmapSD = ones(matrixSize) * 0.01;
        weights = ones(matrixSize);
    end
    
    % store the header the NIfTI files, all following results will have
    % the same header
    outputNiftiTemplate = inputLocalFieldNifti;
    
else
    error('This standalone only reads NIfTI format input data (*.nii or *.nii.gz).');
end

% make sure the L2 norm of B0 direction = 1
B0_dir = B0_dir ./ norm(B0_dir);

% display some header info
disp('Basic DICOM information');
disp(['Voxel size(x,y,z mm^3)   =  ' num2str(voxelSize(1)) 'x' num2str(voxelSize(2)) 'x' num2str(voxelSize(3))]);
disp(['matrix size(x,y,z)       =  ' num2str(matrixSize(1)) 'x' num2str(matrixSize(2)) 'x' num2str(matrixSize(3))]);
disp(['B0 direction(x,y,z)      =  ' num2str(B0_dir(:)')]);
disp(['Field strength(T)        =  ' num2str(B0)]);

%% get brain mask
% look for qsm mask first
maskList = dir([inputDir '/*mask-qsm*']);
% if no final mask then just look for normal mask
if isempty(maskList)
    maskList = dir([inputDir '/*mask*']);
end

if ~isempty(maskFullName)
    % Option 1: mask file is provided
    maskFinal = load_nii_img_only(maskFullName) > 0;
    
elseif ~isempty(maskList) 
    % Option 2: input directory contains NIfTI file with name '*mask*'
    inputMaskNii = load_untouch_nii([inputDir filesep maskList(1).name]);
	maskFinal = inputMaskNii.img > 0;
    
else
    % display error message if nothing is found
    error('No mask file is loaded. Pleasee specific your mask file or put it in the input directory.');
    
end

% make sure all variables are double
localField	= double(localField);
maskFinal   = double(maskFinal);
voxelSize   = double(voxelSize);
matrixSize  = double(matrixSize);
if exist('wmap','var'); weights = double(weights);  end
if exist('magn','var'); magn = double(magn);        end

% create weighting map based on final mask
% for weighting map: higher SNR -> higher weighting
weights = weights .* maskFinal;

%% qsm
disp('Computing QSM...');

% prepare all essential components for individual algorithm
switch lower(QSM_method)
    case 'closedforml2'
        
    case 'ilsqr'
        
    case 'stisuiteilsqr'
        
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
            weights = (magn./max(magn(:))) .* (maskFinal); 
        end
        % if nothing is loaded
        if ~isWeightLoad && ~isMagnLoad
            warning('Providing a weighing map or magnitude images can potentially improve the QSM map quality.');
        end
        
    case 'ssvsharp'
        % not support yet
        
    case 'star'
%         
    case 'medi_l1'
        % zero reference using CSF requires CSF mask
        if QSM_isLambdaCSF && isMagnLoad
            disp('Extracting CSF mask....');
            sepia_addpath('medi_l1');
            % R2* mapping
            r2s = arlo(TE,magn);
            maskCSF = extract_CSF(r2s,maskFinal,voxelSize)>0;
            magn = sqrt(sum(magn.^2,4));
        end
        
    case 'ndi'
        % if both data are loaded
        if isWeightLoad && isMagnLoad
            disp('Both weighting map and magnitude images are loaded.');
            disp('Only the weighing map will be used.');
        end
        % if only magnitude images are loaded
        if ~isWeightLoad && isMagnLoad
            disp('The normalised RMS magnitude image will be used as the weighting map.');
            magn = sqrt(mean(magn.^2,4));
            weights = (magn./max(magn(:))) .* (maskFinal); 
        end
        % if nothing is loaded
        if ~isWeightLoad && ~isMagnLoad
            warning('Providing a weighing map or magnitude images can potentially improve the QSM map quality.');
        end
        
end

% core of QSM
if isGPU
    chi = cuQSMMacro(localField,maskFinal,matrixSize,voxelSize,...
                     'method',QSM_method,'threshold',QSM_threshold,'lambda',QSM_lambda,...
                     'optimise',QSM_optimise,'tol',QSM_tol,'iteration',QSM_maxiter,'weight',weights,...
                     'b0dir',B0_dir,'tol_step1',QSM_tol1,'tol_step2',QSM_tol2,'TE',delta_TE,'B0',B0,...
                     'padsize',QSM_padsize,'mu',QSM_mu1,'mu2',QSM_mu2,QSM_solver,QSM_constraint,...
                     'noisestd',weights,'magnitude',magn,'data_weighting',QSM_wData,...
                     'gradient_weighting',QSM_wGradient,'merit',QSM_merit,'smv',QSM_isSMV,'zeropad',QSM_zeropad,...
                     'lambda_CSF',QSM_lambdaCSF,'CF',CF,'radius',QSM_radius,'Mask_CSF',maskCSF);
else
    chi = QSMMacro(localField,maskFinal,matrixSize,voxelSize,...
                   'method',QSM_method,'threshold',QSM_threshold,'lambda',QSM_lambda,...
                   'optimise',QSM_optimise,'tol',QSM_tol,'iteration',QSM_maxiter,'weight',weights,...
                   'b0dir',B0_dir,'tol_step1',QSM_tol1,'tol_step2',QSM_tol2,'TE',delta_TE,'B0',B0,...
                   'padsize',QSM_padsize,'mu',QSM_mu1,'mu2',QSM_mu2,QSM_solver,QSM_constraint,...
                   'noisestd',weights,'magnitude',magn,'data_weighting',QSM_wData,...
                   'gradient_weighting',QSM_wGradient,'merit',QSM_merit,'smv',QSM_isSMV,'zeropad',QSM_zeropad,...
                   'lambda_CSF',QSM_lambdaCSF,'CF',CF,'radius',QSM_radius,'Mask_CSF',maskCSF,...
                   'stepsize',QSM_stepSize,'percentage',QSM_percentage,'tmp_output_dir',outputDir);
end

% save results
fprintf('Saving susceptibility map...');

save_nii_quick(outputNiftiTemplate, chi, [outputDir filesep prefix 'QSM.nii.gz']);

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
