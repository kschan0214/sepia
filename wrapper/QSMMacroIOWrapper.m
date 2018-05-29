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
% Date last modified: 27 May 2018
%
%
function chi = QSMMacroIOWrapper(inputDir,outputDir,varargin)
%% add general Path
qsm_hub_AddMethodPath % qsm_hub_AddPath;

%% define variables
prefix = 'squirrel_';
gyro = 42.57747892;
% make sure the input only load once (first one)
isLotalFieldLoad    = false;
isFieldmapSDLoad    = false;
isMagnLoad          = false; 
maskCSF             = [];


%% Check output directory exist or not
if exist(outputDir,'dir') ~= 7
    % if not then create the directory
    mkdir(outputDir);
end

%% Parse input argument
[~,maskFullName,~,~,~,~,~,~,~,~,~,~,~,~,~,QSM_method,QSM_threshold,QSM_lambda,...
QSM_optimise,QSM_tol,QSM_maxiter,QSM_tol1,QSM_tol2,QSM_padsize,QSM_mu1,QSM_mu2,QSM_solver,QSM_constraint,...
~,QSM_radius,QSM_zeropad,QSM_wData,QSM_wGradient,QSM_isLambdaCSF,QSM_lambdaCSF,QSM_isSMV,QSM_merit,~,isGPU] = parse_varargin_QSMHub(varargin);

%% Read input
disp('Reading data...');

% look for nifti files first
inputNiftiList = dir([inputDir '/*.nii*']);
if ~isempty(inputNiftiList)
    
    % look for magnitude,localField and fieldmapSD files
    for klist = 1:length(inputNiftiList)
        
        if ContainName(lower(inputNiftiList(klist).name),'localfield') && ~isLotalFieldLoad
            inputLocalFieldNifti = load_untouch_nii([inputDir filesep inputNiftiList(klist).name]);
            localField = double(inputLocalFieldNifti.img);
            isLotalFieldLoad = true;
        end
        
        if ContainName(lower(inputNiftiList(klist).name),'magn') && ~ContainName(lower(inputNiftiList(klist).name),'brain') && ~isMagnLoad
            inputMagnNifti = load_untouch_nii([inputDir filesep inputNiftiList(klist).name]);
            isMagnLoad = true;
            magn = double(inputMagnNifti.img);
        end
        
        if ContainName(lower(inputNiftiList(klist).name),'fieldmapsd') && ~isFieldmapSDLoad
            inputFieldMapSDNifti = load_untouch_nii([inputDir filesep inputNiftiList(klist).name]);
            fieldmapSD = double(inputFieldMapSDNifti.img);
            isFieldmapSDLoad = true;
        end
        
    end
    
    % if no files matched the name format then displays error message
    if ~isLotalFieldLoad
        error('No files loaded. Please make sure the input directory contains files with name *localfield*');
    end
    
    % look for header file
    if ~isempty(dir([inputDir '/*header*']))
        headerList = dir([inputDir '/*header*']);
        
        % load header
        load([inputDir filesep headerList(1).name]);
        
    else
        disp('No header for qsm_hub is found. Creating synthetic header based on NIfTI header...');
        
        % create synthetic header in case no qsm_hub's header is found
        [B0,B0_dir,voxelSize,matrixSize,TE,delta_TE,CF]=SyntheticQSMHubHeader(inputLocalFieldNifti);
        
    end
    
    % if no magnitude found then creates one with all voxels have the same value
    if ~isMagnLoad && strcmpi(QSM_method,'medi_l1')
        error('MEDI requires magnitude data. Please put the magnitude multi-echo data to input directory or use other algorithm');
    elseif ~isMagnLoad
        disp('No magnitude data is found.');
        magn = ones(matrixSize);
    end
    
    % if no fieldmapSD found then creates one with all voxels have the same value
    if ~isFieldmapSDLoad
        disp('No fieldMapSD file is found.');
        fieldmapSD = ones(matrixSize) * 0.01;
    end
    
    % store the header the NIfTI files, all following results will have
    % the same header
    outputNiftiTemplate = inputLocalFieldNifti;
    % make sure the class of output datatype is double
    outputNiftiTemplate.hdr.dime.datatype = 64;
    % remove the time dimension info
    outputNiftiTemplate.hdr.dime.dim(5) = 1;
    
    
else
    error('This standalone only reads NIfTI format input data (*.nii or *.nii.gz).');
end

% display some header info
disp('Basic DICOM information');
disp(['Voxel size(x,y,z mm^3) =  ' num2str(voxelSize(1)) 'x' num2str(voxelSize(2)) 'x' num2str(voxelSize(3))]);
disp(['matrix size(x,y,z) =  ' num2str(matrixSize(1)) 'x' num2str(matrixSize(2)) 'x' num2str(matrixSize(3))]);
disp(['B0 direction(x,y,z) =  ' num2str(B0_dir(:)')]);
disp(['Field strength(T) =  ' num2str(B0)]);

%% get brain mask
% look for final mask first
maskList = dir([inputDir '/*mask_final*']);
% if no final mask then just look for normal mask
if isempty(maskList)
    maskList = dir([inputDir '/*mask*']);
end

if ~isempty(maskFullName)
% first read mask if file is provided
    maskFinal = load_nii_img_only(maskFullName) > 0;
    
elseif ~isempty(maskList) 
% read mask if input directory contains *mask*
    inputMaskNii = load_untouch_nii([inputDir filesep maskList(1).name]);
	maskFinal = inputMaskNii.img > 0;
    
else
    % display error message if nothing is found
    error('No mask is found. Pleasee specific your mask file or put it inside the input directory.');
    
end

% create weighting map based on final mask
% for weighting map: higher SNR -> higher weighting
% wmap = fieldmapSD./norm(fieldmapSD(maskFinal==1));    
wmap = 1./fieldmapSD;
wmap(isinf(wmap)) = 0;
wmap(isnan(wmap)) = 0;
wmap = wmap./max(wmap(maskFinal>0));

%% qsm
qsm_hub_AddMethodPath(QSM_method);
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

nii_chi = make_nii_quick(outputNiftiTemplate,chi);

save_untouch_nii(nii_chi,[outputDir filesep prefix 'QSM.nii.gz']);

disp('Done!');

end

% handy function to save result to nifti format
function nii = make_nii_quick(template,img)
    nii = template;
    nii.img = img;
    nii.hdr.dime.datatype = 64;
end

% return boolean value to check if the input name contains certain string
function bool = ContainName(name,string)
    bool= ~isempty(strfind(lower(name),string));
end