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
% Date last modified: 18 April 2018
%
%
function [chi,localField,totalField,fieldmapSD]=QSMHub(inputDir,outputDir,varargin)
%% add general Path
qsm_hub_AddMethodPath % qsm_hub_AddPath;

%% define variables
prefix = 'squirrel_';
gyro = 42.57747892;
% make sure the input only load once (first one)
isMagnLoad = false;
isPhaseLoad = false;
maskCSF = [];

%% Check output directory exist or not
if exist(outputDir,'dir') ~= 7
    % if not then create the directory
    mkdir(outputDir);
end

%% Parse input argument
% [isBET,maskFullName,unwrap,unit,subsampling,BFR,refine,BFR_tol,BFR_depth,BFR_peel,BFR_iteration,...
% BFR_CGdefault,BFR_radius,BFR_alpha,BFR_threshold,QSM_method,QSM_threshold,QSM_lambda,...
% QSM_optimise,QSM_tol,QSM_maxiter,QSM_tol1,QSM_tol2,QSM_padsize,QSM_mu1,QSM_solver,QSM_constraint,exclude_threshold] = parse_varargin_QSMHub(varargin);
[isBET,maskFullName,unwrap,subsampling,BFR,refine,BFR_tol,BFR_depth,BFR_peel,BFR_iteration,...
BFR_padSize,BFR_radius,BFR_alpha,BFR_threshold,QSM_method,QSM_threshold,QSM_lambda,...
QSM_optimise,QSM_tol,QSM_maxiter,QSM_tol1,QSM_tol2,QSM_padsize,QSM_mu1,QSM_mu2,QSM_solver,QSM_constraint,...
exclude_threshold,QSM_radius,QSM_zeropad,QSM_wData,QSM_wGradient,QSM_isLambdaCSF,QSM_lambdaCSF,QSM_isSMV,QSM_merit,isEddyCorrect] = parse_varargin_QSMHub(varargin);

%% Read input
disp('Reading data...');
% look for nifti files first
inputNiftiList = dir([inputDir '/*.nii*']);
if ~isempty(inputNiftiList)
    % look for magnitude and phase files
    for klist = 1:length(inputNiftiList)
        if ContainName(lower(inputNiftiList(klist).name),'magn') && ~ContainName(lower(inputNiftiList(klist).name),'brain') && ~isMagnLoad
            inputMagnNifti = load_untouch_nii([inputDir filesep inputNiftiList(klist).name]);
            magn = double(inputMagnNifti.img);
            isMagnLoad = true;
        end
        if ContainName(lower(inputNiftiList(klist).name),'phase') && ~isPhaseLoad
            inputPhaseNifti = load_untouch_nii([inputDir filesep inputNiftiList(klist).name]);
            fieldMap = double(inputPhaseNifti.img);
            
            % if input fieldmap is directly converted from nifti converter
            % then converts the fieldmap in rad and save to output dir
            if max(fieldMap(:))>1000
                fieldMap = DICOM2Phase(inputPhaseNifti);
                nii_fieldMap = make_nii_quick(inputPhaseNifti,fieldMap);
                nii_fieldMap.hdr.dime.scl_inter = 0;
                nii_fieldMap.hdr.dime.scl_slope = 1;
                save_untouch_nii(nii_fieldMap,[outputDir filesep prefix 'phase.nii.gz']);
            end
            isPhaseLoad = true;
        end
    end
    
    % if no files matched the name format then displays error message
    if ~isMagnLoad
        error('No files loaded. Please make sure the input directory contains files with name *magn*');
    end
    % if no files matched the name format then displays error message
    if ~isPhaseLoad
        error('No files loaded. Please make sure the input directory contains files with name *phase*');
    end
    
    % look for header file
    if ~isempty(dir([inputDir '/*header*']))
        disp('Reading header for qsm_hub...');
        
        % load header
        headerList = dir([inputDir '/*header*']);
        load([inputDir filesep headerList(1).name]);
    else
        disp('No header for qsm_hub is found. Creating synthetic header based on NIfTI header...');
        
        % create synthetic header in case no qsm_hub's header is found
        [B0,B0_dir,voxelSize,matrixSize,TE,delta_TE,CF]=SyntheticQSMHubHeader(inputMagnNifti);
        
        % if no header file then save the synthetic header in output dir
        save([outputDir filesep 'SyntheticQSMhub_header'],'voxelSize','matrixSize','CF','delta_TE',...
        'TE','B0_dir','B0');
    end
    
    % store the header the NIfTI files, all following results will have
    % the same header
    outputNiftiTemplate = inputMagnNifti;
    % make sure the class of output datatype is double
    outputNiftiTemplate.hdr.dime.datatype = 64;
    % remove the time dimension info
    outputNiftiTemplate.hdr.dime.dim(5) = 1;
    
else
    % if no nifti file then check for DICOM files
%     [iField,voxelSize,matrixSize,CF,delta_TE,TE,B0_dir]=Read_Siemens_DICOM_old(inputDir);
    [iField,voxelSize,matrixSize,CF,delta_TE,TE,B0_dir]=Read_DICOM(inputDir);
    
    B0 = CF/(gyro*1e6);

    fieldMap = angle(iField);
    magn = abs(iField);
    
    isMagnLoad = true;
    isPhaseLoad = true;

    % save magnitude and phase images as nifti files
    disp('Saving DICOM data into NIfTI...');
    
    outputNiftiTemplate = make_nii(zeros(size(iField)), voxelSize);
    % make sure the class of output datatype is double
    outputNiftiTemplate.hdr.dime.datatype = 64;
    
    nii_fieldMap = make_nii_quick(outputNiftiTemplate,fieldMap);
    nii_magn = make_nii_quick(outputNiftiTemplate,magn);
    
    % save magnitude and phase data as NIfTI_GZ format
    save_nii(nii_fieldMap,[outputDir filesep prefix 'phase.nii.gz']);
    save_nii(nii_magn,[outputDir filesep prefix 'magn.nii.gz']);
    % save important header in .mat format
    save([outputDir filesep 'qsmhub_header.mat'],'voxelSize','matrixSize','CF','delta_TE',...
        'TE','B0_dir','B0');
    
    % reload the NIfTI template so that later can use save_untouch_nii for
    % all results
    outputNiftiTemplate = load_untouch_nii([outputDir filesep prefix 'magn.nii.gz']);
    % remove the time dimension info
    outputNiftiTemplate.hdr.dime.dim(5) = 1;
    
end
% display some header info
disp('Basic DICOM information');
disp(['Voxel size(x,y,z mm^3) =  ' num2str(voxelSize(1)) 'x' num2str(voxelSize(2)) 'x' num2str(voxelSize(3))]);
disp(['matrix size(x,y,z) =  ' num2str(matrixSize(1)) 'x' num2str(matrixSize(2)) 'x' num2str(matrixSize(3))]);
disp(['B0 direction(x,y,z) =  ' num2str(B0_dir(:)')]);
disp(['Field strength(T) =  ' num2str(B0)]);

%% get brain mask
mask = [];
maskList = dir([inputDir '/*mask*nii*']);
if ~isempty(maskFullName)
% first read mask if file is provided
    mask = load_nii_img_only(maskFullName) > 0;
elseif ~isempty(maskList) 
% read mask if input directory contains 'mask'
    inputMaskNii = load_untouch_nii([inputDir filesep maskList(1).name]);
	mask = inputMaskNii.img > 0;
end
    
% if BET is checked or no mask is found, run FSL's bet
if isempty(mask) || isBET
    qsm_hub_AddMethodPath('bet');
    disp('Performing FSL BET...');
    
%     tempDir = [outputDir filesep 'qsmhub_temp.nii.gz'];
%     brianDir = [outputDir filesep 'temp_brain.nii.gz'];
%     
%     % save the 1st echo for bet
%     nii_temp = make_nii_quick(outputNiftiTemplate,magn(:,:,:,1));
%     
%     save_untouch_nii(nii_temp,tempDir);
%     
%     % run bet here
%     system(['bet ' tempDir ' ' brianDir ' -R']);
%     try 
%         mask = load_nii_img_only(brianDir) > 0;
%     catch
%         setenv( 'FSLDIR', '/usr/local/fsl');
%         fsldir = getenv('FSLDIR');
%         fsldirmpath = sprintf('%s/etc/matlab',fsldir);
%         path(path, fsldirmpath);
%         oldPATH = getenv('PATH');
%         setenv('PATH',[oldPATH ':' fsldir '/bin']);
%         call_fsl(['bet ' tempDir ' ' brianDir ' -R']);
%         mask = load_nii_img_only(brianDir) > 0;
%     end
% 
%     system(['rm ' tempDir]);

    mask = BET(magn(:,:,:,1),matrixSize,voxelSize);

end

%% total field and Laplacian phase unwrap
% add 'unwrap' method PATH
qsm_hub_AddMethodPath(unwrap);

% Eddy current correction for bipolar readout
if isEddyCorrect
    disp('Correcting eddy current effect on bipolar readout data');
    
    % BipolarEddyCorrect requries complex-valeud input
    imgCplx = BipolarEddyCorrect(magn.*exp(1i*fieldMap),mask,unwrap);
    fieldMap = angle(imgCplx);
    magn = abs(imgCplx);
    
    % save the eddy current corrected output
    nii_fieldMap = make_nii_quick(outputNiftiTemplate,fieldMap);
    nii_fieldMap.hdr.dime.dim(5) = size(fieldMap,4);
    nii_magn = make_nii_quick(outputNiftiTemplate,magn); 
    nii_magn.hdr.dime.dim(5) = size(magn,4);
    save_untouch_nii(nii_fieldMap,[outputDir filesep prefix 'phase_EC.nii.gz']);
    save_untouch_nii(nii_magn,[outputDir filesep prefix 'magn_EC.nii.gz']);
end

disp('Calculating field map...');

% fix the output of field map in Hz
unit = 'Hz';

% core of phase unwrapping
try 
    [totalField,fieldmapSD] = estimateTotalField(fieldMap,magn,matrixSize,voxelSize,...
                            'Unwrap',unwrap,'TE',TE,'B0',B0,'unit',unit,...
                            'Subsampling',subsampling,'mask',mask);
catch
    % if the selected method is not working then do Laplacian 
    disp('The selected method is not supported in this system. Using Laplacian algorithm for phase unwrapping...')
    qsm_hub_AddMethodPath('laplacian');
    
    [totalField,fieldmapSD] = estimateTotalField(fieldMap,magn,matrixSize,voxelSize,...
                            'Unwrap','laplacian','TE',TE,'B0',B0,'unit',unit,...
                            'Subsampling',subsampling,'mask',mask);
end
             
% save the output                           
disp('Saving unwrapped field map...');

nii_totalField = make_nii_quick(outputNiftiTemplate,totalField);
nii_fieldmapSD = make_nii_quick(outputNiftiTemplate,fieldmapSD);
                    
save_untouch_nii(nii_totalField,[outputDir filesep prefix 'totalField.nii.gz']);
save_untouch_nii(nii_fieldmapSD,[outputDir filesep prefix 'fieldMapSD.nii.gz']);

maskReliable = fieldmapSD < exclude_threshold;
mask = and(mask,maskReliable);

%% Background field removal
qsm_hub_AddMethodPath(BFR);
disp('Recovering local field...');

% core of background field removal
localField = BackgroundRemovalMacro(totalField,mask,matrixSize,voxelSize,...
      'method',BFR,'refine',refine,'tol',BFR_tol,'depth',BFR_depth,...
      'peel',BFR_peel,'b0dir',B0_dir,'iteration',BFR_iteration,'padsize',...
      BFR_padSize,'noisestd',fieldmapSD,'radius',BFR_radius,'alpha',BFR_alpha,...
      'threshold',BFR_threshold); 
  
% generate new mask based on backgroudn field removal result
maskFinal = localField ~=0;
  
% save results
disp('Saving local field map...');

nii_localField = make_nii_quick(outputNiftiTemplate,localField);
nii_maskFinal = make_nii_quick(outputNiftiTemplate,maskFinal);
                    
save_untouch_nii(nii_localField,[outputDir filesep prefix 'localField.nii.gz']);
save_untouch_nii(nii_maskFinal,[outputDir filesep prefix 'mask_final.nii.gz']);
            
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
    case 'fansi'
        % FANSI works better with ppm
        localField = localField/(B0*gyro);
    case 'ssvsharp'
    case 'star'
        % star work better with radHz
        localField = localField*2*pi;
    case 'medi_l1'
        % zero reference using CSF requires CSF mask
        if QSM_isLambdaCSF && isMagnLoad
            disp('Extracting CSF mask....');
            
            % R2* mapping
            r2s = arlo(TE,magn);
            maskCSF = extract_CSF(r2s,maskFinal,voxelSize)>0;
        end
        
        % MEDI works better with rad
        localField = localField*2*pi*delta_TE;
end

% core of QSM
chi = qsmMacro(localField,maskFinal,matrixSize,voxelSize,...
      'method',QSM_method,'threshold',QSM_threshold,'lambda',QSM_lambda,...
      'optimise',QSM_optimise,'tol',QSM_tol,'iteration',QSM_maxiter,'weight',wmap,...
      'b0dir',B0_dir,'tol_step1',QSM_tol1,'tol_step2',QSM_tol2,'TE',delta_TE,'B0',B0,...
      'padsize',QSM_padsize,'mu',QSM_mu1,'mu2',QSM_mu2,QSM_solver,QSM_constraint,...
      'noisestd',fieldmapSD,'magnitude',sqrt(sum(magn.^2,4)),'data_weighting',QSM_wData,...
      'gradient_weighting',QSM_wGradient,'merit',QSM_merit,'smv',QSM_isSMV,'zeropad',QSM_zeropad,...
      'lambda_CSF',QSM_lambdaCSF,'CF',CF,'radius',QSM_radius,'Mask_CSF',maskCSF);

% convert the susceptibility map into ppm
switch lower(QSM_method)
    case 'tkd'
        chi = chi/(B0*gyro);
    case 'closedforml2'
        chi = chi/(B0*gyro);
    case 'ilsqr'
        chi = chi/(B0*gyro);
    case 'stisuiteilsqr'
        chi = chi/(B0*gyro);
    case 'fansi'
    case 'ssvsharp'
        chi = chi/(B0*gyro);
    case 'star'
        chi = chi/(2*pi*B0*gyro);
    case 'medi_l1'
        chi = chi/(2*pi*B0*gyro*delta_TE);
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