%% [chi,localField,totalField,fieldmapSD]=QSMHub(inputDir,outputDir,varargin)
%
% Input
% --------------
%
% Output
% --------------
%
% Description:
%
% Kwok-shing Chan @ DCCN
% k.chan@donders.ru.nl
% Date created: 14 September 2017
% Date last modified: 9 April 2018
%
%
function [chi,localField,totalField,fieldmapSD]=QSMHub(inputDir,outputDir,varargin)
%% add general Path
qsm_hub_AddMethodPath % qsm_hub_AddPath;

gyro = 42.57747892;

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
QSM_optimise,QSM_tol,QSM_maxiter,QSM_tol1,QSM_tol2,QSM_padsize,QSM_mu1,QSM_solver,QSM_constraint,...
exclude_threshold,QSM_radius,QSM_zeropad,QSM_wData,QSM_wGradient,QSM_lambdaCSF,QSM_isSMV,QSM_merit,isEddyCorrect] = parse_varargin_QSMHub(varargin);

%% Read input
disp('Reading data...');
% look for nifti files first
inputNiftiList = dir([inputDir '/*.nii*']);
if ~isempty(inputNiftiList)
    % look for magnitude and phase files
    for klist = 1:length(inputNiftiList)
        if contains(inputNiftiList(klist).name,'magn') && ~contains(inputNiftiList(klist).name,'brain')
            inputMagnNifti = load_untouch_nii([inputDir filesep inputNiftiList(klist).name]);
            magn = double(inputMagnNifti.img);
        end
        if contains(inputNiftiList(klist).name,'phase')
            inputPhaseNifti = load_untouch_nii([inputDir filesep inputNiftiList(klist).name]);
            fieldMap = double(inputPhaseNifti.img);
            % if input fieldmap is directly converted from nifti converter
            % then converts the fieldmap in rad and save to output dir
            if max(fieldMap(:))>1000
                fieldMap = DICOM2Phase(inputPhaseNifti);
                nii_fieldMap = make_nii_quick(inputPhaseNifti,fieldMap);
                nii_fieldMap.hdr.dime.scl_inter = 0;
                nii_fieldMap.hdr.dime.scl_slope = 1;
                save_untouch_nii(nii_fieldMap,[outputDir filesep 'qsmhub_phase.nii.gz']);
            end
        end
    end
    % create synthetic header in case no header is provided
    [B0,~,~,~,TE,delta_TE,CF]=SyntheticQSMHubHeader(magn);
    a=qGetR([0, inputMagnNifti.hdr.hist.quatern_b,inputMagnNifti.hdr.hist.quatern_c,inputMagnNifti.hdr.hist.quatern_d]);
    B0_dir = -a(3,:);
    voxelSize = inputMagnNifti.hdr.dime.pixdim(2:4);
    matrixSize = inputMagnNifti.hdr.dime.dim(2:4);
    % look for header file
    if ~isempty(dir([inputDir '/*header*']))
        headerList = dir([inputDir '/*header*']);
        load([inputDir filesep headerList(1).name]);
    else
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

    % save magnitude and phase images as nifti files
    disp('Saving DICOM data into NIfTI...');
    
    outputNiftiTemplate = make_nii(zeros(size(iField)), voxelSize);
    % make sure the class of output datatype is double
    outputNiftiTemplate.hdr.dime.datatype = 64;
    
    nii_fieldMap = make_nii_quick(outputNiftiTemplate,fieldMap);
    nii_magn = make_nii_quick(outputNiftiTemplate,magn);
    %% TODO: bug fix for load_nii and save_nii
    save_nii(nii_fieldMap,[outputDir filesep 'qsmhub_phase.nii.gz']);
    save_nii(nii_magn,[outputDir filesep 'qsmhub_magn.nii.gz']);
    save([outputDir filesep 'qsmhub_header.mat'],'voxelSize','matrixSize','CF','delta_TE',...
        'TE','B0_dir','B0');
    
    % remove the time dimension info
    outputNiftiTemplate = load_untouch_nii([outputDir filesep 'qsmhub_magn.nii.gz']);
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
maskList = dir([inputDir '/*mask*']);
% first read mask if file is provided
if ~isempty(maskFullName)
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

if isEddyCorrect
    disp('Correcting eddy current effect on bipolar readout data');
    imgCplx = BipolarEddyCorrect(magn.*exp(1i*fieldMap),mask,unwrap);
    fieldMap = angle(imgCplx);
    magn = abs(imgCplx);
    
    nii_fieldMap = make_nii_quick(outputNiftiTemplate,fieldMap);
    nii_fieldMap.hdr.dime.dim(5) = size(fieldMap,4);
    nii_magn = make_nii_quick(outputNiftiTemplate,magn); 
    nii_magn.hdr.dime.dim(5) = size(magn,4);
    save_untouch_nii(nii_fieldMap,[outputDir filesep 'qsmhub_phase_EC.nii.gz']);
    save_untouch_nii(nii_magn,[outputDir filesep 'qsmhub_magn_EC.nii.gz']);
end

disp('Calculating field map...');

% fix the output of field map in Hz
unit = 'Hz';

[totalField,fieldmapSD] = estimateTotalField(fieldMap,magn,matrixSize,voxelSize,...
                        'Unwrap',unwrap,'TE',TE,'B0',B0,'unit',unit,...
                        'Subsampling',subsampling,'mask',mask);
                                        
disp('Saving unwrapped field map...');

nii_totalField = make_nii_quick(outputNiftiTemplate,totalField);
nii_fieldmapSD = make_nii_quick(outputNiftiTemplate,fieldmapSD);
                    
save_untouch_nii(nii_totalField,[outputDir filesep 'qsmhub_totalField.nii.gz']);
save_untouch_nii(nii_fieldmapSD,[outputDir filesep 'qsmhub_fieldMapSD.nii.gz']);

maskReliable = fieldmapSD < exclude_threshold;
mask = and(mask,maskReliable);

%% Background field removal
qsm_hub_AddMethodPath(BFR);
disp('Recovering local field...');

localField = BackgroundRemovalMacro(totalField,mask,matrixSize,voxelSize,...
      'method',BFR,'refine',refine,'tol',BFR_tol,'depth',BFR_depth,...
      'peel',BFR_peel,'b0dir',B0_dir,'iteration',BFR_iteration,'padsize',...
      BFR_padSize,'noisestd',fieldmapSD,'radius',BFR_radius,'alpha',BFR_alpha,...
      'threshold',BFR_threshold); 
  
maskFinal = localField ~=0;
  
disp('Saving local field map...');

nii_localField = make_nii_quick(outputNiftiTemplate,localField);
nii_maskFinal = make_nii_quick(outputNiftiTemplate,maskFinal);
                    
save_untouch_nii(nii_localField,[outputDir filesep 'qsmhub_localField.nii.gz']);
save_untouch_nii(nii_maskFinal,[outputDir filesep 'qsmhub_finalMask.nii.gz']);

% % make sure the input of QSM is in rad
% localField = localField;
            
%% create weight map
wmap = fieldmapSD./norm(fieldmapSD(maskFinal==1));    

%% qsm
qsm_hub_AddMethodPath(QSM_method);
disp('Computing QSM...');

switch lower(QSM_method)
    case 'closedforml2'
    case 'ilsqr'
    case 'stisuiteilsqr'
    case 'fansi'
        % FANSI parameter is for ppm
        localField = localField/(B0*gyro);
    case 'ssvsharp'
    case 'star'
        localField = localField*2*pi;
    case 'medi_l1'
end

chi = qsmMacro(localField,maskFinal,matrixSize,voxelSize,...
      'method',QSM_method,'threshold',QSM_threshold,'lambda',QSM_lambda,...
      'optimise',QSM_optimise,'tol',QSM_tol,'iteration',QSM_maxiter,'weight',wmap,...
      'b0dir',B0_dir,'tol_step1',QSM_tol1,'tol_step2',QSM_tol2,'TE',delta_TE,'B0',B0,...
      'padsize',QSM_padsize,'mu',QSM_mu1,QSM_solver,QSM_constraint,...
      'noisestd',fieldmapSD,'magnitude',sqrt(sum(magn.^2,4)),'data_weighting',QSM_wData,...
      'gradient_weighting',QSM_wGradient,'merit',QSM_merit,'smv',QSM_isSMV,'zeropad',QSM_zeropad,...
      'lambda_CSF',QSM_lambdaCSF,'CF',CF,'radius',QSM_radius);

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
        chi = chi/(B0*gyro);
end
  
disp('Saving susceptibility map...');

nii_chi = make_nii_quick(outputNiftiTemplate,chi);

save_untouch_nii(nii_chi,[outputDir filesep 'qsmhub_QSM.nii.gz']);

disp('Done!');
          
end

% handy function to save result to nifti format
function nii = make_nii_quick(template,img)
    nii = template;
    nii.img = img;
    nii.hdr.dime.datatype = 64;
end