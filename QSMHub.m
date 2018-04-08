%% function [chi,localField,totalField,fieldmapSD]=QSMHub(inputDir,outputDir,varargin)
%
% Usage: [chi,localField,totalField,fieldmapSD]=QSMHub(inputDir,outputDir,varargin)
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
% Date last modified: 29 September 2017
%
%
function [chi,localField,totalField,fieldmapSD]=QSMHub(inputDir,outputDir,varargin)
%% add Path
qsm_hub_AddPath;
%% Check output directory exist or not
if exist(outputDir,'dir') ~= 7
    mkdir(outputDir);
end

%% Parse input argument
[isBET,maskFullName,unwrap,unit,subsampling,BFR,refine,BFR_tol,BFR_depth,BFR_peel,BFR_iteration,...
BFR_CGdefault,BFR_radius,BFR_alpha,BFR_threshold,QSM_method,QSM_threshold,QSM_lambda,...
QSM_optimise,QSM_tol,QSM_maxiter,QSM_tol1,QSM_tol2,QSM_padsize,QSM_mu1,QSM_solver,QSM_constraint] = parse_varargin_QSMHub(varargin);

%% Read input
disp('Reading data...');
% look for nifti files first
inputNiiList = dir([inputDir '/*.nii*']);
if ~isempty(inputNiiList)
    % look for magnitude and phase files
    for klist = 1:length(inputNiiList)
        if contains(inputNiiList(klist).name,'magn')
            inputMagnNii = load_nii([inputDir filesep inputNiiList(klist).name]);
            magn = inputMagnNii.img;
        end
        if contains(inputNiiList(klist).name,'phase')
            inputPhaseNii = load_nii([inputDir filesep inputNiiList(klist).name]);
            fieldMap = inputPhaseNii.img;
            voxelSize = inputPhaseNii.hdr.dime.pixdim(2:4);
            matrixSize = inputPhaseNii.hdr.dime.dim(2:4);
            B0_dir = [0 0 1];
            B0 = 3;
            TE = linspace(1e-3,30e-3,inputPhaseNii.hdr.dime.dim(5));
        end
    end
    % look for header file
    if ~isempty(dir([inputDir '/*header*']))
        headerList = dir([inputDir '/*header*']);
        load([inputDir filesep headerList(1).name]);
    end
else
    % if no nifti files than check for DICOM files
    [iField,voxelSize,matrixSize,CF,delta_TE,TE,B0_dir]=Read_Siemens_DICOM_old(inputDir);
    gyro = 42.57747892;
    B0 = CF/(gyro*1e6);

    fieldMap = angle(iField);
    magn = abs(iField);

    % save magnitude and phase images as nifti files
    disp('Saving DICOM data...');
    
    nii_fieldMap = make_nii(fieldMap, voxelSize);
    nii_magn = make_nii(magn, voxelSize);
    save_nii(nii_fieldMap,[outputDir filesep 'qsmhub_phase.nii.gz']);
    save_nii(nii_magn,[outputDir filesep 'qsmhub_magn.nii.gz']);
    save([outputDir filesep 'qsmhub_header.mat'],'voxelSize','matrixSize','CF','delta_TE',...
        'TE','B0_dir','B0');
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
    inputMaskNii = load_nii([inputDir filesep maskList(1).name]);
	mask = inputMaskNii.img > 0;
end
    
% is BET is checked or no mask is found, run FSL's bet
if isempty(mask) || isBET
    disp('Performing FSL BET...');
    nii_temp = make_nii(magn(:,:,:,1), voxelSize);
    tempDir = [outputDir filesep 'qsmhub_temp.nii.gz'];
    brianDir = [outputDir filesep 'temp_brain.nii.gz'];
    
    save_nii(nii_temp,tempDir);
    
    system(['bet ' tempDir ' ' brianDir ' -R']);
    try 
        mask = load_nii_img_only(brianDir) > 0;
    catch
        setenv( 'FSLDIR', '/usr/local/fsl');
        fsldir = getenv('FSLDIR');
        fsldirmpath = sprintf('%s/etc/matlab',fsldir);
        path(path, fsldirmpath);
        oldPATH = getenv('PATH');
        setenv('PATH',[oldPATH ':' fsldir '/bin']);
        call_fsl(['bet ' tempDir ' ' brianDir ' -R']);
        mask = load_nii_img_only(brianDir) > 0;
    end

    system(['rm ' tempDir]);
end

%% total field and Laplacian phase unwrap
disp('Calculating field map...');

[totalField,fieldmapSD] = estimateTotalField(fieldMap,magn,matrixSize,voxelSize,...
                        'Unwarp',unwrap,'TE',TE,'B0',B0,'unit',unit,...
                        'Subsampling',subsampling,'mask',mask);
                                        
disp('Saving unwrapped field map...');

nii_totalField = make_nii(totalField, voxelSize);
nii_fieldmapSD = make_nii(fieldmapSD, voxelSize);
                    
save_nii(nii_totalField,[outputDir filesep 'qsmhub_totalField.nii.gz']);
save_nii(nii_fieldmapSD,[outputDir filesep 'qsmhub_fieldMapSD.nii.gz']);

%% Background field removal
disp('Recovering local field...');

localField = BackgroundRemovalMacro(totalField,mask,matrixSize,voxelSize,...
      'method',BFR,'refine',refine,'tol',BFR_tol,'depth',BFR_depth,...
      'peel',BFR_peel,'b0dir',B0_dir,'iteration',BFR_iteration,'CGsolver',...
      BFR_CGdefault,'noisestd',fieldmapSD,'radius',BFR_radius,'alpha',BFR_alpha,...
      'threshold',BFR_threshold); 
  
maskFinal = localField ~=0;
  
disp('Saving local field map...');

nii_localField = make_nii(localField, voxelSize);
nii_maskFinal = make_nii(uint8(maskFinal), voxelSize);
                    
save_nii(nii_localField,[outputDir filesep 'qsmhub_localField.nii.gz']);
save_nii(nii_maskFinal,[outputDir filesep 'qsmhub_finalMask.nii.gz']);
            
%% create weight map
wmap = fieldmapSD./norm(fieldmapSD(maskFinal==1));    

%% qsm
disp('Computing QSM...');

chi = qsmMacro(localField,maskFinal,matrixSize,voxelSize,...
      'method',QSM_method,'threshold',QSM_threshold,'lambda',QSM_lambda,...
      'optimise',QSM_optimise,'tol',QSM_tol,'iteration',QSM_maxiter,'weight',wmap,...
      'b0dir',B0_dir,'tol_step1',QSM_tol1,'tol_step2',QSM_tol2,'TE',1,'B0',B0,...
      'padsize',QSM_padsize,'mu',QSM_mu1,QSM_solver,QSM_constraint);
  
disp('Saving susceptibility map...');

nii_chi = make_nii(chi, voxelSize);

save_nii(nii_chi,[outputDir filesep 'qsmhub_QSM.nii.gz']);
          
end