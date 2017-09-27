%% function [chi,localField,totalField,fieldmapSD]=QSMHub(inputDir,outputDir,varargin)
%
% Usage:
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
% Date last modified:
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
[isBET,mask,unwrap,unit,subsampling,BFR,refine,BFR_tol,BFR_depth,BFR_peel,BFR_iteration,...
BFR_CGdefault,BFR_radius,BFR_alpha,BFR_threshold,QSM_method,QSM_threshold,QSM_lambda,...
QSM_optimise,QSM_tol,QSM_maxiter,QSM_tol1,QSM_tol2,QSM_padsize,QSM_mu1,QSM_solver,QSM_constraint] = parse_varargin_QSMHub(varargin);

%% Read DICOM
disp('Reading data...');
[iField,voxelSize,matrixSize,CF,delta_TE,TE,B0_dir]=Read_Siemens_DICOM_old(inputDir);
gyro = 42.58;
B0 = CF/(gyro*1e6);
    
fieldMap = angle(iField);
magn = abs(iField);

disp('Saving DICOM data...');

nii_fieldMap = make_nii(fieldMap, voxelSize);
nii_magn = make_nii(magn, voxelSize);

save_nii(nii_fieldMap,[outputDir filesep 'qsmhub_phase.nii.gz']);
save_nii(nii_magn,[outputDir filesep 'qsmhub_magn.nii.gz']);
save([outputDir filesep 'qsmhub_header.mat'],'voxelSize','matrixSize','CF','delta_TE',...
    'TE','B0_dir','B0');
%% FSL's BET
if isempty(mask) || isBET
    disp('Performing FSL BET...');
    nii_temp = make_nii(magn(:,:,:,1), voxelSize);
    tempDir = [outputDir filesep 'qsmhub_temp.nii.gz'];
    brianDir = [outputDir filesep 'temp_brain.nii.gz'];
    save_nii(nii_temp,tempDir);
    try 
        system(['bet ' tempDir ' ' brianDir ' -R']);
    catch
        setenv('PATH', [getenv('PATH'),':','/usr/local/fsl/bin']);
        system(['bet ' tempDir brianDir ' -R']);
    end
    mask = load_nii_img_only(brianDir) > 0;
    system(['rm ' tempDir]);
end

%% total field and Laplacian phase unwrap
disp('Calculating field map...');

[totalField,fieldmapSD] = estimateTotalField(fieldMap,magn,matrixSize,voxelSize,...
                        'Unwarp',unwrap,'TE',TE,'B0',B0,'unit',unit,...
                        'Subsampling',subsampling,'mask',mask);
                                        
disp('Saving total field map...');

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
chi = qsmMacro(localField,maskFinal,matrixSize,voxelSize,...
      'method',QSM_method,'threshold',QSM_threshold,'lambda',QSM_lambda,...
      'optimise',QSM_optimise,'tol',QSM_tol,'iteration',QSM_maxiter,'weight',wmap,...
      'b0dir',B0_dir,'tol_step1',QSM_tol1,'tol_step2',QSM_tol2,'TE',1,'B0',B0,...
      'padsize',QSM_padsize,'mu',QSM_mu1,QSM_solver,QSM_constraint);
  
nii_chi = make_nii(chi, voxelSize);

save_nii(nii_chi,[outputDir filesep 'qsmhub_QSM.nii.gz']);
          
end