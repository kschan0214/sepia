%% function output = function_name(input)
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
% Date created: 
% Date last modified:
%
%
function [totalField,fieldmapSD]=UnwrapPhaseIOMacro(inputDir,outputDir,varargin)
%% add general Path
qsm_hub_AddMethodPath % qsm_hub_AddPath;

gyro = 42.57747892;

%% Check output directory exist or not
if exist(outputDir,'dir') ~= 7
    % if not then create the directory
    mkdir(outputDir);
end

%% Parse input argument
[isBET,maskFullName,unwrap,subsampling,~,~,~,~,~,~,...
~,~,~,~,~,~,~,~,~,~,~,~,~,~,~,~,exclude_threshold,~,~,~,~,~,~,~,isEddyCorrect] = parse_varargin_QSMHub(varargin);

%% Read input
disp('Reading data...');
% look for nifti files first
inputNiftiList = dir([inputDir '/*.nii*']);
if ~isempty(inputNiftiList)
    % look for magnitude and phase files
    for klist = 1:length(inputNiftiList)
        if contains(lower(inputNiftiList(klist).name),'magn') && ~contains(lower(inputNiftiList(klist).name),'brain')
            inputMagnNifti = load_untouch_nii([inputDir filesep inputNiftiList(klist).name]);
            magn = double(inputMagnNifti.img);
        end
        if contains(lower(inputNiftiList(klist).name),'phase')
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

    mask = BET(magn(:,:,:,1),matrixSize,voxelSize);

end
%% total field and phase unwrap
% add 'unwrap' method PATH
qsm_hub_AddMethodPath(unwrap);

% Eddy current correction for bipolar readout
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
                    
maskReliable = fieldmapSD < exclude_threshold;
mask = and(mask,maskReliable);
                                        
disp('Saving unwrapped field map...');

nii_totalField = make_nii_quick(outputNiftiTemplate,totalField);
nii_fieldmapSD = make_nii_quick(outputNiftiTemplate,fieldmapSD);
nii_newMask = make_nii_quick(outputNiftiTemplate,mask);
                    
save_untouch_nii(nii_totalField,[outputDir filesep 'qsmhub_totalField.nii.gz']);
save_untouch_nii(nii_fieldmapSD,[outputDir filesep 'qsmhub_fieldMapSD.nii.gz']);
save_untouch_nii(nii_newMask,[outputDir filesep 'qsmhub_mask_new.nii.gz']);

disp('Done!');

end

% handy function to save result to nifti format
function nii = make_nii_quick(template,img)
    nii = template;
    nii.img = img;
    nii.hdr.dime.datatype = 64;
end