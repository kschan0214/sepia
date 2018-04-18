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
function [localField,maskFinal] = BackgroundRemovalMacroIOWrapper(inputDir,outputDir,varargin)
%% add general Path
qsm_hub_AddMethodPath % qsm_hub_AddPath;

gyro = 42.57747892;

%% Check output directory exist or not
if exist(outputDir,'dir') ~= 7
    % if not then create the directory
    mkdir(outputDir);
end

%% Parse input argument
[~,maskFullName,~,~,BFR,refine,BFR_tol,BFR_depth,BFR_peel,BFR_iteration,...
BFR_padSize,BFR_radius,BFR_alpha,BFR_threshold,~,~,~,~,~,~,~,~,~,~,~,~,~,~,~,~,~,~,~,~,~] = parse_varargin_QSMHub(varargin);

%% Read input
disp('Reading data...');
% look for nifti files first
inputNiftiList = dir([inputDir '/*.nii*']);
if ~isempty(inputNiftiList)
    % look for total field map and fieldmap SD files
    fieldmapSD = [];
    for klist = 1:length(inputNiftiList)
        if contains(lower(inputNiftiList(klist).name),'totalfield') 
            inputTotalFieldNifti = load_untouch_nii([inputDir filesep inputNiftiList(klist).name]);
            totalField = double(inputTotalFieldNifti.img);
        end
        if contains(lower(inputNiftiList(klist).name),'fieldmapsd')
            inputFieldMapSDNifti = load_untouch_nii([inputDir filesep inputNiftiList(klist).name]);
            fieldmapSD = double(inputFieldMapSDNifti.img);
        end
    end
    % create synthetic header in case no header is provided
    B0=3; % T
%     [B0,~,~,~,TE,delta_TE,CF]=SyntheticQSMHubHeader(totalField);
    a=qGetR([0, inputTotalFieldNifti.hdr.hist.quatern_b,inputTotalFieldNifti.hdr.hist.quatern_c,inputTotalFieldNifti.hdr.hist.quatern_d]);
    B0_dir = -a(3,:);
    voxelSize = inputTotalFieldNifti.hdr.dime.pixdim(2:4);
    matrixSize = inputTotalFieldNifti.hdr.dime.dim(2:4);
    
    % look for header file
    if ~isempty(dir([inputDir '/*header*']))
        headerList = dir([inputDir '/*header*']);
        load([inputDir filesep headerList(1).name]);
    else
        disp('No QSMHub header file detected. Using synthetic header parameters...');
    end
    if isempty(fieldmapSD)
        fieldmapSD = ones(matrixSize) * 0.01;
    end
    % store the header the NIfTI files, all following results will have
    % the same header
    outputNiftiTemplate = inputTotalFieldNifti;
    % make sure the class of output datatype is double
    outputNiftiTemplate.hdr.dime.datatype = 64;
    % remove the time dimension info
    outputNiftiTemplate.hdr.dime.dim(5) = 1;
    
else
    error('This standalone only reads NIfTI format input data (nii or nii.gz).');
end

% display some header info
disp('Basic DICOM information');
disp(['Voxel size(x,y,z mm^3) =  ' num2str(voxelSize(1)) 'x' num2str(voxelSize(2)) 'x' num2str(voxelSize(3))]);
disp(['matrix size(x,y,z) =  ' num2str(matrixSize(1)) 'x' num2str(matrixSize(2)) 'x' num2str(matrixSize(3))]);
disp(['B0 direction(x,y,z) =  ' num2str(B0_dir(:)')]);
disp(['Field strength(T) =  ' num2str(B0)]);

%% get brain mask
maskList = dir([inputDir '/*mask*']);
% first read mask if file is provided
if ~isempty(maskFullName)
    mask = load_nii_img_only(maskFullName) > 0;
elseif ~isempty(maskList) 
    % read mask if input directory contains 'mask'
    inputMaskNii = load_untouch_nii([inputDir filesep maskList(1).name]);
	mask = inputMaskNii.img > 0;
else
    error('No mask is found. Pleasee specific your mask file or put it inside the input directory.');
end

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
save_untouch_nii(nii_maskFinal,[outputDir filesep 'qsmhub_mask_final.nii.gz']);

disp('Done!');

end

% handy function to save result to nifti format
function nii = make_nii_quick(template,img)
    nii = template;
    nii.img = img;
    nii.hdr.dime.datatype = 64;
end