%% [localField,maskFinal] = BackgroundRemovalMacroIOWrapper(inputDir,outputDir,varargin)
%
% Input
% --------------
% inputDir              : input directory contains NIfTI (*totalfield* and *fieldmapsd*) files 
% outputDir             : output directory that stores the output (local field and final mask)
% varargin ('Name','Value' pair)
% ---------
% 'mask'                : mask file (in NIfTI format) full name 
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
%
% Output
% --------------
% localField            : local field (or tissue field) (in Hz)
% maskFinal             : final mask used for QSM
%
% Description: This is a wrapper of BackgroundRemovalMacro.m which has the following objeectives:
%               (1) matches the input format of qsm_hub.m
%               (2) save the results in NIfTI format
%
% Kwok-shing Chan @ DCCN
% k.chan@donders.ru.nl
% Date created: 17 April 2018
% Date last modified: 18 April 2018
%
%
function [localField,maskFinal] = BackgroundRemovalMacroIOWrapper(inputDir,outputDir,varargin)
%% add general Path
qsm_hub_AddMethodPath;

%% define variables
prefix = 'squirrel_';
gyro = 42.57747892;
% make sure the input only load once (first one)
isTotalFieldLoad = false;
isFieldmapSDLoad = false;

%% Check output directory exist or not
if exist(outputDir,'dir') ~= 7
    % if not then create the directory
    mkdir(outputDir);
end

%% Parse input argument using parse_varargin_QSMHub.m
[~,maskFullName,~,~,BFR,refine,BFR_tol,BFR_depth,BFR_peel,BFR_iteration,...
BFR_padSize,BFR_radius,BFR_alpha,BFR_threshold,~,~,~,~,~,~,~,~,~,~,~,~,~,~,~,~,~,~,~,~,~,~,~] = parse_varargin_QSMHub(varargin);

%% Read input
disp('Reading data...');
% look for nifti files 
inputNiftiList = dir([inputDir '/*.nii*']);
if ~isempty(inputNiftiList)
    % look for total field map and fieldmap SD (optional) NIfTI files
    for klist = 1:length(inputNiftiList)

        if ContainName(inputNiftiList(klist).name,'totalfield') && ~isTotalFieldLoad
            inputTotalFieldNifti = load_untouch_nii([inputDir filesep inputNiftiList(klist).name]);
            totalField = double(inputTotalFieldNifti.img);
            isTotalFieldLoad = true;
        end
        
        if ContainName(inputNiftiList(klist).name,'fieldmapsd') && ~isFieldmapSDLoad
            inputFieldMapSDNifti = load_untouch_nii([inputDir filesep inputNiftiList(klist).name]);
            fieldmapSD = double(inputFieldMapSDNifti.img);
            isFieldmapSDLoad = true;
        end
    end
    
    % if no files matched the name format then displays error message
    if ~isTotalFieldLoad
        error('No files loaded. Please make sure the input directory contains files with name *totalfield*');
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
        [B0,B0_dir,voxelSize,matrixSize,TE,delta_TE,CF]=SyntheticQSMHubHeader(inputTotalFieldNifti);
        
    end
    
    % if no fieldmapSD found then creates one with all voxels have the same value
    if ~isFieldmapSDLoad
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
    error('This standalone only reads NIfTI format input data (*.nii or *.nii.gz).');
end

% display some header info
disp('Basic DICOM information');
disp(['Voxel size(x,y,z mm^3) =  ' num2str(voxelSize(1)) 'x' num2str(voxelSize(2)) 'x' num2str(voxelSize(3))]);
disp(['matrix size(x,y,z) =  ' num2str(matrixSize(1)) 'x' num2str(matrixSize(2)) 'x' num2str(matrixSize(3))]);
disp(['B0 direction(x,y,z) =  ' num2str(B0_dir(:)')]);
disp(['Field strength(T) =  ' num2str(B0)]);

%% get brain mask
maskList = dir([inputDir '/*mask*']);

if ~isempty(maskFullName)
% first read mask if file is provided
    mask = load_nii_img_only(maskFullName) > 0;
    
elseif ~isempty(maskList) 
% read mask if input directory contains *mask*
    inputMaskNii = load_untouch_nii([inputDir filesep maskList(1).name]);
	mask = inputMaskNii.img > 0;
    
else
    % display error message if nothing is found
    error('No mask is found. Pleasee specific your mask file or put it inside the input directory.');
    
end

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