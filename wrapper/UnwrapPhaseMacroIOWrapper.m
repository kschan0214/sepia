%% [totalField,fieldmapSD]=UnwrapPhaseMacroIOWrapper(inputDir,outputDir,varargin)
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
%
% Output
% --------------
% totalField            : unwrapped field map (in Hz)
% fieldmapSD            : relative standard deviation of field map,
%                         esimated using Eq. 11 of Robinson et al. NMR Biomed 2017 (doi:10.1002/nbm.3601)
% mask                  : new mask with thresholded unreliable voxels
%
% Description: This is a wrapper of estimateTotalField.m which has the following objeectives:
%               (1) matches the input format of qsm_hub.m
%               (2) save the results in NIfTI format
%
% Kwok-shing Chan @ DCCN
% k.chan@donders.ru.nl
% Date created: 17 April 2018
% Date last modified: 18 April 2018
%
%
function [totalField,fieldmapSD,mask]=UnwrapPhaseMacroIOWrapper(inputDir,outputDir,varargin)
%% add general Path
qsm_hub_AddMethodPath;

%% define variables
prefix = 'squirrel_';
gyro = 42.57747892;
% make sure the input only load once (first one)
isMagnLoad = false;
isPhaseLoad = false;

%% Check output directory exist or not
if exist(outputDir,'dir') ~= 7
    % if not then create the directory
    mkdir(outputDir);
end

%% Parse input argument using parse_varargin_QSMHub.m
[isBET,maskFullName,phaseCombMethod,unwrap,subsampling,~,~,~,~,~,~,...
~,~,~,~,~,~,~,~,~,~,~,~,~,~,~,~,~,exclude_threshold,~,~,~,~,~,~,~,~,isEddyCorrect] = parse_varargin_QSMHub(varargin);

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
    
    % store the header of the NIfTI files as template, all following results will have
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

    % after testing with a couple of dataset it seems to me that the field
    % is inverted with DICOM input, so apply conjugate here
    fieldMap = angle(conj(iField));
    magn = abs(iField);
    
    isMagnLoad = true;
    isPhaseLoad = true;

    % save magnitude and phase images as nifti files
    disp('Saving DICOM data into NIfTI...');
    
    % create NIfTI header template
    outputNiftiTemplate = make_nii(zeros(size(iField)), voxelSize);
    % make sure the class of output datatype is double
    outputNiftiTemplate.hdr.dime.datatype = 64;
    
    % create NIfTI structure for magnitude and phase data
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
maskList = dir([inputDir '/*mask*']);

if ~isempty(maskFullName)
% first read mask if file is provided
    mask = load_nii_img_only(maskFullName) > 0;
    
elseif ~isempty(maskList) 
% read mask if input directory has NIfTi file contains *mask* in name
    inputMaskNii = load_untouch_nii([inputDir filesep maskList(1).name]);
	mask = inputMaskNii.img > 0;
    
end
    
% if BET is checked or no mask is found, run FSL's bet
if isempty(mask) || isBET
    qsm_hub_AddMethodPath('bet');
    disp('Performing FSL BET...');
    
    % this is the BET functino provided with MEDI toolbox
    mask = BET(magn(:,:,:,1),matrixSize,voxelSize);

end

%% total field and phase unwrap
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
    switch phaseCombMethod
        % optimum weight method from Robinson et al. NMR Biomed 2017
        case 'optimum_weights'
            [totalField,fieldmapSD] = estimateTotalField(fieldMap,magn,...
                                    matrixSize,voxelSize,...
                                    'Unwrap',unwrap,'TE',TE,'B0',B0,'unit',unit,...
                                    'Subsampling',subsampling,'mask',mask);
                                
        % nonlinear ffitting method from MEDI toolbox
        case 'nonlinear_fit'
            [totalField,fieldmapSD] = estimateTotalField_MEDI(fieldMap,magn,...
                                    matrixSize,voxelSize,...
                                    'Unwrap',unwrap,'TE',TE,'B0',B0,'unit',unit,...
                                    'Subsampling',subsampling,'mask',mask);
    end
    
catch
    % if the above selected method is not working then do Laplacian with
    % optimum weights
    disp('The selected method is not supported in this system. Using Laplacian algorithm for phase unwrapping...')
    qsm_hub_AddMethodPath('laplacian');
    [totalField,fieldmapSD] = estimateTotalField(fieldMap,magn,matrixSize,voxelSize,...
                        'Unwrap','laplacian','TE',TE,'B0',B0,'unit',unit,...
                        'Subsampling',subsampling,'mask',mask);
end
                    
% exclude unreliable voxel, based on monoexponential decay model with
% single freuqnecy shift
r2s = arlo(TE,magn);
s0 = magn(:,:,:,1).*exp(r2s*TE(1));
shat = zeros([matrixSize length(TE)]);
for kt=1:length(TE)
    shat(:,:,:,kt) = s0.*exp(TE(kt)*(-r2s+1i*2*pi*totalField));
end
shat = bsxfun(@times,shat,conj(shat(:,:,:,1)));
s = bsxfun(@times,magn.*exp(1i*fieldMap),conj(magn(:,:,:,1).*exp(1i*fieldMap(:,:,:,1))));
relativeResidual=sum(abs(shat-s).^2,4)./sum(abs(s).^2,4);
relativeResidual(isnan(relativeResidual)) = exclude_threshold;
relativeResidual(isinf(relativeResidual)) = exclude_threshold;
maskReliable = relativeResidual < exclude_threshold;

% maskReliable = (fieldmapSD./norm(fieldmapSD(fieldmapSD~=0))) < exclude_threshold;
mask = and(mask,maskReliable);
                                        
disp('Saving unwrapped field map...');

% save the output
nii_totalField = make_nii_quick(outputNiftiTemplate,totalField);
nii_fieldmapSD = make_nii_quick(outputNiftiTemplate,fieldmapSD);
nii_newMask = make_nii_quick(outputNiftiTemplate,mask);
                    
save_untouch_nii(nii_totalField,[outputDir filesep prefix 'totalField.nii.gz']);
save_untouch_nii(nii_fieldmapSD,[outputDir filesep prefix 'fieldMapSD.nii.gz']);
save_untouch_nii(nii_newMask,[outputDir filesep prefix 'mask_new.nii.gz']);

if relativeResidual ~= 100
    nii_relativeResidual = make_nii_quick(outputNiftiTemplate,relativeResidual);
    save_untouch_nii(nii_relativeResidual,[outputDir filesep prefix 'relative_residual.nii.gz']);
end

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