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
function chi = qsmMacroIOWrapper(inputDir,outputDir,varargin)
%% add general Path
qsm_hub_AddMethodPath % qsm_hub_AddPath;

gyro = 42.57747892;

%% Check output directory exist or not
if exist(outputDir,'dir') ~= 7
    % if not then create the directory
    mkdir(outputDir);
end

%% Parse input argument
[~,maskFullName,~,~,~,~,~,~,~,~,~,~,~,~,QSM_method,QSM_threshold,QSM_lambda,...
QSM_optimise,QSM_tol,QSM_maxiter,QSM_tol1,QSM_tol2,QSM_padsize,QSM_mu1,QSM_solver,QSM_constraint,...
~,QSM_radius,QSM_zeropad,QSM_wData,QSM_wGradient,QSM_lambdaCSF,QSM_isSMV,QSM_merit,~] = parse_varargin_QSMHub(varargin);

%% Read input
disp('Reading data...');
% look for nifti files first
inputNiftiList = dir([inputDir '/*.nii*']);
if ~isempty(inputNiftiList)
    % look for magnitude,localField and fieldmapSD files
    fieldmapSD = []; magn = [];
    for klist = 1:length(inputNiftiList)
        if contains(lower(inputNiftiList(klist).name),'localfield')
            inputLocalFieldNifti = load_untouch_nii([inputDir filesep inputNiftiList(klist).name]);
            localField = double(inputLocalFieldNifti.img);
        end
        if contains(lower(inputNiftiList(klist).name),'magn') && ~contains(lower(inputNiftiList(klist).name),'brain')
            inputMagnNifti = load_untouch_nii([inputDir filesep inputNiftiList(klist).name]);
            magn = double(inputMagnNifti.img);
        end
        if contains(lower(inputNiftiList(klist).name),'fieldmapsd')
            inputFieldMapSDNifti = load_untouch_nii([inputDir filesep inputNiftiList(klist).name]);
            fieldmapSD = double(inputFieldMapSDNifti.img);
        end
    end
    % create synthetic header in case no header is provided
    [B0,~,~,~,TE,delta_TE,CF]=SyntheticQSMHubHeader(magn);
    a=qGetR([0, inputLocalFieldNifti.hdr.hist.quatern_b,inputLocalFieldNifti.hdr.hist.quatern_c,inputLocalFieldNifti.hdr.hist.quatern_d]);
    B0_dir = -a(3,:);
    voxelSize = inputLocalFieldNifti.hdr.dime.pixdim(2:4);
    matrixSize = inputLocalFieldNifti.hdr.dime.dim(2:4);
    % look for header file
    if ~isempty(dir([inputDir '/*header*']))
        headerList = dir([inputDir '/*header*']);
        load([inputDir filesep headerList(1).name]);
    else
        disp('No QSMHub header file detected. Using synthetic header parameters...');
    end
    if isempty(magn)
        disp('No magnitude data is found.');
        magn = ones(matrixSize);
    end
    if isempty(fieldmapSD)
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
    error('This standalone only reads NIfTI format input data (nii or nii.gz).');
end

% display some header info
disp('Basic DICOM information');
disp(['Voxel size(x,y,z mm^3) =  ' num2str(voxelSize(1)) 'x' num2str(voxelSize(2)) 'x' num2str(voxelSize(3))]);
disp(['matrix size(x,y,z) =  ' num2str(matrixSize(1)) 'x' num2str(matrixSize(2)) 'x' num2str(matrixSize(3))]);
disp(['B0 direction(x,y,z) =  ' num2str(B0_dir(:)')]);
disp(['Field strength(T) =  ' num2str(B0)]);

%% get brain mask
maskList = dir([inputDir '/*mask_final*']);
if isempty(maskList)
    maskList = dir([inputDir '/*mask*']);
end
% first read mask if file is provided
if ~isempty(maskFullName)
    maskFinal = load_nii_img_only(maskFullName) > 0;
elseif ~isempty(maskList) 
    % read mask if input directory contains 'mask'
    inputMaskNii = load_untouch_nii([inputDir filesep maskList(1).name]);
	maskFinal = inputMaskNii.img > 0;
else
    error('No mask is found. Pleasee specific your mask file or put it inside the input directory.');
end

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