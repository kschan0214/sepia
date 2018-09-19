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
% Date last modified: 16 September 2018
%
%
function [totalField,fieldmapSD,mask]=UnwrapPhaseMacroIOWrapper(input,output,maskFullName,varargin)
%% add general Path
sepia_addpath;

%% define variables
prefix = 'squirrel_';
gyro = 42.57747892;
isInputDir = true;
% make sure the input only load once (first one)
isMagnLoad = false;
isPhaseLoad = false;

%% Check output directory exist or not
output_index = strfind(output, filesep);
outputDir = output(1:output_index(end));
if ~isempty(output(output_index(end)+1:end))
    prefix = [output(output_index(end)+1:end) '_'];
end

if exist(outputDir,'dir') ~= 7
    % if not then create the directory
    mkdir(outputDir);
end

%% Parse input argument using parse_varargin_QSMHub.m
[isInvert,~,isBET,isEddyCorrect,...
    phaseCombMethod,unwrap,subsampling,exclude_threshold,...
    ~,~,~,~,~,~,~,~,~,~,...
    ~,~,~,~,~,~,...
    ~,~,~,~,~,~,~,...
    ~,~,~,~,~,~,...
    ~,~] = parse_varargin_sepia(varargin);

%% Read input
disp('Reading data...');

% Step 1: check input for nifti files first
if isstruct(input)
    % Option 1: input are files
    inputNiftiList = input;
    isInputDir = false;
    
    % take the phase data directory as reference input directory 
    [inputDir,~,~] = fileparts(inputNiftiList(1).name);
else
    % Option 2: input is a directory
    inputDir = input;
    inputNiftiList = dir([inputDir '/*.nii*']);
end

% Step 2: load data
if ~isempty(inputNiftiList)
    
    if ~isInputDir
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %%%%%%%%%%%%%%%% Pathway 1: Input are NIfTI files %%%%%%%%%%%%%%%%%
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        
                        %%%%%%%%%% magnitude data %%%%%%%%%%
        if ~isempty(inputNiftiList(2).name)
            inputMagnNifti = load_untouch_nii([inputNiftiList(2).name]);
            % make sure the data is multi-echo magnitude data
            if size(inputMagnNifti.img,4) > 1
                magn = double(inputMagnNifti.img);
                isMagnLoad = true;
                disp('Magnitude data is loaded.');
            else
                error('QSM Hub only works with 4D data.');
            end
        else
            error('Please specify a 4D magnitude data.');
        end
                         %%%%%%%%%% phase data %%%%%%%%%%
        if ~isempty(inputNiftiList(1).name)
            inputPhaseNifti = load_untouch_nii([inputNiftiList(1).name]);
            fieldMap = double(inputPhaseNifti.img);
            % check whether phase data contains DICOM values or wrapped
            % phase value
            if max(fieldMap(:))>1000
                disp('Converting phase data from DICOM image value to radian unit...')
                fieldMap = DICOM2Phase(inputPhaseNifti);

                disp('Saving phase images in unit of radian...');
                save_nii_quick(inputPhaseNifti,fieldMap, [outputDir filesep prefix 'phase.nii.gz']);

            end
            isPhaseLoad = true;
            disp('Phase data is loaded.')
        else
            error('Please specify a 4D Phase data.');
        end
                        %%%%%%%%%% qsm hub header %%%%%%%%%%
        if ~isempty(inputNiftiList(4).name)
            load([inputNiftiList(4).name]);
            disp('Header data is loaded.');
        else
            error('Please specify a qsm_hub header.');
        end
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    else
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %%%%%%%%%% Pathway 2: Input is a directory with NIfTI %%%%%%%%%%
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        
        % loop all NIfTI files in the directory for magnitude and phase files
        for klist = 1:length(inputNiftiList)
                            %%%%%%%%%% magnitude data %%%%%%%%%%
            if ContainName(lower(inputNiftiList(klist).name),'magn') && ~isMagnLoad
                inputMagnNifti = load_untouch_nii([inputDir filesep inputNiftiList(klist).name]);
                % only load multi-echo magnitude data
                if size(inputMagnNifti.img,4) > 1
                    magn = double(inputMagnNifti.img);
                    isMagnLoad = true;
                    disp('Magnitude data is loaded.')
                end
            end

                            %%%%%%%%%% phase data %%%%%%%%%%
            if ContainName(lower(inputNiftiList(klist).name),'phase') && ~isPhaseLoad
                inputPhaseNifti = load_untouch_nii([inputDir filesep inputNiftiList(klist).name]);
                fieldMap = double(inputPhaseNifti.img);

                % if input fieldmap is directly converted from nifti converter
                % then converts the fieldmap to unit of radian and save to output dir
                if max(fieldMap(:))>1000
                    disp('Converting phase data from DICOM image value to radian unit ...')
                    fieldMap = DICOM2Phase(inputPhaseNifti);

                    disp('Saving phase images in unit of radian...');
                    save_nii_quick(inputPhaseNifti,fieldMap, [outputDir filesep prefix 'phase.nii.gz']);

                end
                isPhaseLoad = true;
                disp('Phase data is loaded.')
            end
        end

        % if no files matched the name format then displays error message
        if ~isMagnLoad
            error('No magnitude data is loaded. Please make sure the input directory contains NIfTI files with name *magn*');
        end
        if ~isPhaseLoad
            error('No phase data is loaded. Please make sure the input directory contains files with name *phase*');
        end

        %%%%%%%%%% qsm hub header file %%%%%%%%%%
        if ~isempty(dir([inputDir '/*header*']))
            % load header
            headerList = dir([inputDir '/*header*']);
            load([inputDir filesep headerList(1).name]);

            disp('Header data is loaded.');

        else
            disp('No header for qsm_hub is found. Creating synthetic header based on NIfTI header...');

            % create synthetic header in case no qsm_hub's header is found
            [B0,B0_dir,voxelSize,matrixSize,TE,delta_TE,CF]=SyntheticQSMHubHeader(inputMagnNifti);

            % look for text file for TEs information
            teTextFullName = dir([inputDir filesep '*txt']);
            if ~isempty(teTextFullName)
                te_ = readTEfromText([inputDir filesep teTextFullName(1).name]);
                te_ = te_(:);
                if ~isempty(te_)
                    TE = te_;
                    if length(TE) > 1
                        delta_TE = TE(2)-TE(1);
                    end
                end
            end

            % if no header file then save the synthetic header in output dir
            save([outputDir filesep 'SyntheticQSMhub_header'],'voxelSize','matrixSize','CF','delta_TE',...
            'TE','B0_dir','B0');

            disp('The synthetic header is saved in output directory.');

        end
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    end
    
    % store the header the NIfTI files, all following results will have
    % the same header
    outputNiftiTemplate = inputMagnNifti;
    
else
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%%%%%%%%%% Pathway 3: Input is a directory with DICOM %%%%%%%%%%%%%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    sepia_addpath('dicom');
    
    % if no nifti file then check for DICOM files
    [iField,voxelSize,matrixSize,CF,delta_TE,TE,B0_dir]=Read_DICOM(inputDir);
    
    % deprecated
    %     [iField,voxelSize,matrixSize,CF,delta_TE,TE,B0_dir]=Read_Siemens_DICOM_old(inputDir);
    
    B0 = CF/(gyro*1e6);
    
    % after testing with a couple of dataset it seems to me that the field
    % is inverted with DICOM input, so apply conjugate here
    fieldMap = angle(conj(iField));
    magn = abs(iField);
    
    isMagnLoad = true;
    isPhaseLoad = true;

    % save magnitude and phase images as nifti files
    disp('Saving DICOM data into NIfTI format...');
    
    % save magnitude and phase data as NIfTI_GZ format
    nii_fieldMap    = make_nii(single(fieldMap),    voxelSize,[],16);
    nii_magn        = make_nii(single(magn),        voxelSize,[],16);
    save_nii(nii_fieldMap,  [outputDir filesep prefix 'phase.nii.gz']);
    save_nii(nii_magn,      [outputDir filesep prefix 'magn.nii.gz']);
    % save important header in .mat format
    save([outputDir filesep prefix 'header.mat'],'voxelSize','matrixSize','CF','delta_TE',...
        'TE','B0_dir','B0');
    
    % reload the NIfTI template so that later can use save_untouch_nii for
    % all results
    outputNiftiTemplate = load_untouch_nii([outputDir filesep prefix 'magn.nii.gz']);
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
end

% in case user want to correct the frequency shift direction
if isInvert
    fieldMap = -fieldMap;
end

% display some header info
disp('Basic DICOM information');
disp(['Voxel size(x,y,z mm^3) =  ' num2str(voxelSize(1)) 'x' num2str(voxelSize(2)) 'x' num2str(voxelSize(3))]);
disp(['matrix size(x,y,z) =  ' num2str(matrixSize(1)) 'x' num2str(matrixSize(2)) 'x' num2str(matrixSize(3))]);
disp(['B0 direction(x,y,z) =  ' num2str(B0_dir(:)')]);
disp(['Field strength(T) =  ' num2str(B0)]);

% make sure the following variables are row vectors
matrixSize = matrixSize(:).';
voxelSize = voxelSize(:).';

%% get brain mask
mask = [];
maskList = dir([inputDir '/*mask*nii*']);
if ~isempty(maskFullName)
    % Option 1: mask file is provided
    mask = load_nii_img_only(maskFullName) > 0;
    
elseif ~isempty(maskList) 
    % Option 2: input directory contains NIfTI file with name '*mask*'
    inputMaskNii = load_untouch_nii([inputDir filesep maskList(1).name]);
	mask = inputMaskNii.img > 0;
    
end

% if no mask is found then display the following message
if isempty(mask) && ~isBET
    disp('No mask data is loaded, using FSL BET to obtain brain mask.');
end
    
% if BET is checked or no mask is found, run FSL's bet
if isempty(mask) || isBET
    sepia_addpath('bet');
    disp('Performing FSL BET...');
    
    % this is the BET functino provided with MEDI toolbox
    mask = BET(magn(:,:,:,1),matrixSize,voxelSize);
    
    disp('Saving brain mask...')
    save_nii_quick(outputNiftiTemplate,mask, [outputDir filesep prefix 'mask.nii.gz']);

end

%% total field and phase unwrap

% Step 0: Eddy current correction for bipolar readout
if isEddyCorrect
    disp('Correcting eddy current effect on bipolar readout data...');
    
    % BipolarEddyCorrect requries complex-valued input
    imgCplx = BipolarEddyCorrect(magn.*exp(1i*fieldMap),mask,unwrap);
    fieldMap = angle(imgCplx);
    
    % save the eddy current corrected output
    save_nii_quick(outputNiftiTemplate,fieldMap,    [outputDir filesep prefix 'phase_eddy-correct.nii.gz']);
    
end

% Step 1: Phase unwrapping and echo phase combination
disp('Calculating field map...');

% fix the output field map unit in Hz
unit = 'Hz';

% core of phase unwrapping
try 
    [totalField,fieldmapSD,fieldmapUnwrapAllEchoes] = estimateTotalField(fieldMap,magn,...
                        matrixSize,voxelSize,...
                        'method',phaseCombMethod,'Unwrap',unwrap,...
                        'TE',TE,'B0',B0,'unit',unit,...
                        'Subsampling',subsampling,'mask',mask);
                    
    if ~isempty(fieldmapUnwrapAllEchoes)
        save_nii_quick(outputNiftiTemplate,fieldmapUnwrapAllEchoes,[outputDir filesep prefix 'unwrapped-phase.nii.gz']);
    end
    
catch
    % if the above selected method is not working then do Laplacian with
    % optimum weights
    disp('The selected method is not supported in this system. Using Laplacian algorithm for phase unwrapping...')
    
    unwrap = 'laplacian'; exclude_threshold = Inf;
    sepia_addpath(unwrap);
    
    [totalField,fieldmapSD] = estimateTotalField(fieldMap,magn,matrixSize,voxelSize,...
                            'Unwrap',unwrap,'TE',TE,'B0',B0,'unit',unit,...
                            'Subsampling',subsampling,'mask',mask);
end
                    
% Step 2: exclude unreliable voxel, based on monoexponential decay model with
% single freuqnecy shift
r2s = R2star_trapezoidal(magn,TE);
relativeResidual = ComputeResidualGivenR2sFieldmap(TE,r2s,totalField,magn.*exp(1i*fieldMap));
maskReliable = relativeResidual < exclude_threshold;

% threshold fieldmapSD with the reliable voxel mask
fieldmapSD = fieldmapSD .* maskReliable;

% 20180815: test with creating weights using relativeResidual
% weightResidual = 1-(relativeResidual./exclude_threshold);
% weightResidual(weightResidual>1) = 1;
% weightResidual(weightResidual<0) = 0;

% deprecated
% mask = and(mask,maskReliable);

% computing weights
wmap = 1./fieldmapSD;
wmap(isinf(wmap)) = 0;
wmap(isnan(wmap)) = 0;
wmap = wmap./max(wmap(and(mask>0,maskReliable>0)));
wmap = wmap .* and(mask>0,maskReliable);
             
% save the output                           
disp('Saving unwrapped field map...');

save_nii_quick(outputNiftiTemplate,totalField,  [outputDir filesep prefix 'total-field.nii.gz']);
save_nii_quick(outputNiftiTemplate,fieldmapSD,  [outputDir filesep prefix 'noise-sd.nii.gz']);
save_nii_quick(outputNiftiTemplate,wmap,  [outputDir filesep prefix 'weights.nii.gz']);

if relativeResidual ~= Inf
    save_nii_quick(outputNiftiTemplate,maskReliable,   	[outputDir filesep prefix 'mask-reliable.nii.gz']);
    save_nii_quick(outputNiftiTemplate,relativeResidual,[outputDir filesep prefix 'relative-residual.nii.gz']);
end

disp('Done!');

end
