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
% Date modified: 26 August 2018
% Date modified: 29 March 2019
%
%
function [localField,maskFinal] = BackgroundRemovalMacroIOWrapper(input,output,maskFullName,algorParam)
%% add general Path
sepia_addpath;

%% define variables
prefix = 'sepia_';
gyro = 42.57747892;
isInputDir = true;
% make sure the input only load once (first one)
isTotalFieldLoad = false;
isFieldmapSDLoad = false;

%% Check if output directory exists 
output_index = strfind(output, filesep);
outputDir = output(1:output_index(end));
% get prefix
if ~isempty(output(output_index(end)+1:end))
    prefix = [output(output_index(end)+1:end) '_'];
end
% if the output directory does not exist then create the directory
if exist(outputDir,'dir') ~= 7
    mkdir(outputDir);
end

%% Check and set default algorithm parameters
algorParam = CheckAndSetDefault(algorParam);
isGPU           = algorParam.general.isGPU;
BFR             = algorParam.bfr.method;
refine          = algorParam.bfr.refine;
BFR_tol         = algorParam.bfr.tol;
BFR_depth       = algorParam.bfr.depth;
BFR_peel        = algorParam.bfr.peel;
BFR_iteration	= algorParam.bfr.iteration;
BFR_padSize     = algorParam.bfr.padSize;
BFR_radius      = algorParam.bfr.radius;
BFR_alpha       = algorParam.bfr.alpha;
BFR_threshold   = algorParam.bfr.threshold;

%% Read input
disp('Reading data...');

% Step 1: check input for nifti files first
if isstruct(input)
    % Option 1: input are files
    inputNiftiList = input;
    isInputDir = false;
    
    % take the total field map directory as reference input directory 
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
        
                        %%%%%%%%%% Total field map %%%%%%%%%%
        if ~isempty(inputNiftiList(1).name)
            inputTotalFieldNifti = load_untouch_nii([inputNiftiList(1).name]);
            totalField = double(inputTotalFieldNifti.img);
            isTotalFieldLoad = true;
            
            disp('Total field map is loaded.')
        else
            error('Please specify a 3D total field map.');
        end
        
                         %%%%%%%%%% Fieldmapsd data %%%%%%%%%%
        if ~isempty(inputNiftiList(3).name)
            inputFieldMapSDNifti = load_untouch_nii([inputNiftiList(3).name]);
            fieldmapSD = double(inputFieldMapSDNifti.img);
            isFieldmapSDLoad = true;
            
            disp('Noise SD data is loaded.')
        else
            disp('No field map standard deviation data is loaded.');
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
        
        % loop all NIfTI files in the directory for total field map and fieldmap SD (optional)
        for klist = 1:length(inputNiftiList)

                        %%%%%%%%%% Total field map %%%%%%%%%%
            if ContainName(inputNiftiList(klist).name,'total-field') && ~isTotalFieldLoad
                inputTotalFieldNifti = load_untouch_nii([inputDir filesep inputNiftiList(klist).name]);
                totalField = double(inputTotalFieldNifti.img);
                isTotalFieldLoad = true;

                disp('Total field map is loaded.')
            end

                         %%%%%%%%%% Fieldmapsd data %%%%%%%%%%
            if ContainName(inputNiftiList(klist).name,'noise-sd') && ~isFieldmapSDLoad
                inputFieldMapSDNifti = load_untouch_nii([inputDir filesep inputNiftiList(klist).name]);
                fieldmapSD = double(inputFieldMapSDNifti.img);
                isFieldmapSDLoad = true;

                disp('Noise SD data is loaded.')
            end
        end

        % if no files matched the name format then displays error message
        if ~isTotalFieldLoad
            error('No total field map is loaded. Please make sure the input directory contains files with name *totalfield*');
        end
        if ~isFieldmapSDLoad
            disp('No field map standard deviation data is loaded.');
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
            [B0,B0_dir,voxelSize,matrixSize,TE,delta_TE,CF]=SyntheticQSMHubHeader(inputTotalFieldNifti);

            % if no header file then save the synthetic header in output dir
            save([outputDir filesep 'SyntheticQSMhub_header'],'voxelSize','matrixSize','CF','delta_TE',...
            'TE','B0_dir','B0');

            disp('The synthetic header is saved in output directory.');

        end

        % if no fieldmapSD found then creates one with all voxels have the same value
        if ~isFieldmapSDLoad
            fieldmapSD = ones(matrixSize) * 0.01;
        end
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    end
    
    % store the header the NIfTI files, all following results will have
    % the same header
    outputNiftiTemplate = inputTotalFieldNifti;
    
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
    % Option 1: mask file is provided
    mask = load_nii_img_only(maskFullName) > 0;
    
elseif ~isempty(maskList) 
    % Option 2: input directory contains NIfTI file with name '*mask*'
    inputMaskNii = load_untouch_nii([inputDir filesep maskList(1).name]);
	mask = inputMaskNii.img > 0;
    
else
    % display error message if nothing is found
    error('No mask file is found. Please specify your mask file or put it in the input directory.');
    
end

%% Background field removal
disp('Recovering local field...');

% core of background field removal
if isGPU
    localField = cuBackgroundRemovalMacro(...
                        totalField,mask,matrixSize,voxelSize,...
                        'method',BFR,'refine',refine,'tol',BFR_tol,'depth',BFR_depth,...
                        'peel',BFR_peel,'b0dir',B0_dir,'iteration',BFR_iteration,'padsize',...
                        BFR_padSize,'noisestd',fieldmapSD,'radius',BFR_radius,'alpha',BFR_alpha,...
                        'threshold',BFR_threshold); 
else
    localField = BackgroundRemovalMacro(...
                        totalField,mask,matrixSize,voxelSize,...
                        'method',BFR,'refine',refine,'tol',BFR_tol,'depth',BFR_depth,...
                        'peel',BFR_peel,'b0dir',B0_dir,'iteration',BFR_iteration,'padsize',...
                        BFR_padSize,'noisestd',fieldmapSD,'radius',BFR_radius,'alpha',BFR_alpha,...
                        'threshold',BFR_threshold); 
end
  
% generate new mask based on backgroudn field removal result
maskFinal = localField ~=0;
  
% save results
disp('Saving local field map...');

save_nii_quick(outputNiftiTemplate,localField, [outputDir filesep prefix 'local-field.nii.gz']);
save_nii_quick(outputNiftiTemplate,maskFinal,  [outputDir filesep prefix 'mask-qsm.nii.gz']);

disp('Done!');

end

function algorParam2 = CheckAndSetDefault(algorParam)
algorParam2 = algorParam;
try algorParam2.general.isGPU  	= algorParam.general.isGPU;	catch; algorParam2.general.isGPU = false;   end
% default background field removal method is VSHARP
try algorParam2.bfr.method      = algorParam.bfr.method;   	catch; algorParam2.bfr.method = 'vsharpsti';end
try algorParam2.bfr.radius      = algorParam.bfr.radius; 	catch; algorParam2.bfr.radius = 10;         end
try algorParam2.bfr.refine      = algorParam.bfr.refine;   	catch; algorParam2.bfr.refine = false;      end
% for the rest, if the parameter does not exist then initiates it with an empty array
try algorParam2.bfr.tol         = algorParam.bfr.tol;     	catch; algorParam2.bfr.tol = [];            end
try algorParam2.bfr.depth       = algorParam.bfr.depth;  	catch; algorParam2.bfr.depth = [];          end
try algorParam2.bfr.peel        = algorParam.bfr.peel;   	catch; algorParam2.bfr.peel = [];           end
try algorParam2.bfr.iteration   = algorParam.bfr.iteration;	catch; algorParam2.bfr.iteration = [];      end
try algorParam2.bfr.padSize     = algorParam.bfr.padSize; 	catch; algorParam2.bfr.padSize = [];        end
try algorParam2.bfr.alpha       = algorParam.bfr.alpha;    	catch; algorParam2.bfr.alpha = [];          end
try algorParam2.bfr.threshold   = algorParam.bfr.threshold;	catch; algorParam2.bfr.threshold = [];      end

end