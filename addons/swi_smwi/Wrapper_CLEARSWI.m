%% Wrapper_CLEARSWI(headerAndExtraData,output_structure,algorParam)
%
% Input
% --------------
% headerAndExtraData : structure contains extra header info/data for the algorithm
% output_structure   : structure contains output file path and output NIfTI header
% algorParam         : structure contains fields with algorithm-specific parameter(s)
%
% Output
% --------------
%
% Description: This is a wrapper function to access CLEARSWI for SEPIA
%
% Korbinian Eckstein @ UQ
% korbinian90@gmail.com
% Date created: 26 September 2022
% Date modified: 
%
%
function Wrapper_CLEARSWI(headerAndExtraData,output_structure,algorParam)
sepia_universal_variables;

% add path
sepia_addpath('MRITOOLS');

% get algorithm parameters
algorParam = check_and_set_algorithm_default(algorParam);
parameters.phase_scaling_type = algorParam.swismwi.phaseScalingType;
parameters.phase_scaling_strength = algorParam.swismwi.phaseScalingStrength;
parameters.filter_size = algorParam.swismwi.filterSize;
parameters.unwrapping_algorithm = algorParam.swismwi.unwrappingAlgorithm;

parameters.mag_combine = algorParam.swismwi.echoCombineMethod;
if ~isempty(algorParam.swismwi.echoCombineMethodAdd)
    parameters.mag_combine = [parameters.mag_combine " " algorParam.swismwi.echoCombineMethodAdd];
end
parameters.echoes = algorParam.swismwi.echoes;
parameters.softplusScaling = 'on';
if ~algorParam.swismwi.softplusScaling
    parameters.softplusScaling = 'off';
end

parameters.sensitivityCorrection = 'on';
if ~algorParam.swismwi.sensitivityCorrection
    parameters.sensitivityCorrection = 'off';
end

ismIP       = algorParam.swismwi.ismIP;
slice_mIP   = algorParam.swismwi.slice_mIP;

% get extra data such as magnitude/weights/B0 direction/TE/etc.
headerAndExtraData                      = check_and_set_SEPIA_header_data(headerAndExtraData);
% check if the phase image is in wrapped phase range
headerAndExtraData.availableFileList	= io_true_phase_value(headerAndExtraData.availableFileList, output_structure.outputFileFullPrefix);
% b0dir = headerAndExtraData.sepia_header.B0_dir;
% matrixSize    = headerAndExtraData.sepia_header.matrixSize;
phase = get_variable_from_headerAndExtraData(headerAndExtraData, 'phasechi');
magn  = get_variable_from_headerAndExtraData(headerAndExtraData, 'magnitude');

outputNiftiTemplate     = output_structure.outputNiftiTemplate;
outputFileFullPrefix    = output_structure.outputFileFullPrefix;
% add path
% sepia_addpath;

%% SWI
disp('Computing CLEARSWI...');

[swi, ~] = CLEARSWI(magn,phase,parameters);

% minimum intensity projection
if ismIP
    mIP = zeros(size(swi_phase));
    for kz = 1:size(swi_phase,3) - slice_mIP
        mIP(:,:,kz,:) = min(swi(:,:,kz:kz+slice_mIP-1,:),[],3);
    end
end

%% save nii
disp('Saving SWI results...');
save_nii_quick(outputNiftiTemplate, swi, [outputFileFullPrefix 'clearswi.nii.gz']);
if ismIP
    save_nii_quick(outputNiftiTemplate, mIP, [outputFileFullPrefix 'clearswi-minIP.nii.gz']);
end

disp('Done!');

end

%% set default parameter if not specified
function algorParam2 = check_and_set_algorithm_default(algorParam)

algorParam2 = algorParam;

try algorParam2.swismwi.phaseScalingType = algorParam.swismwi.phaseScalingType; catch; algorParam.swismwi.phaseScalingType = 'tanh'; end
try algorParam2.swismwi.phaseScalingStrength = algorParam.swismwi.phaseScalingStrength; catch; algorParam.swismwi.phaseScalingStrength = 4; end    
try algorParam2.swismwi.filterSize	= algorParam.swismwi.filterSize;	catch; algorParam2.swismwi.filterSize	= '[4,4,0]';	end
try algorParam2.swismwi.unwrappingAlgorithm = algorParam.swismwi.unwrappingAlgorithm; catch; algorParam.swismwi.unwrappingAlgorithm = 'laplacian'; end

try algorParam2.swismwi.echoCombineMethod = algorParam.swismwi.echoCombineMethod; catch; algorParam.swismwi.echoCombineMethod = 'SNR'; end
try algorParam2.swismwi.echoCombineMethodAdd = algorParam.swismwi.echoCombineMethodAdd; catch; algorParam.swismwi.echoCombineMethodAdd = ''; end
try algorParam2.swismwi.echoes = algorParam.swismwi.echoes; catch; algorParam.swismwi.echoes = 'all'; end
try algorParam2.swismwi.softplusScaling = algorParam.swismwi.softplusScaling; catch; algorParam.swismwi.softplusScaling = true; end
try algorParam2.swismwi.sensitivityCorrection = algorParam.swismwi.sensitivityCorrection; catch; algorParam.swismwi.sensitivityCorrection = true; end

try algorParam2.swismwi.ismIP   	= algorParam.swismwi.ismIP;         catch; algorParam2.swismwi.ismIP        = true;	end
try algorParam2.swismwi.slice_mIP  	= algorParam.swismwi.slice_mIP; 	catch; algorParam2.swismwi.slice_mIP	= 4;	end

end

%% I/O: convert phase to radian unit if required
function availableFileList          = io_true_phase_value(availableFileList, outputFileFullPrefix)

% load phase image to check if the 
phaseNIFTI = load_untouch_nii(availableFileList.phasechi);
phaseIMG   = load_nii_img_only(availableFileList.phasechi);

if abs(max(phaseIMG(:))-pi)>0.1 || abs(min(phaseIMG(:))-(-pi))>0.1 % allow small differences possibly due to data stype conversion or DICOM digitisation
% if abs(max(phaseNIFTI.img(:))-pi)>0.1 || abs(min(phaseNIFTI.img(:))-(-pi))>0.1 % allow small differences possibly due to data stype conversion or DICOM digitisation
    
    outputPhaseRadian_filename = [outputFileFullPrefix 'part-phase_rad.nii.gz'];

    disp('Values of input phase map exceed the range of [-pi,pi]. DICOM value is assumed.')
    fprintf('Rescaling phase data from DICOM image value to wrapped radian unit...')
    phase = DICOM2Phase(phaseNIFTI);
    fprintf('Done.\n')

    fprintf('Saving phase images in unit of radian...');
    save_nii_quick(phaseNIFTI, phase, outputPhaseRadian_filename);
    fprintf('Done.\n')
    
    % update the phase data for subsequent processing
    availableFileList.phasechi = outputPhaseRadian_filename;
    
end

end