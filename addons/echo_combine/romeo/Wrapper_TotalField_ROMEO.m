%% [totalField, N_std, headerAndExtraData] = Wrapper_TotalField_ROMEO(wrappedField,mask,matrixSize,voxelSize,algorParam, headerAndExtraData)
%
% Input
% --------------
% wrappedField  : original single-/multi-echo wrapped phase image, in rad
% mask          : signal mask
% matrixSize    : size of the input image
% voxelSize     : spatial resolution of each dimension of the data, in mm
% algorParam    : structure contains fields with algorithm-specific parameter(s)
% headerAndExtraData : structure contains extra header info/data for the algorithm
%
% Output
% --------------
% totalField    : unwrapped total field, in Hz
% N_std         : noise standard deviation in the field map
% fieldmapUnwrapAllEchoes : unwrapped echo phase, in rad
%
% Description: This is a wrapper function to access ROMEO for SEPIA
%
% Korbinian Eckstein @ HFMR Vienna
% korbinian90@gmail.com
% Date created: 15 Juni 2021
% Date last modified:
%
%
function [totalField, N_std, headerAndExtraData] = Wrapper_TotalField_ROMEO(wrappedField,mask,matrixSize,voxelSize,algorParam, headerAndExtraData)
sepia_universal_variables;
headerAndExtraData = check_and_set_SEPIA_header_data(headerAndExtraData);

% add path
sepia_addpath('ROMEO');

% get algorithm parameters
parameters = check_and_set_algorithm_default(algorParam, headerAndExtraData, mask);
parameters.voxel_size = voxelSize; % for MCPC-3D-S phase offset smoothing
% Should create a suitable temporary directory on every machine
parameters.output_dir = fullfile(tempdir, 'romeo_tmp');
mkdir(parameters.output_dir);

%% main
[fieldmapUnwrapAllEchoes, totalField] = ROMEO(wrappedField, parameters);

% Remove all temp output files and the temp folder
rmdir(parameters.output_dir, 's')

%% Set additional outputs
if algorParam.unwrap.isSaveUnwrappedEcho
    headerAndExtraData.fieldmapUnwrapAllEchoes = fieldmapUnwrapAllEchoes;
end
TEs = reshape(headerAndExtraData.te, 1,1,1,length(headerAndExtraData.te));
mag = headerAndExtraData.magn;
N_std = sqrt(sum(mag .* mag .* (TEs .* TEs), 4)); % TODO compare with SEPIA version
       
end

%% set default parameters if not specified
function parameters = check_and_set_algorithm_default(algorParam, headerAndExtraData, mask)

parameters.TE = headerAndExtraData.te;
parameters.no_unwrapped_output = ~algorParam.unwrap.isSaveUnwrappedEcho;
parameters.calculate_B0 = true;
parameters.mag = headerAndExtraData.magn;

if contains(lower(algorParam.unwrap.offsetCorrect), 'bipolar')
    parameters.phase_offset_correction = 'bipolar';
elseif contains(lower(algorParam.unwrap.offsetCorrect), 'on')
    parameters.phase_offset_correction = 'on';
else
    parameters.phase_offset_correction = 'off';
end

if contains(lower(algorParam.unwrap.mask), 'sepia')
    parameters.mask = mask;
elseif contains(lower(algorParam.unwrap.mask), 'romeo')
    parameters.mask = 'robustmask';
else
    parameters.mask = 'nomask';
end

end