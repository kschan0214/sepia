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

% get algorithm parameters
parameters      = check_and_set_algorithm_default(algorParam, headerAndExtraData);
parameters.voxelSize = voxelSize; % for MCPC-3D-S phase offset smoothing
method          = algorParam.unwrap.unwrapMethod;

% get magnitude
headerAndExtraData = check_and_set_SEPIA_header_data(headerAndExtraData);
magn = headerAndExtraData.magn;

% add path
sepia_addpath('ROMEO');

%% main
[totalField, fieldmapUnwrapAllEchoes] = ROMEO(wrappedField, magn, mask, parameters);
if algorParam.unwrap.isSaveUnwrappedEcho
    headerAndExtraData.fieldmapUnwrapAllEchoes = fieldmapUnwrapAllEchoes;
end
TEs = resphape(headerAndExtraData.te, 1,1,1,length(headerAndExtraData.te));
N_std = sqrt(sum(magn .* magn .* (TEs .* TEs), 4)); % TODO compare with SEPIA version
       
end

%% set default parameter if not specified
function parameters = check_and_set_algorithm_default(algorParam, headerAndExtraData)

parameters.TE = headerAndExtraData.te;
parameters.isSaveUnwrappedEcho = algorParam.unwrap.isSaveUnwrappedEcho;
parameters.calculateB0 = true;
parameters.useMag = true;

if contains(lower(algorParam.unwrap.offsetCorrect), 'bipolar')
    parameters.phaseOffsetCorrection = 'bipolar';
elseif contains(lower(algorParam.unwrap.offsetCorrect), 'on')
    parameters.phaseOffsetCorrection = 'on';
else
    parameters.phaseOffsetCorrection = 'off';
end

if contains(lower(algorParam.unwrap.mask), 'sepia')
    parameters.mask = 'file';
elseif contains(lower(algorParam.unwrap.mask), 'romeo')
    parameters.mask = 'robustmask';
else
    parameters.mask = 'nomask';
end

end