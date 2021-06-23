%% [unwrappedField] = Wrapper_TotalField_ROMEO(fieldMap,mask,matrixSize,voxelSize,algorParam,headerAndExtraData)
%
% Input
% --------------
% fieldMap      : original single-/multi-echo wrapped phase image, in rad
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
method          = algorParam.unwrap.unwrapMethod;

% get magnitude
headerAndExtraData = check_and_set_SEPIA_header_data(headerAndExtraData);
magn = headerAndExtraData.magn;

% add path
sepia_addpath('ROMEO');

%% main
% TODO mcpc3ds, bipolar correction
totalField  = ROMEO(wrappedField, magn, mask, parameters);
N_std = totalField;
       
end

%% set default parameter if not specified
function parameters = check_and_set_algorithm_default(algorParam, headerAndExtraData)

parameters.TE = headerAndExtraData.te;

try parameters.unwrap.subsampling = algorParam.unwrap.subsampling; catch; parameters.unwrap.subsampling  = 1;  end

end