%% [unwrappedField] = Wrapper_Unwrap_ROMEO(wrappedField,mask,matrixSize,voxelSize,algorParam, headerAndExtraData)
%
% Input
% --------------
% wrappedField  : wrapped field map
% mask          : signal mask
% matrixSize    : size of the input image
% voxelSize     : spatial resolution of each dimension of the data, in mm
% algorParam    : structure contains fields with algorithm-specific parameter(s)
% headerAndExtraData : structure contains extra header info/data for the algorithm
%
% Output
% --------------
% unwrappedField: local field map
%
% Description: This is a wrapper function to access ROMEO for SEPIA
%
% Korbinian Eckstein @ HFMR Vienna
% korbinian90@gmail.com
% Date created: 1 Juni 2021
% Date last modified:
%
%
function [unwrappedField] = Wrapper_Unwrap_ROMEO(wrappedField,mask,matrixSize,voxelSize,algorParam, headerAndExtraData)
sepia_universal_variables;

% get algorithm parameters
algorParam      = check_and_set_algorithm_default(algorParam);
method          = algorParam.unwrap.unwrapMethod;

% get magnitude
headerAndExtraData = check_and_set_SEPIA_header_data(headerAndExtraData);
magn = headerAndExtraData.magn(:,:,:,1);

% add path
sepia_addpath('ROMEO');

%% main
unwrappedField  = ROMEO(wrappedField, magn, mask);
       
end

%% set default parameter if not specified
function algorParam2 = check_and_set_algorithm_default(algorParam)

algorParam2 = algorParam;

try algorParam2.unwrap.subsampling = algorParam.unwrap.subsampling; catch; algorParam2.unwrap.subsampling  = 1;  end

end