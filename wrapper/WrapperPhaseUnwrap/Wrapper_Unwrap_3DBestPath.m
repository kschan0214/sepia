%% [RDF] = Wrapper_BFR_PDF(totalField,mask,matrixSize,voxelSize,algorParam, headerAndExtraData)
%
% Input
% --------------
% totalField    : total field map (background + tissue fields), in Hz
% mask          : signal mask
% matrixSize    : size of the input image
% voxelSize     : spatial resolution of each dimension of the data, in mm
% algorParam    : structure contains fields with algorithm-specific parameter(s)
% headerAndExtraData : structure contains extra header info/data for the algorithm
%
% Output
% --------------
% RDF           : local field map
%
% Description: This is a wrapper function to access PDF for SEPIA
%
% Kwok-shing Chan @ DCCN
% k.chan@donders.ru.nl
% Date created: 8 March 2020 (v0.8.0)
% Date last modified:
%
%
function [unwrappedField] = Wrapper_Unwrap_3DBestPath(wrappedField,mask,matrixSize,voxelSize,algorParam, headerAndExtraData)
sepia_universal_variables;

% get algorithm parameters
algorParam      = check_and_set_algorithm_default(algorParam);
method          = algorParam.unwrap.unwrapMethod;

% add path
sepia_addpath('MEDI');
addpath(fullfile(SEPIA_HOME,'misc','phase_unwrap','unwrapBestpath3D'));

%% main
unwrappedField = UnwrapPhase_3DBestPath(wrappedField,mask,matrixSize);
       
end

%% set default parameter if not specified
function algorParam2 = check_and_set_algorithm_default(algorParam)

algorParam2 = algorParam;

% try algorParam2.bfr.tol         = algorParam.bfr.tol;       catch; algorParam2.bfr.tol      = 0.1;  end

end