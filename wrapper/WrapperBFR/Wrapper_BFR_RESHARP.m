%% [RDF] = Wrapper_BFR_LBV(totalField,mask,matrixSize,voxelSize,algorParam, headerAndExtraData)
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
% Description: This is a wrapper function to access RESHARP for SEPIA
%
% Kwok-shing Chan @ DCCN
% k.chan@donders.ru.nl
% Date created: 8 March 2020 (v0.8.0)
% Date last modified:
%
%
function [RDF] = Wrapper_BFR_RESHARP(totalField,mask,matrixSize,voxelSize,algorParam, headerAndExtraData)
sepia_universal_variables;

% get algorithm parameters
algorParam  = check_and_set_algorithm_default(algorParam);
method      = algorParam.bfr.method;
radius     	= algorParam.bfr.radius;
alpha       = algorParam.bfr.alpha;

% add path
sepia_addpath('MEDI');
addpath(fullfile(SEPIA_HOME,'misc','background_removal','RESHARP'));

%% Display algorithm parameters
disp('The following parameters are being used...');
disp(['Radius(voxel) = ' num2str(radius)]);
disp(['Lambda        = ' num2str(alpha)]);

%% main
RDF         = RESHARP(totalField, mask, matrixSize, voxelSize, radius, alpha);
mask_RDF    = SMV(mask, matrixSize, voxelSize, radius)>0.999;
RDF         = RDF .* mask_RDF;
       
end

%% set default parameter if not specified
function algorParam2 = check_and_set_algorithm_default(algorParam)

algorParam2 = algorParam;

try algorParam2.bfr.radius	= algorParam.bfr.radius; 	catch; algorParam2.bfr.radius	= 4;  end
try algorParam2.bfr.alpha   = algorParam.bfr.alpha;     catch; algorParam2.bfr.alpha    = 0.1;   end

end