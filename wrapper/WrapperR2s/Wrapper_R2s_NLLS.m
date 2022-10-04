%% [r2s,t2s,m0] = Wrapper_R2s_Trapezoidal(magn,mask,matrixSize,voxelSize,algorParam,headerAndExtraData)
%
% Input
% --------------
% localField    : local field map (tissue fields), in Hz
% mask          : signal mask
% matrixSize    : size of the input image
% voxelSize     : spatial resolution of each dimension of the data, in mm
% algorParam    : structure contains fields with algorithm-specific parameter(s)
% headerAndExtraData : structure contains extra header info/data for the algorithm
%
% Output
% --------------
% chi           : magnetic susceptibility map, in ppm
%
% Description: This is a wrapper function to access closed-form solution for SEPIA
%
% Kwok-shing Chan @ DCCN
% k.chan@donders.ru.nl
% Date created: 8 March 2020
% Date modified: 13 August 2021 (v1.0)
%
%
function [r2s,t2s,m0] = Wrapper_R2s_NLLS(magn,te,mask,algorParam,headerAndExtraData)
sepia_universal_variables;

% get algorithm parameters
algorParam  = check_and_set_algorithm_default(algorParam);
method      = algorParam.r2s.method;
isParallel  = algorParam.r2s.isParallel;
% NUM_MAGN    = algorParam.r2s.NUM_MAGN;

% get extra data such as magnitude/weights/B0 direction/TE/etc.
% headerAndExtraData = check_and_set_SEPIA_header_data(headerAndExtraData);
% te = headerAndExtraData.headerAndExtraData.TE;

% add path
sepia_addpath;

%% Display algorithm parameters
% disp('The following parameters are being used...');
% disp(['Regularisation lambda = ' num2str(lambda)]);

%% main
% if isempty(NUM_MAGN)
NUM_MAGN = length(te);
% end

[r2s,t2s,m0] = R2star_NLLS(magn,te,mask,isParallel,NUM_MAGN);
        
r2s = r2s.*mask;
t2s = t2s.*mask;
m0  = m0.*mask;

end

%% set default parameter if not specified
function algorParam2 = check_and_set_algorithm_default(algorParam)

algorParam2 = algorParam;

try algorParam2.r2s.isParallel      = algorParam.r2s.isParallel;    catch; algorParam2.r2s.isParallel = false; end
% try algorParam2.r2s.NUM_MAGN        = algorParam.r2s.NUM_MAGN;      catch; algorParam2.r2s.NUM_MAGN   = []; end

end