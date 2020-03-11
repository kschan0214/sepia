%% [chi] = Wrapper_QSM_CFSL2(localField,mask,matrixSize,voxelSize,algorParam,headerAndExtraData)
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
% Date last modified:
%
%
function [chi] = Wrapper_QSM_CFSL2(localField,mask,matrixSize,voxelSize,algorParam,headerAndExtraData)
sepia_universal_variables;

% get algorithm parameters
algorParam  = check_and_set_algorithm_default(algorParam);
method     = algorParam.qsm.method;
lambda      = algorParam.qsm.lambda;
optimise    = algorParam.qsm.optimise;

% get extra data such as magnitude/weights/B0 direction/TE/etc.
headerAndExtraData = check_and_set_SEPIA_header_data(headerAndExtraData);
b0dir = headerAndExtraData.b0dir;
b0    = headerAndExtraData.b0;

% add path
sepia_addpath(method);

%% Display algorithm parameters
disp('The following parameters are being used...');
disp(['Regularisation lambda = ' num2str(lambda)]);
disp(['Self-optimisation     = ' num2str(optimise)]);

%% main
[chi, lamdaOptimal] = qsmClosedFormL2(localField,mask,matrixSize,voxelSize,...
            'lambda',lambda,'optimise',optimise,'b0dir',b0dir);
        
        
% convert from Hz to ppm
chi = chi/(b0*gyro);

end

%% set default parameter if not specified
function algorParam2 = check_and_set_algorithm_default(algorParam)

algorParam2 = algorParam;

try algorParam2.qsm.lambda      = algorParam.qsm.lambda;    catch; algorParam2.qsm.lambda   = 0.1; end
try algorParam2.qsm.optimise    = algorParam.qsm.optimise;  catch; algorParam2.qsm.optimise = false; end

end