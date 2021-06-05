%% [chi] = Wrapper_QSM_STISuiteiLSQRv3(localField,mask,matrixSize,voxelSize,algorParam, headerAndExtraData)
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
% Description: This is a wrapper function to access iLSQR from STI suite for SEPIA
%
% Kwok-shing Chan @ DCCN
% k.chan@donders.ru.nl
% Date created: 8 March 2020
% Date last modified:
%
%
function [chi] = Wrapper_QSM_STISuiteiLSQRv3(localField,mask,matrixSize,voxelSize,algorParam, headerAndExtraData)
sepia_universal_variables;

% get algorithm parameters
algorParam = check_and_set_algorithm_default(algorParam);
method     = algorParam.qsm.method;
params.Kthreshold   = algorParam.qsm.threshold;
params.niter        = algorParam.qsm.maxiter;
params.tol_step1    = algorParam.qsm.tol1;
params.tol_step2    = algorParam.qsm.tol2;
params.padsize      = algorParam.qsm.padsize;

% get extra data such as magnitude/weights/B0 direction/TE/etc.
headerAndExtraData = check_and_set_SEPIA_header_data(headerAndExtraData);
params.H            = headerAndExtraData.b0dir(:).';
params.TE           = headerAndExtraData.delta_TE*1e3;
params.B0           = headerAndExtraData.b0;
params.voxelsize    = double(voxelSize(:).');

% add path
sepia_addpath('STISuite');

%% Display algorithm parameters
disp('The following parameters are being used...');
disp(['Threshold      	= ' num2str(params.Kthreshold)]);
disp(['Max. iteration 	= ' num2str(params.niter)]);
disp(['Tolerance 1    	= ' num2str(params.tol_step1)]);
disp(['Tolerance 2     	= ' num2str(params.tol_step2)]);
disp(['Pas size         = [' num2str(params.padsize(1)) ', ' num2str(params.padsize(2)) ', ' num2str(params.padsize(3)) ']']);

%% main
% The order of local field values doesn't affect the result of chi  
% in STI suite v3 implementation, i.e. 
% chi = method(localField,...) = method(localField*C,...)/C, where
% C is a constant.
% Therefore, because of the scaling factor in their implementation,
% the local field map is converted to rad
localField = localField * 2*pi * (params.TE*1e-3); 

% double precision is requried for this function
chi = QSM_iLSQR(double(localField),double(mask),'params',params);
        
end

%% set default parameter if not specified
function algorParam2 = check_and_set_algorithm_default(algorParam)

algorParam2 = algorParam;

try algorParam2.qsm.threshold   = algorParam.qsm.threshold; catch; algorParam2.qsm.threshold    = 0.25; end
try algorParam2.qsm.tol1        = algorParam.qsm.tol1;      catch; algorParam2.qsm.tol1         = 0.01; end
try algorParam2.qsm.tol2        = algorParam.qsm.tol2;      catch; algorParam2.qsm.tol2         = 0.001; end
try algorParam2.qsm.maxiter     = algorParam.qsm.maxiter;   catch; algorParam2.qsm.maxiter      = 100; end
try algorParam2.qsm.padsize     = algorParam.qsm.padsize;   catch; algorParam2.qsm.padsize      = [4,4,4]; end

end