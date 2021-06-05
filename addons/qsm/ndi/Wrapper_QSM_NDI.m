%% [chi] = Wrapper_QSM_NDI(localField,mask,matrixSize,voxelSize,algorParam, headerAndExtraData)
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
% Description: This is a wrapper function to access NDI for SEPIA
%
% Kwok-shing Chan @ DCCN
% k.chan@donders.ru.nl
% Date created: 8 March 2020
% Date last modified:
%
%
function [chi] = Wrapper_QSM_NDI(localField,mask,matrixSize,voxelSize,algorParam, headerAndExtraData)
sepia_universal_variables;

% get algorithm parameters
algorParam  = check_and_set_algorithm_default(algorParam);
method      = algorParam.qsm.method;
tol         = algorParam.qsm.tol;
stepSize    = algorParam.qsm.stepSize;
maxiter     = algorParam.qsm.maxiter;

% get extra data such as magnitude/weights/B0 direction/TE/etc.
headerAndExtraData = check_and_set_SEPIA_header_data(headerAndExtraData);
b0dir = headerAndExtraData.b0dir;
b0    = headerAndExtraData.b0;
wmap  = headerAndExtraData.weights;
magn  = headerAndExtraData.magn;

% add path
sepia_addpath;

%% preparation
% if both data are loaded
if ~isempty(magn) && ~isempty(wmap)
    disp('Both weighting map and magnitude images are loaded.');
    disp('Only the weighing map will be used.');
end
% if only magnitude images are loaded
if ~isempty(magn) && isempty(wmap)
    disp('The normalised RMS in time dimension of magnitude image will be used as the weighting map.');
    tmp     = sqrt(mean(magn.^2,4));
    wmap    = (tmp./max(tmp(:))) .* (mask); 
end
% if nothing is loaded
if ~isempty(magn) && isempty(wmap)
    warning('Providing a weighing map or magnitude images can potentially improve the QSM map quality.');
end

%% Display algorithm parameters
disp('The following parameters are being used...');
disp(['Tolerance    	= ' num2str(tol)]);
disp(['Max. iteration	= ' num2str(maxiter)]);
disp(['Step size        = ' num2str(stepSize)]);

%% main
% NDI default parameters are relative so okay for local field map in ppm
localField = localField/(b0*gyro);

chi = NDI(localField,mask,voxelSize,'b0dir',b0dir,'weight',wmap,...
          'iteration',maxiter,'stepsize',stepSize,'tol',tol);

end

%% set default parameter if not specified
function algorParam2 = check_and_set_algorithm_default(algorParam)

algorParam2 = algorParam;

try algorParam2.qsm.tol         = algorParam.qsm.tol;       catch; algorParam2.qsm.tol      = 1; end
try algorParam2.qsm.stepSize    = algorParam.qsm.stepSize;  catch; algorParam2.qsm.stepSize = 1; end
try algorParam2.qsm.maxiter     = algorParam.qsm.maxiter;   catch; algorParam2.qsm.maxiter  = 200; end

end