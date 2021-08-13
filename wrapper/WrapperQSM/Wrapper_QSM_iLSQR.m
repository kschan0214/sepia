%% [chi] = Wrapper_QSM_iLSQR(localField,mask,matrixSize,voxelSize,algorParam, headerAndExtraData)
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
% Description: This is a wrapper function to access iLSQR for SEPIA
%
% Kwok-shing Chan @ DCCN
% k.chan@donders.ru.nl
% Date created: 8 March 2020
% Date modified: 13 August 2021 (v1.0)
%
%
function [chi] = Wrapper_QSM_iLSQR(localField,mask,matrixSize,voxelSize,algorParam, headerAndExtraData)
sepia_universal_variables;

% get algorithm parameters
algorParam  = check_and_set_algorithm_default(algorParam);
method      = algorParam.qsm.method;
tol         = algorParam.qsm.tol;
maxiter     = algorParam.qsm.maxiter;
lambda      = algorParam.qsm.lambda;
optimise    = algorParam.qsm.optimise;

% get extra data such as magnitude/weights/B0 direction/TE/etc.
headerAndExtraData = check_and_set_SEPIA_header_data(headerAndExtraData);
b0dir       = headerAndExtraData.sepia_header.B0_dir;
b0          = headerAndExtraData.sepia_header.B0;
wmap        = get_variable_from_headerAndExtraData(headerAndExtraData, 'weights', matrixSize);
initGuess   = get_variable_from_headerAndExtraData(headerAndExtraData, 'initGuess', matrixSize);

if isempty(wmap)
    wmap = ones(matrixSize);
end
if isempty(initGuess)
    initGuess = zeros(matrixSize);
end

% add path
sepia_addpath;
addpath(fullfile(SEPIA_HOME,'misc','qsm_algorithm','closedFormL2'));
addpath(fullfile(SEPIA_HOME,'misc','qsm_algorithm','iLSQR_qsmhub'));

%% Display algorithm parameters
disp('The following parameters are being used...');
disp(['Tolerance                = ' num2str(tol)]);
disp(['Max. iteration           = ' num2str(maxiter)]);
disp(['Regularisation lambda	= ' num2str(lambda)]);
disp(['Self-optimisation?       = ' num2str(optimise)]);

%% main
chi = qsmIterativeLSQR(localField,mask,matrixSize,voxelSize,...
            'lambda',lambda,'tol',tol,'iteration',maxiter,'weight',wmap,...
            'initGuess',initGuess,'optimise',optimise,'b0dir',b0dir);
        
% convert from Hz to ppm
chi = chi/(b0*gyro);

end

%% set default parameter if not specified
function algorParam2 = check_and_set_algorithm_default(algorParam)

algorParam2 = algorParam;

try algorParam2.qsm.tol         = algorParam.qsm.tol;       catch; algorParam2.qsm.tol      = 1e-3; end
try algorParam2.qsm.lambda      = algorParam.qsm.lambda;    catch; algorParam2.qsm.lambda   = 1e-1; end
try algorParam2.qsm.maxiter     = algorParam.qsm.maxiter;   catch; algorParam2.qsm.maxiter  = 50; end
try algorParam2.qsm.optimise    = algorParam.qsm.optimise;  catch; algorParam2.qsm.optimise = false; end

end