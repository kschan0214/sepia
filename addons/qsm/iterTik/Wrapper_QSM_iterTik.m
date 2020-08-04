%% [chi] = Wrapper_QSM_myQSM(localField,mask,matrixSize,voxelSize,algorParam, headerAndExtraData)
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
% Description: This is a wrapper function to access TKD for SEPIA
%
% Kwok-shing Chan @ DCCN
% k.chan@donders.ru.nl
% Date created: 4 August 2020
% Date modified:
%
% You can change the name of the function but DO NOT change the input/output variables
function [chi] = Wrapper_QSM_iterTik(localField,mask,matrixSize,voxelSize,algorParam,headerAndExtraData)
% load some constants 
sepia_universal_variables;

% get algorithm parameters, if user doesn't specify them then set some default values
algorParam = check_and_set_algorithm_default(algorParam);
solver     = algorParam.qsm.solver;
thre_tkd   = algorParam.qsm.threshold;  % here you can define how SEPIA will store the user input in the 'algorParam' variable
alpha      = algorParam.qsm.lambda;
tolerance  = algorParam.qsm.tolerance;

% get extra data such as magnitude/weights/B0 direction/TE/etc.
headerAndExtraData = check_and_set_SEPIA_header_data(headerAndExtraData);
b0dir = headerAndExtraData.b0dir;
b0    = headerAndExtraData.b0;
weights = headerAndExtraData.weights;
% magn  = headerAndExtraData.magn;  % you can access the magnitude and/or other data from the 'headerAndExtraData' variable

% add path
% sepia_addpath(fullfile(SEPIA_HOME,'external','MRI_susceptibility_calculation','MATLAB'));
addpath(fullfile(SEPIA_HOME,'external','MRI_susceptibility_calculation','MATLAB'))

%% Display algorithm parameters + main
Parameters.FieldMap     = localField/(b0*gyro); % in ppm!
Parameters.Mask         = mask; % default = 1
Parameters.Resolution   = voxelSize; % default = [1,1,1]
Parameters.B0direction  = b0dir; % default = [0,0,1] 

disp('The following parameter is being used...');
switch solver
    case 'Truncated kspace division'
        disp(['K-space threshold value  = ' num2str(thre_tkd)]);
        
        Parameters.Threshold    = thre_tkd; % delfault = 2/3
        
        chi = TKD(Parameters);
        
    case 'Direct Tikhonov'
        disp(['Regularisation value  = ' num2str(alpha)]);
        
        Parameters.Alpha = alpha; % delfault = 0.05
        
        chi = dirTik(Parameters);
        
    case 'Iterative Tikhonov'
        disp(['Regularisation value  = ' num2str(alpha)]);
        disp(['Conjugate gradient stopping threshold  = ' num2str(tolerance)]);
        
        Parameters.Noise = 1./weights; % up to a scaling factor
        Parameters.Alpha = alpha; % delfault = 0.05
        Parameters.StoppingThreshold = 0.03; % delfault = 0.03

        chi = iterTik(Parameters);
        
end

end

%% set default parameter if not specified
function algorParam2 = check_and_set_algorithm_default(algorParam)

algorParam2 = algorParam;

try algorParam2.qsm.solver    = algorParam.qsm.solver;    catch; algorParam2.qsm.solver = 'Truncated kspace division'; end
try algorParam2.qsm.threshold = algorParam.qsm.threshold; catch; algorParam2.qsm.threshold = 2/3; end
try algorParam2.qsm.lambda    = algorParam.qsm.lambda;    catch; algorParam2.qsm.lambda = 0.05; end
try algorParam2.qsm.tolerance = algorParam.qsm.tolerance; catch; algorParam2.qsm.tolerance = 0.03; end

end