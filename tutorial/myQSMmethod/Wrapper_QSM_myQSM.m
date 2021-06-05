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
% Date created: 6 April 2020
% Date modified:
%
% You can change the name of the function but DO NOT change the input/output variables
function [chi] = Wrapper_QSM_myQSM(localField,mask,matrixSize,voxelSize,algorParam,headerAndExtraData)
% load some constants 
sepia_universal_variables;

% get algorithm parameters, if user doesn't specify them then set some default values
algorParam = check_and_set_algorithm_default(algorParam);
thre_tkd   = algorParam.qsm.threshold;  % here you can define how SEPIA will store the user input in the 'algorParam' variable

% get extra data such as magnitude/weights/B0 direction/TE/etc.
headerAndExtraData = check_and_set_SEPIA_header_data(headerAndExtraData);
b0dir = headerAndExtraData.b0dir;
b0    = headerAndExtraData.b0;
% magn  = headerAndExtraData.magn;  % you can access the magnitude and/or other data from the 'headerAndExtraData' variable

% add path
sepia_addpath;

%% Display algorithm parameters
disp('The following parameter is being used...');
disp(['K-space threshold value  = ' num2str(thre_tkd)]);

%% main
% you can change the unit before your method if you wish
% localField = localField/(b0*gyro); % convert from Hz to ppm

chi = myQSM(localField,mask,matrixSize,voxelSize,thre_tkd,b0dir);
        
% make sure the output susceptibility map is in 'ppm' which is the default
% unit in SEPIA
chi = chi/(b0*gyro); % convert from Hz to ppm

end

%% set default parameter if not specified
function algorParam2 = check_and_set_algorithm_default(algorParam)

algorParam2 = algorParam;

try algorParam2.qsm.threshold = algorParam.qsm.threshold; catch; algorParam2.qsm.threshold = 0.15; end

end