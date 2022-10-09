%% [chi] = Wrapper_QSM_TKD(localField,mask,matrixSize,voxelSize,algorParam, headerAndExtraData)
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
% Date created: 8 March 2020
% Date modified: 13 August 2021 (v1.0)
%
%
function [chi] = Wrapper_QSM_TKD(localField,mask,matrixSize,voxelSize,algorParam, headerAndExtraData)
sepia_universal_variables;

% get algorithm parameters
algorParam = check_and_set_algorithm_default(algorParam);
method     = algorParam.qsm.method;
thre_tkd   = algorParam.qsm.threshold;

% get extra data such as magnitude/weights/B0 direction/TE/etc.
headerAndExtraData = check_and_set_SEPIA_header_data(headerAndExtraData);
b0dir = headerAndExtraData.sepia_header.B0_dir;
b0    = headerAndExtraData.sepia_header.B0;

% add path
sepia_addpath;
addpath(fullfile(SEPIA_HOME,'misc','qsm_algorithm','TKD'));

%% Display algorithm parameters
disp('The following parameter is being used...');
disp(['K-space threshold value  = ' num2str(thre_tkd)]);

%% main
chi = qsmTKD(localField,mask,matrixSize,voxelSize,'threshold',thre_tkd,'b0dir',b0dir);
        
% convert from Hz to ppm
chi = chi/(b0*gyro);

end

%% set default parameter if not specified
function algorParam2 = check_and_set_algorithm_default(algorParam)

algorParam2 = algorParam;

try algorParam2.qsm.threshold = algorParam.qsm.threshold; catch; algorParam2.qsm.threshold = 0.15; end

end