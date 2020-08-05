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
function [RDF] = Wrapper_BFR_iRSHARP(totalField,mask,matrixSize,voxelSize,algorParam, headerAndExtraData)
% load some constants 
sepia_universal_variables;

% get algorithm parameters
algorParam  = check_and_set_algorithm_default(algorParam);
method      = algorParam.bfr.method;
radius     	= algorParam.bfr.radius;
thresh_tsvd = algorParam.bfr.threshold;
iRSHARP_C   = algorParam.bfr.constant;

% get extra data such as magnitude/weights/B0 direction/TE/etc.
headerAndExtraData = check_and_set_SEPIA_header_data(headerAndExtraData);
b0    = headerAndExtraData.b0;
phase = headerAndExtraData.phase;

% add path
% sepia_addpath(fullfile(SEPIA_HOME,'external','MRI_susceptibility_calculation','MATLAB'));
addpath(genpath(fullfile(SEPIA_HOME,'external','JHUKKI_QSM_Toolbox','JHUKKI_QSM_Toolbox_v3','JHUKKI_QSM_Toolbox_v3.0')));

%% Display algorithm parameters
disp('The following parameters are being used...');
disp(['Radius(voxel)            = ' num2str(radius)]);
disp(['Truncated SVD threshold  = ' num2str(thresh_tsvd)]);
disp(['Adjusting constant       = ' num2str(iRSHARP_C)]);

%% main
Params.iRSHARP_C    = iRSHARP_C;    % add in Params to pass to iRHSARPv1 
Params.SHARPradius  = radius;
Params.thresh_tsvd  = thresh_tsvd;
Params.sizeVol      = matrixSize;
Params.fov          = matrixSize .* voxelSize;
Params.voxSize      = voxelSize;
% Params.B0           = b0;
% Params.gamma        = gyro*1e6;
Params.echoNums     = 1;

RDF         = iRSHARPv1(totalField, phase, mask, Params);

multiWaitbar('CloseAll'); 

end

%% set default parameter if not specified
function algorParam2 = check_and_set_algorithm_default(algorParam)

algorParam2 = algorParam;

try algorParam2.bfr.radius      = algorParam.bfr.radius;  	catch; algorParam2.bfr.radius       = 8;  end
try algorParam2.bfr.threshold   = algorParam.bfr.threshold;	catch; algorParam2.bfr.threshold	= 0.05;   end
try algorParam2.bfr.constant    = algorParam.bfr.constant;	catch; algorParam2.bfr.constant     = 0.25;   end

end