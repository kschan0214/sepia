%% [chi] = Wrapper_QSM_xQSM(localField,mask,matrixSize,voxelSize,algorParam, headerAndExtraData)
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
% Date created: 24 September 2022
% Date modified: 
%
%
function [chi] = Wrapper_QSM_xQSM(localField,mask,matrixSize,voxelSize,algorParam, headerAndExtraData)
sepia_universal_variables;

% get algorithm parameters, if user doesn't specify them then set some default values
algorParam      = check_and_set_algorithm_default(algorParam);
recon_method	= algorParam.qsm.solver;

% get extra data such as magnitude/weights/B0 direction/TE/etc.
headerAndExtraData = check_and_set_SEPIA_header_data(headerAndExtraData);
% b0dir = headerAndExtraData.sepia_header.B0_dir;
b0    = headerAndExtraData.sepia_header.B0;

% add path
sepia_addpath;
setup_xQSM_environment;
addpath(fullfile(deepMRI_HOME,'xQSM','matlab','eval'));
addpath(checkpoint_dir)

%% preparation
% xQSM works in ppm
localField = localField/(b0*gyro);

% note the size of the field map input needs to be divisibel by 8
% otherwise 0 padding should be done first
imSize = size(localField);
if any(mod(imSize, 8))
    [localField, pos] = ZeroPadding(field, 8);
end

%% Display algorithm parameters
% recon_methods_list = {'xQSM_invivo', 'xQSM_syn', 'xQSM_invivo_withNoiseLayer', 'Unet_invivo', 'Unet_syn'};
fprintf('Reconstructing QSM using %s\n', recon_method);

% disp('The following parameters are being used...');

%% main
if canUseGPU()
    % (1) if your MATLAB is configured with CUDA GPU acceleration
    chi = Eval(localField, recon_method, 'gpu');
else
    % (2) otherwise if CUDA is not available, use CPU instead, this is much slower
    chi = Eval(localField, recon_method, 'cpu');
end

% if zeropadding was performed, then do zero-removing before next step;
if any(mod(imSize, 8))
    chi = ZeroRemoving(chi, pos);
end

chi = chi .* mask;

end

%% set default parameter if not specified
function algorParam2 = check_and_set_algorithm_default(algorParam)

algorParam2 = algorParam;
% 
try algorParam2.qsm.solver         = algorParam.qsm.solver;       catch; algorParam2.qsm.solver      = 'xQSM_invivo_withNoiseLayer'; end
% try algorParam2.qsm.stepSize    = algorParam.qsm.stepSize;  catch; algorParam2.qsm.stepSize = 1; end
% try algorParam2.qsm.maxiter     = algorParam.qsm.maxiter;   catch; algorParam2.qsm.maxiter  = 200; end

end