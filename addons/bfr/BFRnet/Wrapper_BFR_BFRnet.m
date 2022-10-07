%% [RDF] = Wrapper_BFR_BFRnet(totalField,mask,~,~,~, headerAndExtraData)
%
% Input
% --------------
% totalField    : total field map (background + tissue fields), in Hz
% mask          : signal mask
% matrixSize    : size of the input image
% voxelSize     : spatial resolution of each dimension of the data, in mm
% algorParam    : structure contains fields with algorithm-specific parameter(s)
% headerAndExtraData : structure contains extra header info/data for the algorithm
%
% Output
% --------------
% RDF           : local field map
%
% Description: This is a wrapper function to access BFRnet in SEPIA
%
% Kwok-shing Chan @ DCCN
% k.chan@donders.ru.nl
% Date created: 24 September 2022
% Date modified: 
%
%
function [RDF] = Wrapper_BFR_BFRnet(totalField,mask,matrixSize,voxelSize,algorParam, headerAndExtraData)
sepia_universal_variables;

% get algorithm parameters, if user doesn't specify them then set some default values
% algorParam      = check_and_set_algorithm_default(algorParam);
% recon_method	= algorParam.qsm.solver;

% get extra data such as magnitude/weights/B0 direction/TE/etc.
headerAndExtraData = check_and_set_SEPIA_header_data(headerAndExtraData);
% b0dir = headerAndExtraData.sepia_header.B0_dir;
b0    = headerAndExtraData.sepia_header.B0;

% add path
sepia_addpath;
setup_BFRnet_environment;
addpath(fullfile(deepMRI_HOME,'BFRnet','Eval'));

%% preparation
% BFRnet works in ppm
totalField = double(totalField/(b0*gyro) .* mask);

% note the size of the field map input needs to be divisibel by 8
% otherwise 0 padding should be done first
imSize = size(totalField);
if any(mod(imSize, 8))
    [totalField, pos] = ZeroPadding(totalField, 8);
end

%% Display algorithm parameters

%% main
% Load the BFRnet and process reconstruction
if canUseGPU()
    [bkg] = MyPredictGPU(totalField, checkpoints); % Recon using GPU
else
    [bkg] = MyPredictCPU(totalField, checkpoints); % Recon using CPU
end

% bkg = double(bkg .* mask); % The reconstructio  result is background field map.
RDF = totalField - bkg;

% if zeropadding was performed, then do zero-removing before next step;
if any(mod(imSize, 8))
    RDF = ZeroRemoving(RDF, pos);
end

RDF = RDF .* mask;

RDF = RDF * b0 * gyro; % ppm to Hz

end

% %% set default parameter if not specified
% function algorParam2 = check_and_set_algorithm_default(algorParam)
% 
% algorParam2 = algorParam;
% % 
% try algorParam2.qsm.solver         = algorParam.qsm.solver;       catch; algorParam2.qsm.solver      = 'xQSM_invivo_withNoiseLayer'; end
% % try algorParam2.qsm.stepSize    = algorParam.qsm.stepSize;  catch; algorParam2.qsm.stepSize = 1; end
% % try algorParam2.qsm.maxiter     = algorParam.qsm.maxiter;   catch; algorParam2.qsm.maxiter  = 200; end
% 
% end