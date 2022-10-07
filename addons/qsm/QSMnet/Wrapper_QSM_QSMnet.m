%% [chi] = Wrapper_QSM_QSMnet(localField,mask,matrixSize,voxelSize,algorParam, headerAndExtraData)
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
function [chi] = Wrapper_QSM_QSMnet(localField,mask,matrixSize,voxelSize,algorParam, headerAndExtraData)
sepia_universal_variables;

% get algorithm parameters
% algorParam  = check_and_set_algorithm_default(algorParam);
% method      = algorParam.qsm.method;

% get extra data such as magnitude/weights/B0 direction/TE/etc.
headerAndExtraData = check_and_set_SEPIA_header_data(headerAndExtraData);
% b0dir = headerAndExtraData.sepia_header.B0_dir;
b0    = headerAndExtraData.sepia_header.B0;

addon_qsmnet_dir = fileparts(mfilename('fullpath'));

% add path
sepia_addpath;

setup_qsmnet_environment;
inference_template  = fullfile(addon_qsmnet_dir, 'inference_template.py');

temp_dir        = fullfile(pwd,'temp_qsmnet_dir');
temp_output_dir = fullfile(temp_dir, 'Prediction');
temp_INPUT      = fullfile(temp_dir, 'tissue_phase_qsmnet.mat');
inference_fn    = fullfile(temp_dir, 'inference.py');

mkdir(temp_dir)
mkdir(temp_output_dir)
%% preparation
% TODO: up-/down-sample inut data to match 1 mm isotropic resolution
% TODO: rotate image such that B0_dir = [0;0;1]
% clear GPU
gpuDevice([]);

%% main
% QSMnet+ works with ppm
phs_tissue = localField/(b0*gyro);

save(temp_INPUT,'phs_tissue');

copyfile(inference_template,inference_fn);
% should work for at least unix system
system(sprintf('sed -i "2i QSMnet_HOME = ''%s''" %s',       QSMnet_HOME,inference_fn));
system(sprintf('sed -i "35i dir_net = ''%s''" %s',          dir_net,inference_fn));
system(sprintf('sed -i "40i FILE_PATH_INPUT = ''%s''" %s',  temp_INPUT,inference_fn));
system(sprintf('sed -i "41i FILE_PATH_PRED = ''%s''" %s',   temp_output_dir,inference_fn));

% run QSMnet
qsmnet_cmd = [python_interpreter ' ' inference_fn];
system(qsmnet_cmd);

load(fullfile(temp_output_dir,'subject_QSMnet+_64_25.mat'),'sus');

chi = sus .* mask;

rmdir(temp_dir, 's')
rmdir(fullfile(pwd,'..','Checkpoints'), 's')

end

% %% set default parameter if not specified
% function algorParam2 = check_and_set_algorithm_default(algorParam)
% 
% algorParam2 = algorParam;
% % 
% % try algorParam2.qsm.tol         = algorParam.qsm.tol;       catch; algorParam2.qsm.tol      = 1; end
% % try algorParam2.qsm.stepSize    = algorParam.qsm.stepSize;  catch; algorParam2.qsm.stepSize = 1; end
% % try algorParam2.qsm.maxiter     = algorParam.qsm.maxiter;   catch; algorParam2.qsm.maxiter  = 200; end
% 
% end