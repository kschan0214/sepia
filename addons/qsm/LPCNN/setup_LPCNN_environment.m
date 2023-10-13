%% This file specifies the Python environment and Checkpoint of QSMnet+
% qsmnet_dir          = '/project/3015069.05/bids/code/QSMnet/';
% dir_net             = fullfile(qsmnet_dir,'Checkpoints/');

% Specify LPCNN_HOME corresponding to the code from Github
LPCNN_HOME          = '/home/common/matlab/sepia/external/LPCNN/LPCNN_v1.0';
% Specify the Python environment that has QSMnet+ installed
python_interpreter  = '/home/common/matlab/sepia/external/LPCNN/LPCNN_v1.0/lpcnn-env/bin/python';
% Specify the directory that contains the training parameters
checkpoint_fn       = fullfile(LPCNN_HOME,'checkpoints','lpcnn_test_Bmodel.pkl');
% checkpoint_fn       = '/project/3015069.05/bids/code/LPCNN/checkpoints/lpcnn_test_Bmodel.pkl';

