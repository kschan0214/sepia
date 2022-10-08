%% get_CIT168_reinf_learn_labels_in_native_space(input,mask,output_dir,algorParam)
%
% Usage:
%
% Input
% --------------
% input         : structure that contains input nifti file names
% mask          : structure that contains mask nifti file names
% output_dir    : directiry where the output files will be stored
% algorParam    : algorithm parameters that are specific to this function
%   .mode       : mode of operation: 1 - perform registration; 
%                                    2 - no regitration (transformation
%                                    files required)
%   .isBiasFieldCorr : correct bias field
%
% Description: Using ANTs to obtain subcortical labels vias atlas
%
% Kwok-shing Chan @ DCCN
% k.chan@donders.ru.nl
% Date created: 6 October 2022
% Date last modified:
%
%
function get_CIT168_reinf_learn_labels_in_native_space(input,output_dir,algorParam)

sepia_universal_variables;

algorParam      = check_and_set_algorithm_default(algorParam);
isBiasFieldCorr = algorParam.isBiasFieldCorr;
mode            = algorParam.mode;

% setup ANTs in environment
setup_ANTs;
% test if antsRegistration exists
[status,~] = system('antsRegistration');
if status == 127
    setenv('PATH', [getenv('PATH') ':' ANTs_HOME]);
end
SEPIA_ANALYSIS_SEGMENTATION_dir = fullfile(SEPIA_HOME,'analysis','segmentation');
shell_script = fullfile(SEPIA_ANALYSIS_SEGMENTATION_dir,'ANTs_get_CIT168_reinf_learn_labels.sh');

output_tmp_dir = fullfile(output_dir,'CIT168_reinf_learn_intermediate_files',filesep);
if ~exist(output_dir,'dir')
    mkdir(output_dir);
    mkdir(output_tmp_dir);
end

switch mode
    case 1      % registration is required
        GRE_nii      = input(1).name;
        GRE_mask_nii = input(2).name;
        T1w_nii      = input(3).name;
        T1w_mask_nii = input(4).name;
        
        cmd = ['sh ' shell_script ' ' num2str(mode) ' ' output_tmp_dir ' ' GRE_nii ' ' GRE_mask_nii ' ' T1w_nii ' ' T1w_mask_nii ' ' num2str(isBiasFieldCorr)];

    case 2      % transformation is provided
        gre_2_t1w_mat                   = input(1).name;
        t1_2_mni2009c_mat               = input(2).name;
        t1_2_mni2009c_inverseWrap_nii   = input(3).name;
        T1w_nii                         = input(4).name;
        
        cmd = ['sh ' shell_script ' ' num2str(mode) ' ' output_tmp_dir ' ' gre_2_t1w_mat ' ' t1_2_mni2009c_mat ' ' t1_2_mni2009c_inverseWrap_nii ' ' T1w_nii];
end

system(cmd);

copyfile(fullfile(output_tmp_dir,'CIT168toMNI152-2009c_det_2gre.nii.gz'),fullfile(output_dir));
copyfile(fullfile(output_tmp_dir,'CIT168toMNI152-2009c_prob_2gre.nii.gz'),fullfile(output_dir));
copyfile(fullfile(output_tmp_dir,'labels.txt'),fullfile(output_dir));

disp('Segmentation is done!')

end

function algorParam2 = check_and_set_algorithm_default(algorParam)

try algorParam2.isBiasFieldCorr = algorParam.isBiasFieldCorr;   catch; algorParam2.isBiasFieldCorr = false; end
try algorParam2.mode = algorParam.mode;                         catch; algorParam2.mode = 1; end

end