%% get_CIT168_reinf_learn_labels(input,mask,output_dir,algorParam)
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
% Date modified: 9 October 2023
%
%
function get_CIT168_reinf_learn_labels(input,output_dir,algorParam)

sepia_universal_variables;

algorParam      = check_and_set_algorithm_default(algorParam);
isBiasFieldCorr = algorParam.isBiasFieldCorr;
mode            = algorParam.mode;

% setup ANTs in environment
SpecifyToolboxesDirectory;
% test if antsRegistration exists
[status,~] = system('antsRegistration');
if status == 127
    setenv('PATH', [getenv('PATH') ':' ANTS_HOME]);
end
% get atlas directory
SpecifyAtlasDirectory;
% SEPIA_ANALYSIS_SEGMENTATION_dir = fullfile(SEPIA_HOME,'analysis','segmentation');
% SEPIA_ATLAS_dir                 = fullfile(SEPIA_HOME,'atlas');
% CIT168_reinf_learn_ATLAS_HOME   = fullfile(SEPIA_ATLAS_dir,'CIT168_Reinf_Learn_v1.1.0');

template_nii = fullfile(CIT168_reinf_learn_ATLAS_HOME,'MNI152-Nonlin-Asym-2009c','CIT168toMNI152-2009c_T1w_brain.nii.gz');

output_tmp_dir = fullfile(output_dir,'CIT168_reinf_learn_intermediate_files',filesep);
if ~exist(output_dir,'dir')
    mkdir(output_dir);
    mkdir(output_tmp_dir);
end

switch mode
    case 1      % registration is required
        GRE_nii      = input.gre;
        GRE_mask_nii = input.greMask;
        T1w_nii      = input.t1w;
        T1w_mask_nii = input.t1wMask;
        
        % Step 1: GRE to T1w     
        shell_script    = fullfile(SEPIA_ANALYSIS_SEGMENTATION_dir,'ANTs_gre_2_t1w.sh');
        cmd             = ['sh ' shell_script ' ' output_tmp_dir ' ' GRE_nii ' ' GRE_mask_nii ' ' T1w_nii ' ' T1w_mask_nii ' ' num2str(isBiasFieldCorr)];
        system(cmd);
        isBiasFieldCorr = 0;    % this should be done after the previous line 
        
        % Step 2: register T1w image to template space, non linear transform
        shell_script = fullfile(SEPIA_ANALYSIS_SEGMENTATION_dir,'ANTs_t1w_2_t1w_atlas_template.sh');
        cmd = ['sh ' shell_script ' ' output_tmp_dir ' ' T1w_nii ' ' T1w_mask_nii ' ' template_nii ' ' num2str(isBiasFieldCorr)];
        system(cmd);

        % Step 3: Apply tranformation
        shell_script                        = fullfile(SEPIA_ANALYSIS_SEGMENTATION_dir,'ANTs_gre_2_t1wAtlas_applyTransform.sh');
        gre_2_T1w_mat                       = fullfile(output_tmp_dir, 'GRE_2_T1w_0GenericAffine.mat');
        t1_2_t1wTemplate_mat                = fullfile(output_tmp_dir, 'T1w_2_CIT168toMNI152-2009c_T1w_brain_0GenericAffine.mat');
        t1_2_t1wTemplate_inverseWrap_nii    = fullfile(output_tmp_dir, 'T1w_2_CIT168toMNI152-2009c_T1w_brain_1InverseWarp.nii.gz');
        
        % deterministic labels
        mode_interp = 1; % 1=GenericLabel; 2=linear
        mode_4D     = 0; % 0=3D; 3=4D; equivalent to option -e in antsApplyTransforms see ANTs doc
        label_nii   = fullfile(CIT168_reinf_learn_ATLAS_HOME,'MNI152-Nonlin-Asym-2009c','CIT168toMNI152-2009c_det.nii.gz');
        cmd = ['sh ' shell_script ' ' output_tmp_dir ' ' num2str(mode_interp) ' ' num2str(mode_4D) ' ' label_nii ' ' GRE_nii ' ' gre_2_T1w_mat ' ' t1_2_t1wTemplate_mat ' ' t1_2_t1wTemplate_inverseWrap_nii];
        system(cmd);
        
        % probabilistic labels
        mode_interp = 2; % 1=GenericLabel; 2=linear
        mode_4D     = 3; % 0=3D; 3=4D; equivalent to option -e in antsApplyTransforms see ANTs doc
        label_nii   = fullfile(CIT168_reinf_learn_ATLAS_HOME,'MNI152-Nonlin-Asym-2009c','CIT168toMNI152-2009c_prob.nii.gz');
        cmd = ['sh ' shell_script ' ' output_tmp_dir ' ' num2str(mode_interp) ' ' num2str(mode_4D) ' ' label_nii ' ' GRE_nii ' ' gre_2_T1w_mat ' ' t1_2_t1wTemplate_mat ' ' t1_2_t1wTemplate_inverseWrap_nii];
        system(cmd);

    case 2      % transformation is provided
        gre_2_T1w_mat                       = input.gre2T1wMat;
        t1_2_t1wTemplate_mat                = input.t1w2TemplateMat;
        t1_2_t1wTemplate_inverseWrap_nii    = input.t1w2TemplateiWrap;
        GRE_nii                             = input.gre;
        
        % Step 1: Apply tranformation
        shell_script                        = fullfile(SEPIA_ANALYSIS_SEGMENTATION_dir,'ANTs_gre_2_t1wAtlas_applyTransform.sh');
        % deterministic labels
        mode_interp = 1; % 1=GenericLabel; 2=linear
        mode_4D     = 0; % 0=3D; 3=4D; equivalent to option -e in antsApplyTransforms see ANTs doc
        label_nii   = fullfile(CIT168_reinf_learn_ATLAS_HOME,'MNI152-Nonlin-Asym-2009c','CIT168toMNI152-2009c_det.nii.gz');
        cmd = ['sh ' shell_script ' ' output_tmp_dir ' ' num2str(mode_interp) ' ' num2str(mode_4D) ' ' label_nii ' ' GRE_nii ' ' gre_2_T1w_mat ' ' t1_2_t1wTemplate_mat ' ' t1_2_t1wTemplate_inverseWrap_nii];
        system(cmd);
        
        % probabilistic labels
        mode_interp = 2; % 1=GenericLabel; 2=linear
        mode_4D     = 3; % 0=3D; 3=4D; equivalent to option -e in antsApplyTransforms see ANTs doc
        label_nii   = fullfile(CIT168_reinf_learn_ATLAS_HOME,'MNI152-Nonlin-Asym-2009c','CIT168toMNI152-2009c_prob.nii.gz');
        cmd = ['sh ' shell_script ' ' output_tmp_dir ' ' num2str(mode_interp) ' ' num2str(mode_4D) ' ' label_nii ' ' GRE_nii ' ' gre_2_T1w_mat ' ' t1_2_t1wTemplate_mat ' ' t1_2_t1wTemplate_inverseWrap_nii];
        system(cmd);

end
 
copyfile(fullfile(output_tmp_dir,'CIT168toMNI152-2009c_det_2gre.nii.gz'),fullfile(output_dir));
copyfile(fullfile(output_tmp_dir,'CIT168toMNI152-2009c_prob_2gre.nii.gz'),fullfile(output_dir));
copyfile(fullfile(CIT168_reinf_learn_ATLAS_HOME,'labels.txt'),fullfile(output_dir,'CIT168_reinf_learn_labels.txt'));

disp('Segmentation is done!')

end

function algorParam2 = check_and_set_algorithm_default(algorParam)

try algorParam2.isBiasFieldCorr = algorParam.isBiasFieldCorr;   catch; algorParam2.isBiasFieldCorr = false; end
try algorParam2.mode = algorParam.mode;                         catch; algorParam2.mode = 1; end

end