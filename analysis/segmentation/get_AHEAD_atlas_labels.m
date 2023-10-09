%% get_AHEAD_atlas_labels(input,mask,output_dir,algorParam)
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
function get_AHEAD_atlas_labels(input,output_dir,algorParam)

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
% AHEAD_ATLAS_HOME                = fullfile(SEPIA_ATLAS_dir,'AHEAD_atlas');

output_tmp_dir = fullfile(output_dir,'AHEAD_intermediate_files',filesep);
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
        Chi_nii      = input.chi;
        
        % get Chimap basename (just in case is compressed file)
        [~,Chi_t1w_nii,~] = fileparts(Chi_nii);
        [~,Chi_t1w_nii,~] = fileparts(Chi_t1w_nii);
        
        % Step 1: GRE to T1w     
        shell_script    = fullfile(SEPIA_ANALYSIS_SEGMENTATION_dir,'ANTs_gre_2_t1w.sh');
        cmd             = ['sh ' shell_script ' ' output_tmp_dir ' ' GRE_nii ' ' GRE_mask_nii ' ' T1w_nii ' ' T1w_mask_nii ' ' num2str(isBiasFieldCorr)];
        system(cmd);
        isBiasFieldCorr = 0;    % this should be done after the previous line 

        % Step 2: Bring Chimap to T1w space
        shell_script    = fullfile(SEPIA_ANALYSIS_SEGMENTATION_dir,'ANTs_gre_2_t1w_applyTransform.sh');
        cmd             = ['sh ' shell_script ' ' output_tmp_dir ' ' Chi_nii ' ' T1w_nii];
        system(cmd);
        
        % Step 3.1: create hybrid T1w+QSM image on the provided data
        Chi_t1w_nii         = fullfile(output_tmp_dir,[Chi_t1w_nii '_2T1w.nii.gz']);
        img_hybrid_nii      = fullfile(output_tmp_dir,'Hybrid_image.nii.gz');
        template_hybrid_nii = fullfile(output_tmp_dir,'Template_hybrid_image.nii.gz');
        
        t1w         = load_nii_img_only(T1w_nii); 
        mask_t1     = load_nii_img_only(T1w_mask_nii);
        Chi_t1w     = load_nii_img_only(Chi_t1w_nii);
        img_hybrid  = compute_hybrid_t1w_chi(t1w, Chi_t1w, mask_t1, 400, 0.5 );
        save_nii_img_only(T1w_nii,img_hybrid_nii, img_hybrid );

        % create hybrid T1w+QSM image on the atlas data
        t1w         = load_nii_img_only(fullfile(AHEAD_ATLAS_HOME,'Templates_mni09b','ahead_final_med_r1map_n104_mni09b.nii.gz'));
        Chi_t1w     = load_nii_img_only(fullfile(AHEAD_ATLAS_HOME,'Templates_mni09b','ahead_final_med_qsm_n104_mni09b.nii.gz'));
        mask_t1     = t1w ~= 0;
        img_hybrid  = compute_hybrid_t1w_chi(t1w, Chi_t1w, mask_t1, 400, 0.5 );
        save_nii_img_only(fullfile(AHEAD_ATLAS_HOME,'Templates_mni09b','ahead_final_med_r1map_n104_mni09b.nii.gz'),template_hybrid_nii, img_hybrid);
        clear img_hybrid Chi_t1w t1w mask_t1
        
        % Step 4: register hybrid image to template space, non linear transform
        shell_script = fullfile(SEPIA_ANALYSIS_SEGMENTATION_dir,'ANTs_t1w_2_t1w_atlas_template.sh');
        cmd = ['sh ' shell_script ' ' output_tmp_dir ' ' img_hybrid_nii ' ' T1w_mask_nii ' ' template_hybrid_nii ' ' num2str(isBiasFieldCorr)];
        system(cmd);

        % Step 5: Apply tranformation
        shell_script = fullfile(SEPIA_ANALYSIS_SEGMENTATION_dir,'ANTs_gre_2_t1wAtlas_applyTransform.sh');
        gre_2_T1w_mat                       = fullfile(output_tmp_dir, 'GRE_2_T1w_0GenericAffine.mat');
        t1_2_t1wTemplate_mat                = fullfile(output_tmp_dir, 'T1w_2_Template_hybrid_image_0GenericAffine.mat');
        t1_2_t1wTemplate_inverseWrap_nii    = fullfile(output_tmp_dir, 'T1w_2_Template_hybrid_image_1InverseWarp.nii.gz');
        
        mask_list = dir(fullfile(AHEAD_ATLAS_HOME,'structures_mni09b','*_mask-*nii*'));
        mode_interp = 1; % 1=GenericLabel; 2=linear
        mode_4D     = 0; % 0=3D; 3=4D; equivalent to option -e in antsApplyTransforms see ANTs doc
        for k = 1:numel(mask_list)
            label_nii   = fullfile(AHEAD_ATLAS_HOME,'structures_mni09b',mask_list(k).name);
            cmd = ['sh ' shell_script ' ' output_tmp_dir ' ' num2str(mode_interp) ' ' num2str(mode_4D) ' ' label_nii ' ' Chi_nii ' ' gre_2_T1w_mat ' ' t1_2_t1wTemplate_mat ' ' t1_2_t1wTemplate_inverseWrap_nii];
            system(cmd);
        end

        prob_list = dir(fullfile(AHEAD_ATLAS_HOME,'structures_mni09b','*_proba-*nii*'));
        mode_interp = 2; % 1=GenericLabel; 2=linear
        mode_4D     = 0; % 0=3D; 3=4D; equivalent to option -e in antsApplyTransforms see ANTs doc
        for k = 1:numel(prob_list)
            label_nii   = fullfile(AHEAD_ATLAS_HOME,'structures_mni09b',prob_list(k).name);
            cmd = ['sh ' shell_script ' ' output_tmp_dir ' ' num2str(mode_interp) ' ' num2str(mode_4D) ' ' label_nii ' ' Chi_nii ' ' gre_2_T1w_mat ' ' t1_2_t1wTemplate_mat ' ' t1_2_t1wTemplate_inverseWrap_nii];
            system(cmd);
        end

    case 2      % transformation is provided
        gre_2_T1w_mat                       = input.gre2T1wMat;
        t1_2_t1wTemplate_mat                = input.t1w2TemplateMat;
        t1_2_t1wTemplate_inverseWrap_nii    = input.t1w2TemplateiWrap;
        GRE_nii                             = input.gre;
        
        % Step 1: Apply tranformation
        shell_script = fullfile(SEPIA_ANALYSIS_SEGMENTATION_dir,'ANTs_gre_2_t1wAtlas_applyTransform.sh');

        mask_list = dir(fullfile(AHEAD_ATLAS_HOME,'structures_mni09b','*_mask-*nii*'));
        mode_interp = 1; % 1=GenericLabel; 2=linear
        mode_4D     = 0; % 0=3D; 3=4D; equivalent to option -e in antsApplyTransforms see ANTs doc
        for k = 1:numel(mask_list)
            label_nii   = fullfile(AHEAD_ATLAS_HOME,'structures_mni09b',mask_list(k).name);
            cmd = ['sh ' shell_script ' ' output_tmp_dir ' ' num2str(mode_interp) ' ' num2str(mode_4D) ' ' label_nii ' ' GRE_nii ' ' gre_2_T1w_mat ' ' t1_2_t1wTemplate_mat ' ' t1_2_t1wTemplate_inverseWrap_nii];
            system(cmd);
        end

        prob_list = dir(fullfile(AHEAD_ATLAS_HOME,'structures_mni09b','*_proba-*nii*'));
        mode_4D     = 0; % 0=3D; 3=4D; equivalent to option -e in antsApplyTransforms see ANTs doc
        mode_interp = 2; % 1=GenericLabel; 2=linear
        for k = 1:numel(prob_list)
            label_nii   = fullfile(AHEAD_ATLAS_HOME,'structures_mni09b',prob_list(k).name);
            cmd = ['sh ' shell_script ' ' output_tmp_dir ' ' num2str(mode_interp) ' ' num2str(mode_4D) ' ' label_nii ' ' GRE_nii ' ' gre_2_T1w_mat ' ' t1_2_t1wTemplate_mat ' ' t1_2_t1wTemplate_inverseWrap_nii];
            system(cmd);
        end

end


copyfile(fullfile(output_tmp_dir, '*_mask-*.nii*'),fullfile(output_dir));
copyfile(fullfile(output_tmp_dir, '*_proba-*.nii*'),fullfile(output_dir));

disp('Segmentation is done!')

end

function algorParam2 = check_and_set_algorithm_default(algorParam)

try algorParam2.isBiasFieldCorr = algorParam.isBiasFieldCorr;   catch; algorParam2.isBiasFieldCorr = false; end
try algorParam2.mode = algorParam.mode;                         catch; algorParam2.mode = 1; end

end