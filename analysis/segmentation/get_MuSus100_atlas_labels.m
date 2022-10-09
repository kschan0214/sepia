%% get_MuSus100_atlas_labels(input,mask,output_dir,algorParam)
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
function get_MuSus100_atlas_labels(input,output_dir,algorParam)

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
SEPIA_ANALYSIS_SEGMENTATION_dir = fullfile(SEPIA_HOME,'analysis','segmentation');
SEPIA_ATLAS_dir                 = fullfile(SEPIA_HOME,'atlas');
MuSus100_ATLAS_dir              = fullfile(SEPIA_ATLAS_dir,'MuSus-100_Atlas');

output_tmp_dir = fullfile(output_dir,'MuSus100_intermediate_files',filesep);
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
        
        % TODO: chi only registration
        
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
        
        % Step 3: create hybrid T1w+QSM image on the provided data
        Chi_t1w_nii         = fullfile(output_tmp_dir,[Chi_t1w_nii '_2T1w.nii.gz']);
        img_hybrid_nii      = fullfile(output_tmp_dir,'Hybrid_image.nii.gz');
        template_hybrid_nii = fullfile(MuSus100_ATLAS_dir,'atlas','hybrid.nii.gz');
        
        t1w         = load_nii_img_only(T1w_nii); 
        mask_t1     = load_nii_img_only(T1w_mask_nii);
        Chi_t1w     = load_nii_img_only(Chi_t1w_nii);
        img_hybrid  = compute_hybrid_t1w_chi(t1w, Chi_t1w, mask_t1, 400, 0.5 );
        save_nii_img_only(T1w_nii,img_hybrid_nii, img_hybrid );
        clear img_hybrid Chi_t1w t1w mask_t1
        
        % Step 4: register hybrid image to template space, non linear transform
        shell_script = fullfile(SEPIA_ANALYSIS_SEGMENTATION_dir,'ANTs_t1w_2_t1w_atlas_template.sh');
        cmd = ['sh ' shell_script ' ' output_tmp_dir ' ' img_hybrid_nii ' ' T1w_mask_nii ' ' template_hybrid_nii ' ' num2str(isBiasFieldCorr)];
        system(cmd);

        % Step 5: Apply tranformation
        shell_script = fullfile(SEPIA_ANALYSIS_SEGMENTATION_dir,'ANTs_gre_2_t1wAtlas_applyTransform.sh');
        gre_2_T1w_mat                       = fullfile(output_tmp_dir, 'GRE_2_T1w_0GenericAffine.mat');
        t1_2_t1wTemplate_mat                = fullfile(output_tmp_dir, 'T1w_2_hybrid_0GenericAffine.mat');
        t1_2_t1wTemplate_inverseWrap_nii    = fullfile(output_tmp_dir, 'T1w_2_hybrid_1InverseWarp.nii.gz');

        mode_interp = 1; % 1=GenericLabel; 2=linear
        mode_4D     = 0; % 0=3D; 3=4D; equivalent to option -e in antsApplyTransforms see ANTs doc
        label_nii   = fullfile(MuSus100_ATLAS_dir,'label','mixed.nii.gz');
        cmd = ['sh ' shell_script ' ' output_tmp_dir ' ' num2str(mode_interp) ' ' num2str(mode_4D) ' ' label_nii ' ' Chi_nii ' ' gre_2_T1w_mat ' ' t1_2_t1wTemplate_mat ' ' t1_2_t1wTemplate_inverseWrap_nii];
        system(cmd);

    case 2      % transformation is provided
        gre_2_T1w_mat                       = input.gre2T1wMat;
        t1_2_t1wTemplate_mat                = input.t1w2TemplateMat;
        t1_2_t1wTemplate_inverseWrap_nii    = input.t1w2TemplateiWrap;
        GRE_nii                             = input.gre;
        
        % Step 1: Apply tranformation
        shell_script = fullfile(SEPIA_ANALYSIS_SEGMENTATION_dir,'ANTs_gre_2_t1wAtlas_applyTransform.sh');
        mode_interp = 1; % 1=GenericLabel; 2=linear
        mode_4D     = 0; % 0=3D; 3=4D; equivalent to option -e in antsApplyTransforms see ANTs doc
        label_nii   = fullfile(MuSus100_ATLAS_dir,'label','mixed.nii.gz');
        cmd = ['sh ' shell_script ' ' output_tmp_dir ' ' num2str(mode_interp) ' ' num2str(mode_4D) ' ' label_nii ' ' GRE_nii ' ' gre_2_T1w_mat ' ' t1_2_t1wTemplate_mat ' ' t1_2_t1wTemplate_inverseWrap_nii];
        system(cmd);

end
 
copyfile(fullfile(output_tmp_dir, 'mixed_2gre.nii.gz'),fullfile(output_dir));
copyfile(fullfile(MuSus100_ATLAS_dir,'label','mixed.txt'),fullfile(output_dir,'MuSus100_mixed.txt'));

disp('Segmentation is done!')

end

function algorParam2 = check_and_set_algorithm_default(algorParam)

try algorParam2.isBiasFieldCorr = algorParam.isBiasFieldCorr;   catch; algorParam2.isBiasFieldCorr = false; end
try algorParam2.mode = algorParam.mode;                         catch; algorParam2.mode = 1; end

end