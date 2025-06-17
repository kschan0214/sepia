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
% Date modified: 16 June 2025
%
%
function get_CIT168_reinf_learn_labels(input,output_dir,algorParam)

sepia_universal_variables;

algorParam          = check_and_set_algorithm_default(algorParam);
isBiasFieldCorr     = algorParam.isBiasFieldCorr;
mode                = algorParam.mode;
isAccelerate        = algorParam.isAccelerate;
saveIntermediate    = algorParam.saveIntermediate;

% setup ANTs in environment
SpecifyToolboxesDirectory;
% test if antsRegistration exists
[status,~] = system('antsRegistration'); if status == 127; setenv('PATH', [getenv('PATH') ':' ANTS_HOME]); end
% get atlas directory
SpecifyAtlasDirectory;

template_nii = fullfile(CIT168_reinf_learn_ATLAS_HOME,'MNI152-Nonlin-Asym-2009c','CIT168toMNI152-2009c_T1w_brain.nii.gz');

output_tmp_dir = fullfile(output_dir,'CIT168_reinf_learn_intermediate_files',filesep);
if ~exist(output_dir,'dir');        mkdir(output_dir);      end
if ~exist(output_tmp_dir,'dir');    mkdir(output_tmp_dir);  end

%%%%%%%%%% Handle 4D data %%%%%%%%%%
%if 4D then get 1st echo for registration
img = load_nii_img_only(input.gre);
if ndims(img) == 4      
    img         = img(:,:,:,1); % only first echo
    echo1_fn    = fullfile(output_tmp_dir,'GRE_echo-1.nii.gz');
    save_nii_img_only(input.gre,echo1_fn,img);
    input.gre = echo1_fn;
end
clear img
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

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

        %%%%%%%%%% Speed up using label mask (and ventricle) only %%%%%%%%%%
        mask_regis_nii = fullfile(output_tmp_dir,'mask_registration.nii.gz');
        if isAccelerate
            label_nii   = fullfile(CIT168_reinf_learn_ATLAS_HOME,'MNI152-Nonlin-Asym-2009c','CIT168toMNI152-2009c_det.nii.gz');
            label       = load_nii_img_only(label_nii);
            load(fullfile(SEPIA_ANALYSIS_SEGMENTATION_dir,'private','CIT_mask_ventricle.mat'));
            mask_ventricle = zeros(size(label)); mask_ventricle(idx) = 1;
            mask_regis = imdilate( or(label >0,mask_ventricle), strel('sphere',10));
            save_nii_img_only(label_nii,mask_regis_nii,mask_regis);
            clear mask_regis mask_ventricle label
        end
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        
        % Step 2: register T1w image to template space, non linear transform
        shell_script = fullfile(SEPIA_ANALYSIS_SEGMENTATION_dir,'ANTs_t1w_2_t1w_atlas_template.sh');
        cmd = ['sh ' shell_script ' ' output_tmp_dir ' ' T1w_nii ' ' T1w_mask_nii ' ' template_nii ' ' num2str(isBiasFieldCorr) ' ' mask_regis_nii];
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

%%%%%%%%%% export csv %%%%%%%%%%
if ~isempty(input.chi)

label_mixed     = load_nii_img_only(fullfile(output_dir,'CIT168toMNI152-2009c_det_2gre.nii.gz')); nii = load_untouch_nii(fullfile(output_dir,'CIT168toMNI152-2009c_det_2gre.nii.gz')); % deterministic label
prob_mixed     = load_nii_img_only(fullfile(output_dir,'CIT168toMNI152-2009c_prob_2gre.nii.gz'));
label_unique    = unique(label_mixed(:)); % Proba label 
label_unique    = label_unique(label_unique>0);
chimap          = load_nii_img_only(input.chi);
roi_mean        = nan(numel(label_unique),1);
roi_sd          = nan(numel(label_unique),1);
roi_median      = nan(numel(label_unique),1);
roi_iqr         = nan(numel(label_unique),1);
Nvoxel          = nan(numel(label_unique),1);
roi_wmean       = nan(numel(label_unique),1);
roi_wsd         = nan(numel(label_unique),1);
roi_wvol        = nan(numel(label_unique),1);
for kIdx = 1:numel(label_unique)
    
    label_curr          = label_unique(kIdx);

    if label_curr > 0
        mask_label          = label_mixed == label_curr;
    
        Nvoxel(kIdx)        = numel(mask_label(mask_label>0));
    
        % deterministic mask
        roi_mean(kIdx)      = mean(chimap(mask_label));
        roi_sd(kIdx)        = std(chimap(mask_label));
        roi_median(kIdx)    = median(chimap(mask_label));
        roi_iqr(kIdx)       = iqr(chimap(mask_label));
    
        % probability mask stats
        prob_label          = prob_mixed(:,:,:,kIdx);
        roi_wmean(kIdx)     = sum(chimap(:).*prob_label(:))/sum(prob_label(:));
        roi_wsd(kIdx)       = sqrt(sum(prob_label(:) .* (chimap(:) - roi_wmean(kIdx)).^2) / sum(prob_label(:)));
        roi_wvol(kIdx)      = sum(prob_label(:))*prod(nii.hdr.dime.pixdim(2:4));
    end
end
metric = {'label index','mean (ppm)','sd (ppm)','median (ppm)','iqr (ppm)','no. voxels','weighted mean (ppm)','weighted sd (ppm)','weighted volume (mm3)'};
T = table(label_unique,roi_mean,roi_sd,roi_median,roi_iqr,Nvoxel,roi_wmean,roi_wsd,roi_wvol,'VariableNames',metric);
writetable(T,fullfile(output_dir,'stats_Chimap_CIT168.csv'));

end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if ~saveIntermediate; rmdir(output_tmp_dir,'s'); end

disp('Segmentation is done!')

end

function algorParam2 = check_and_set_algorithm_default(algorParam)

try algorParam2.isBiasFieldCorr     = algorParam.isBiasFieldCorr;       catch; algorParam2.isBiasFieldCorr = false; end
try algorParam2.mode                = algorParam.mode;                  catch; algorParam2.mode = 1; end
try algorParam2.isAccelerate        = algorParam.isAccelerate;      catch; algorParam2.isAccelerate = false; end
try algorParam2.saveIntermediate    = algorParam.saveIntermediate;  catch; algorParam2.saveIntermediate = false; end

end