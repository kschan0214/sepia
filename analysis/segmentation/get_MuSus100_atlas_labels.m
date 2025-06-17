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
% Date modified: 9 October 2023
% Date modified: 16 June 2025
%
%
function get_MuSus100_atlas_labels(input,output_dir,algorParam)

sepia_universal_variables;

algorParam          = check_and_set_algorithm_default(algorParam);
isBiasFieldCorr     = algorParam.isBiasFieldCorr;
isAccelerate        = algorParam.isAccelerate;
saveIntermediate    = algorParam.saveIntermediate;
mode                = algorParam.mode;
isAutoMatchContrast = algorParam.isAutoMatchContrast;

% setup ANTs in environment
SpecifyToolboxesDirectory;
% test if antsRegistration exists
[status,~] = system('antsRegistration');if status == 127; setenv('PATH', [getenv('PATH') ':' ANTS_HOME]); end
% get atlas directory
SpecifyAtlasDirectory;

output_tmp_dir = fullfile(output_dir,'MuSus100_intermediate_files',filesep);
if ~exist(output_dir,'dir');        mkdir(output_dir);      end
if ~exist(output_tmp_dir,'dir');    mkdir(output_tmp_dir);  end

%%%%%%%%%% Handle 4D data %%%%%%%%%%
% %if 4D then get 1st echo for registration
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
        Chi_nii      = input.chi;
        
        % TODO: chi only registration
        
        % get Chimap basename (just in case is compressed file)
        Chi_t1w_nii = get_nifti_filename(Chi_nii);
        
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
        img_hybrid_nii      = fullfile(output_tmp_dir,'Hybrid_image.nii.gz');
        template_hybrid_nii = fullfile(MuSus100_ATLAS_HOME,'atlas','hybrid.nii.gz');
        Chi_t1w_nii         = fullfile(output_tmp_dir,[Chi_t1w_nii '_2T1w.nii.gz']);

        % compute hybrid image based on predefined contrast
        t1w         = load_nii_img_only(T1w_nii); 
        mask_t1     = load_nii_img_only(T1w_mask_nii);
        Chi_t1w     = load_nii_img_only(Chi_t1w_nii);
        img_hybrid  = compute_hybrid_t1w_chi(t1w, Chi_t1w, mask_t1, 400, 0.5 );
        %%%%%%%%%% automatic contrast matching on hybrid image, see Steps #1&2 in https://doi.org/10.1162/imag_a_00456 %%%%%%%%%%
        if isAutoMatchContrast
            img_hybrid_init_nii      = fullfile(output_tmp_dir,'Hybrid_image_init.nii.gz');
            save_nii_img_only(T1w_nii,img_hybrid_init_nii, img_hybrid );

            % quick non-linear registration, does not need to be precise at this step
            shell_script = fullfile(SEPIA_ANALYSIS_SEGMENTATION_dir,'ANTs_t1w_2_t1w_atlas_quick.sh');
            cmd = ['sh ' shell_script ' ' output_tmp_dir ' ' img_hybrid_init_nii ' ' T1w_mask_nii ' ' template_hybrid_nii];
            system(cmd);

            % apply transform to move chi and t1w image to atlas space
            gre_2_T1w_mat               = fullfile(output_tmp_dir, 'GRE_2_T1w_0GenericAffine.mat');
            t1_2_t1wTemplate_mat        = fullfile(output_tmp_dir, 'T1w_2_hybrid_0GenericAffine.mat');
            t1_2_t1wTemplate_Wrap_nii   = fullfile(output_tmp_dir, 'T1w_2_hybrid_1Warp.nii.gz');
            shell_script = fullfile(SEPIA_ANALYSIS_SEGMENTATION_dir,'ANTs_gre_2_Atlas_applyTransform.sh');
            cmd = ['sh ' shell_script ' ' output_tmp_dir ' 0 ' template_hybrid_nii ' ' Chi_nii ' ' gre_2_T1w_mat ' ' t1_2_t1wTemplate_mat ' ' t1_2_t1wTemplate_Wrap_nii];
            system(cmd);
            shell_script = fullfile(SEPIA_ANALYSIS_SEGMENTATION_dir,'ANTs_t1w_2_Atlas_applyTransform.sh');
            cmd = ['sh ' shell_script ' ' output_tmp_dir ' 0 ' template_hybrid_nii ' ' T1w_nii ' ' t1_2_t1wTemplate_mat ' ' t1_2_t1wTemplate_Wrap_nii];
            system(cmd);

            % compute contrast matched hybrid image
            Chi_t1w2atlas_nii = get_nifti_filename(Chi_nii); Chi_t1w2atlas_nii = fullfile(output_tmp_dir,[Chi_t1w2atlas_nii '_2atlas.nii.gz']);
            T1w_t1w2atlas_nii = get_nifti_filename(T1w_nii); T1w_t1w2atlas_nii = fullfile(output_tmp_dir,[T1w_t1w2atlas_nii '_2atlas.nii.gz']);
            img_hybrid = compute_MuSus100_hybrid_contrastmatch(Chi_t1w_nii,Chi_t1w2atlas_nii,T1w_nii,T1w_t1w2atlas_nii);
        end
        save_nii_img_only(T1w_nii,img_hybrid_nii, img_hybrid );
        clear img_hybrid Chi_t1w t1w mask_t1
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

        %%%%%%%%%% if speed up is enable %%%%%%%%%%
        mask_regis_nii = fullfile(output_tmp_dir,'mask_registration.nii.gz');
        if isAccelerate
            label_nii   = fullfile(MuSus100_ATLAS_HOME,'label','mixed.nii.gz');
            label       = load_nii_img_only(label_nii);
            load(fullfile(SEPIA_ANALYSIS_SEGMENTATION_dir,'private','MuSus100_ventricle_index.mat'));
            masK_ventricle = zeros(size(label)); masK_ventricle(idx) = 1;
            mask_regis = imdilate( or(label >0,masK_ventricle), strel('sphere',10));
            save_nii_img_only(label_nii,mask_regis_nii,mask_regis);
            clear mask_regis masK_ventricle label
        end
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        
        % Step 4: register hybrid image to template space, non linear transform
        shell_script = fullfile(SEPIA_ANALYSIS_SEGMENTATION_dir,'ANTs_t1w_2_t1w_atlas_template.sh');
        cmd = ['sh ' shell_script ' ' output_tmp_dir ' ' img_hybrid_nii ' ' T1w_mask_nii ' ' template_hybrid_nii ' ' num2str(isBiasFieldCorr) ' ' mask_regis_nii];
        system(cmd);

        % Step 5: Apply tranformation
        shell_script = fullfile(SEPIA_ANALYSIS_SEGMENTATION_dir,'ANTs_gre_2_t1wAtlas_applyTransform.sh');
        gre_2_T1w_mat                       = fullfile(output_tmp_dir, 'GRE_2_T1w_0GenericAffine.mat');
        t1_2_t1wTemplate_mat                = fullfile(output_tmp_dir, 'T1w_2_hybrid_0GenericAffine.mat');
        t1_2_t1wTemplate_inverseWrap_nii    = fullfile(output_tmp_dir, 'T1w_2_hybrid_1InverseWarp.nii.gz');

        mode_interp = 1; % 1=GenericLabel; 2=linear
        mode_4D     = 0; % 0=3D; 3=4D; equivalent to option -e in antsApplyTransforms see ANTs doc
        label_nii   = fullfile(MuSus100_ATLAS_HOME,'label','mixed.nii.gz');
        cmd = ['sh ' shell_script ' ' output_tmp_dir ' ' num2str(mode_interp) ' ' num2str(mode_4D) ' ' label_nii ' ' Chi_nii ' ' gre_2_T1w_mat ' ' t1_2_t1wTemplate_mat ' ' t1_2_t1wTemplate_inverseWrap_nii];
        system(cmd);

    case 2      % transformation is provided
        gre_2_T1w_mat                       = input.gre2T1wMat;
        t1_2_t1wTemplate_mat                = input.t1w2TemplateMat;
        t1_2_t1wTemplate_inverseWrap_nii    = input.t1w2TemplateiWrap;
        Chi_nii                             = input.gre;                input.chi = input.gre;
        
        % Step 1: Apply tranformation
        shell_script = fullfile(SEPIA_ANALYSIS_SEGMENTATION_dir,'ANTs_gre_2_t1wAtlas_applyTransform.sh');
        mode_interp = 1; % 1=GenericLabel; 2=linear
        mode_4D     = 0; % 0=3D; 3=4D; equivalent to option -e in antsApplyTransforms see ANTs doc
        label_nii   = fullfile(MuSus100_ATLAS_HOME,'label','mixed.nii.gz');
        cmd = ['sh ' shell_script ' ' output_tmp_dir ' ' num2str(mode_interp) ' ' num2str(mode_4D) ' ' label_nii ' ' Chi_nii ' ' gre_2_T1w_mat ' ' t1_2_t1wTemplate_mat ' ' t1_2_t1wTemplate_inverseWrap_nii];
        system(cmd);

end
 
copyfile(fullfile(output_tmp_dir, 'mixed_2gre.nii.gz'),fullfile(output_dir));
copyfile(fullfile(MuSus100_ATLAS_HOME,'label','mixed.txt'),fullfile(output_dir,'MuSus100_mixed.txt'));

%%%%%%%%%% export csv %%%%%%%%%%
label_mixed     = load_nii_img_only(fullfile(output_dir,'mixed_2gre.nii.gz'));
label_unique    = unique(label_mixed(:));
chimap          = load_nii_img_only(input.chi);
roi_mean        = nan(numel(label_unique),1);
roi_sd          = nan(numel(label_unique),1);
roi_median      = nan(numel(label_unique),1);
roi_iqr         = nan(numel(label_unique),1);
Nvoxel          = nan(numel(label_unique),1);
for kIdx = 1:numel(label_unique)
    label_curr = label_unique(kIdx);
    if label_curr > 0
        mask_label          = label_mixed == label_curr;

        Nvoxel(kIdx)        = numel(mask_label(mask_label>0));

        roi_mean(kIdx)      = mean(chimap(mask_label));
        roi_sd(kIdx)        = std(chimap(mask_label));
        roi_median(kIdx)    = median(chimap(mask_label));
        roi_iqr(kIdx)       = iqr(chimap(mask_label));
    end
end
metric = {'label index','mean (ppm)','sd (ppm)','median (ppm)','iqr (ppm)','no. voxels'};
T = table(label_unique,roi_mean,roi_sd,roi_median,roi_iqr,Nvoxel,'VariableNames',metric);
writetable(T,fullfile(output_dir,'stats_Chimap_MuSus100.csv'));
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if ~saveIntermediate; rmdir(output_tmp_dir,'s'); end

disp('Segmentation is done!')

end

function algorParam2 = check_and_set_algorithm_default(algorParam)

try algorParam2.isBiasFieldCorr     = algorParam.isBiasFieldCorr;       catch; algorParam2.isBiasFieldCorr = false; end
try algorParam2.mode                = algorParam.mode;                  catch; algorParam2.mode = 1; end
try algorParam2.isAccelerate        = algorParam.isAccelerate;          catch; algorParam2.isAccelerate = false; end
try algorParam2.saveIntermediate    = algorParam.saveIntermediate;      catch; algorParam2.saveIntermediate = false; end
try algorParam2.isAutoMatchContrast = algorParam.isAutoMatchContrast;   catch; algorParam2.isAutoMatchContrast = true; end

end

function fn = get_nifti_filename(fn)
    [~,fn,~] = fileparts(fn);
    [~,fn,~] = fileparts(fn);
end