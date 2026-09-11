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
% Date modified: 16 June 2025
%
%
function get_AHEAD_atlas_labels(input,output_dir,algorParam)

sepia_universal_variables;

algorParam          = check_and_set_algorithm_default(algorParam);
isBiasFieldCorr     = algorParam.isBiasFieldCorr;
mode                = algorParam.mode;
isDownsample        = algorParam.isDownsample;
isAccelerate        = algorParam.isAccelerate;
saveIntermediate    = algorParam.saveIntermediate;
targetResolution    = algorParam.targetResolution;
isAutoMatchContrast = algorParam.isAutoMatchContrast;

% setup ANTs in environment
SpecifyToolboxesDirectory;
% test if antsRegistration exists
[status,~] = system('antsRegistration');if status == 127; setenv('PATH', [getenv('PATH') ':' ANTS_HOME]); end
% get atlas directory
SpecifyAtlasDirectory;

output_tmp_dir = fullfile(output_dir,'AHEAD_intermediate_files',filesep);
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

%%%%%%%%%% Atlas data %%%%%%%%%%
t1w_atlas_nii       = fullfile(AHEAD_ATLAS_HOME,'Templates_mni09b','ahead_final_med_r1map_n104_mni09b.nii.gz');
chi_atlas_nii       = fullfile(AHEAD_ATLAS_HOME,'Templates_mni09b','ahead_final_med_qsm_n104_mni09b.nii.gz');
mask_ventricle_mat  = fullfile(SEPIA_ANALYSIS_SEGMENTATION_dir,'private','AHEAD_ventricle_index.mat');
label_dir           = fullfile(fullfile(AHEAD_ATLAS_HOME,'structures_mni09b'));
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%% Check downsample %%%%%%%%%%
voxelSize = [0.5,0.5,0.5];
if isDownsample
    disp('Down-sampling AHEAD atlas using cubic interpolation...')

    scaleFactor = voxelSize(1)/targetResolution; % allow this to change
    voxelSize   = voxelSize*(1/scaleFactor);

    t1w_atlas_ds_nii = fullfile(output_tmp_dir,strcat('ahead_final_med_r1map_n104_mni09b.nii.gz'));
    chi_atlas_ds_nii = fullfile(output_tmp_dir,strcat('ahead_final_med_qsm_n104_mni09b.nii.gz'));
    label_ds_dir     = fullfile(fullfile(output_tmp_dir,'structures_mni09b'));
    mkdir(label_ds_dir)

    fprintf('Target resolution [x,y,z]: %.2fmm x %.2fmm x %.2fmm\n',voxelSize(1),voxelSize(2),voxelSize(3));
    
    % R1
    tmp = load_nii_img_only(t1w_atlas_nii); 
    tmp = imresize3(tmp, scaleFactor, 'cubic'); 
    save_nii_img_only_newres(t1w_atlas_nii,t1w_atlas_ds_nii,tmp,16,voxelSize);

    % Chimap
    tmp = load_nii_img_only(chi_atlas_nii); 
    tmp = imresize3(tmp, scaleFactor, 'cubic'); 
    save_nii_img_only_newres(chi_atlas_nii,chi_atlas_ds_nii,tmp,16,voxelSize);

    % mask
    mask_list = dir(fullfile(label_dir,'*_mask-*nii*'));
    for k = 1:numel(mask_list)
        label_nii       = fullfile(label_dir, mask_list(k).name);
        label_ds_nii    = fullfile(label_ds_dir,mask_list(k).name);

        tmp             = imresize3( load_nii_img_only(label_nii), scaleFactor, "nearest"); 
        save_nii_img_only_newres(label_nii,label_ds_nii,tmp,16,voxelSize);
    end

    % probablistic mask
    prob_list = dir(fullfile(label_dir,'*_proba-*nii*'));
    for k = 1:numel(prob_list)
        label_nii       = fullfile(label_dir,prob_list(k).name);
        label_ds_nii    = fullfile(label_ds_dir,prob_list(k).name);

        tmp             = imresize3( load_nii_img_only(label_nii), scaleFactor, "cubic");  tmp(tmp<0) = 0; tmp(tmp>1) = 1; % make sure no negative probability after downsampling
        save_nii_img_only_newres(label_nii,label_ds_nii,tmp,16,voxelSize);
    end

    load(mask_ventricle_mat);
    mask_ventricle          = zeros(size(load_nii_img_only(t1w_atlas_nii))); mask_ventricle(idx) = 1; mask_ventricle = imresize3( mask_ventricle, scaleFactor, "nearest"); 
    idx                     = find(mask_ventricle>0);
    mask_ventricle_ds_mat   = fullfile(output_tmp_dir,'private''AHEAD_ventricle_index.mat'); save(mask_ventricle_ds_mat,'idx')

    % update atlas details
    t1w_atlas_nii       = t1w_atlas_ds_nii;
    chi_atlas_nii       = chi_atlas_ds_nii;
    label_dir           = label_ds_dir;
    mask_ventricle_mat  = mask_ventricle_ds_mat;

    fprintf('Done!\n')

    clear tmp
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%



switch mode
    case 1      % registration is required
        GRE_nii      = input.gre;
        GRE_mask_nii = input.greMask;
        T1w_nii      = input.t1w;
        T1w_mask_nii = input.t1wMask;
        Chi_nii      = input.chi;
        
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
        
        % Step 3.1: reate hybrid T1w+QSM image on the atlas data
        % t1w         = load_nii_img_only(t1w_atlas_nii);
        % Chi_t1w     = load_nii_img_only(chi_atlas_nii);
        % mask_t1     = t1w ~= 0;
        % img_hybrid  = compute_hybrid_t1w_chi(t1w, Chi_t1w, mask_t1, 400, 0.5 );
        % create a hybrid image that matched MuSus-100 atlas hybrid, which has better contrast in basal ganglia
        load('private/AHEAD_2_MuSus100_hybrid.mat');
        template_hybrid_nii = fullfile(output_tmp_dir,'Template_hybrid_image.nii.gz');
        img_hybrid          = compute_MuSus100_hybrid_given_fit(chi_atlas_nii,t1w_atlas_nii,fitres);
        save_nii_img_only(t1w_atlas_nii,template_hybrid_nii, img_hybrid);
        clear img_hybrid Chi_t1w t1w mask_t1

        % Step 3.2: create hybrid T1w+QSM image on the provided data
        Chi_t1w_nii         = fullfile(output_tmp_dir,[Chi_t1w_nii '_2T1w.nii.gz']);
        img_hybrid_nii      = fullfile(output_tmp_dir,'Hybrid_image.nii.gz');
        
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
            t1_2_t1wTemplate_mat        = fullfile(output_tmp_dir, 'T1w_2_Template_hybrid_image_0GenericAffine.mat');
            t1_2_t1wTemplate_Wrap_nii   = fullfile(output_tmp_dir, 'T1w_2_Template_hybrid_image_1Warp.nii.gz');
            shell_script = fullfile(SEPIA_ANALYSIS_SEGMENTATION_dir,'ANTs_gre_2_Atlas_applyTransform.sh');
            cmd = ['sh ' shell_script ' ' output_tmp_dir ' 0 ' template_hybrid_nii ' ' Chi_nii ' ' gre_2_T1w_mat ' ' t1_2_t1wTemplate_mat ' ' t1_2_t1wTemplate_Wrap_nii];
            system(cmd);
            shell_script = fullfile(SEPIA_ANALYSIS_SEGMENTATION_dir,'ANTs_t1w_2_Atlas_applyTransform.sh');
            cmd = ['sh ' shell_script ' ' output_tmp_dir ' 0 ' template_hybrid_nii ' ' T1w_nii ' ' t1_2_t1wTemplate_mat ' ' t1_2_t1wTemplate_Wrap_nii];
            system(cmd);

            % compute contrast matched hybrid image
            AHEAD_input_struct.hybrid_nii           = template_hybrid_nii;
            AHEAD_input_struct.t1w_nii              = t1w_atlas_nii;
            AHEAD_input_struct.chi_nii              = chi_atlas_nii;
            AHEAD_input_struct.label_dir            = label_dir;
            AHEAD_input_struct.mask_ventricle_mat   = mask_ventricle_mat;

            Chi_t1w2atlas_nii = get_nifti_filename(Chi_nii); Chi_t1w2atlas_nii = fullfile(output_tmp_dir,[Chi_t1w2atlas_nii '_2atlas.nii.gz']);
            T1w_t1w2atlas_nii = get_nifti_filename(T1w_nii); T1w_t1w2atlas_nii = fullfile(output_tmp_dir,[T1w_t1w2atlas_nii '_2atlas.nii.gz']);
            img_hybrid = compute_AHEAD_hybrid_contrastmatch(AHEAD_input_struct,Chi_t1w_nii,Chi_t1w2atlas_nii,T1w_nii,T1w_t1w2atlas_nii);
        end
        save_nii_img_only(T1w_nii,img_hybrid_nii, img_hybrid );
        clear img_hybrid Chi_t1w t1w mask_t1
        
        %%%%%%%%%% Speed up using label mask (and ventricle) only %%%%%%%%%%
        mask_regis_nii = fullfile(output_tmp_dir,'mask_registration.nii.gz');
        if isAccelerate
            mask_list = dir(fullfile(label_dir,'*_mask-*nii*'));
            for k = 1:numel(mask_list)
                label_nii       = fullfile(label_dir, mask_list(k).name);
                if k == 1; label = load_nii_img_only(label_nii); else; label = label + load_nii_img_only(label_nii); end
            end

            load(mask_ventricle_mat);
            mask_ventricle = zeros(size(label)); mask_ventricle(idx) = 1;

            mask_regis = imdilate( or(label >0,mask_ventricle), strel('sphere',ceil(10/voxelSize(1))));   % 10mm radius
            save_nii_img_only(label_nii,mask_regis_nii,mask_regis);
            clear mask_regis mask_ventricle label
        end
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        
        % Step 4: register hybrid image to template space, non linear transform
        shell_script = fullfile(SEPIA_ANALYSIS_SEGMENTATION_dir,'ANTs_t1w_2_t1w_atlas_template.sh');
        cmd = ['sh ' shell_script ' ' output_tmp_dir ' ' img_hybrid_nii ' ' T1w_mask_nii ' ' template_hybrid_nii ' ' num2str(isBiasFieldCorr) ' ' mask_regis_nii];
        system(cmd);

        % Step 5: Apply tranformation
        shell_script = fullfile(SEPIA_ANALYSIS_SEGMENTATION_dir,'ANTs_gre_2_t1wAtlas_applyTransform.sh');
        gre_2_T1w_mat                       = fullfile(output_tmp_dir, 'GRE_2_T1w_0GenericAffine.mat');
        t1_2_t1wTemplate_mat                = fullfile(output_tmp_dir, 'T1w_2_Template_hybrid_image_0GenericAffine.mat');
        t1_2_t1wTemplate_inverseWrap_nii    = fullfile(output_tmp_dir, 'T1w_2_Template_hybrid_image_1InverseWarp.nii.gz');
        
        mask_list = dir(fullfile(label_dir,'*_mask-*nii*'));
        mode_interp = 1; % 1=GenericLabel; 2=linear
        mode_4D     = 0; % 0=3D; 3=4D; equivalent to option -e in antsApplyTransforms see ANTs doc
        for k = 1:numel(mask_list)
            label_nii   = fullfile(AHEAD_ATLAS_HOME,'structures_mni09b',mask_list(k).name);
            cmd = ['sh ' shell_script ' ' output_tmp_dir ' ' num2str(mode_interp) ' ' num2str(mode_4D) ' ' label_nii ' ' Chi_nii ' ' gre_2_T1w_mat ' ' t1_2_t1wTemplate_mat ' ' t1_2_t1wTemplate_inverseWrap_nii];
            system(cmd);
        end

        prob_list = dir(fullfile(label_dir,'*_proba-*nii*'));
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

%%%%%%%%%% export csv %%%%%%%%%%
mask_list = dir(fullfile(output_dir,'*_mask-*nii*'));
% label_unique    = unique(label_mixed(:));
chimap          = load_nii_img_only(input.chi); nii = load_untouch_nii(chimap);
roi_mean        = nan(numel(mask_list),1);
roi_sd          = nan(numel(mask_list),1);
roi_median      = nan(numel(mask_list),1);
roi_iqr         = nan(numel(mask_list),1);
Nvoxel          = nan(numel(mask_list),1);
roi_wmean       = nan(numel(mask_list),1);
roi_wsd         = nan(numel(mask_list),1);
roi_wvol        = nan(numel(mask_list),1);
label_name      = {};
for kIdx = 1:numel(mask_list)

        label_name(kIdx)    = extractBetween(mask_list(kIdx).name,'atlas-','_mask-');
    
        label_nii           = fullfile(output_dir, mask_list(kIdx).name);
        proba_nii           = fullfile(output_dir, replace(mask_list(kIdx).name,'mask','proba'));

        mask_label          = load_nii_img_only(label_nii) > 0;
        prob_label          = load_nii_img_only(proba_nii);

        % mask stats
        Nvoxel(kIdx)        = numel(mask_label(mask_label>0));
        roi_mean(kIdx)      = mean(chimap(mask_label));
        roi_sd(kIdx)        = std(chimap(mask_label));
        roi_median(kIdx)    = median(chimap(mask_label));
        roi_iqr(kIdx)       = iqr(chimap(mask_label));
        
        % probability mask stats
        roi_wmean(kIdx)     = sum(chimap(:).*prob_label(:))/sum(prob_label(:));
        roi_wsd(kIdx)       = sqrt(sum(prob_label(:) .* (chimap(:) - roi_wmean(kIdx)).^2) / sum(prob_label(:)));
        roi_wvol(kIdx)      = sum(prob_label(:))*prod(nii.hdr.dime.pixdim(2:4));

end
metric = {'label name','mean (ppm)','sd (ppm)','median (ppm)','iqr (ppm)','no. voxels','weighted mean (ppm)','weighted sd (ppm)','weighted volume (mm3)'};
T = table(label_name(:),roi_mean,roi_sd,roi_median,roi_iqr,Nvoxel,roi_wmean,roi_wsd,roi_wvol,'VariableNames',metric);
writetable(T,fullfile(output_dir,'stats_Chimap_AHEAD.csv'));
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if ~saveIntermediate; rmdir(output_tmp_dir,'s'); end

disp('Segmentation is done!')

end

function algorParam2 = check_and_set_algorithm_default(algorParam)

try algorParam2.isBiasFieldCorr     = algorParam.isBiasFieldCorr;   catch; algorParam2.isBiasFieldCorr = false; end
try algorParam2.mode                = algorParam.mode;              catch; algorParam2.mode = 1; end
try algorParam2.isDownsample        = algorParam.isDownsample;      catch; algorParam2.isDownsample = false; end
try algorParam2.isAccelerate        = algorParam.isAccelerate;      catch; algorParam2.isAccelerate = false; end
try algorParam2.saveIntermediate    = algorParam.saveIntermediate;  catch; algorParam2.saveIntermediate = false; end
try algorParam2.targetResolution    = algorParam.targetResolution;  catch; algorParam2.targetResolution = 0.5; end
try algorParam2.isAutoMatchContrast = algorParam.isAutoMatchContrast;   catch; algorParam2.isAutoMatchContrast = true; end


end

function fn = get_nifti_filename(fn)
    [~,fn,~] = fileparts(fn);
    [~,fn,~] = fileparts(fn);
end