%% Creating imrpoved hybrid image (T1w+QSM) by matching DGM contrast accurately to the hybrid atlas
% Dependencies: (1) SEPIA v1.2.2.4
% Data: (1) ABRIM dataset (https://doi.org/10.34973/7q0a-vj19) (2) MuSus-100 (https://doi.org/10.1007/s00429-022-02547-1)
% Created by Kwok-Shing Chan
% Date created: 25 April 2023
%
function hybrid_t1w_match = compute_MuSus100_hybrid_contrastmatch(chi_t1w_nii,chi_t1w2atlas_nii,t1w_nii,t1w_t1w2atlas_nii)
%% Main
% current master branch is SEPIA v1.2.2.4
% addpath('/home/common/matlab/sepia/sepia_1.2.2.4/'); % addpath('/project/3015046.06/QSM/tools/sepia');
% sepia_addpath;
% 
sepia_universal_variables;

% input filename
% chi_t1w_nii         = fullfile(proj_ABRIM_SEPIA_dir,subj_label,'anat','MuSus100_segmentation','MuSus100_intermediate_files',[subj_label '_run-1_MEGRE_QSM-LPCNN_Chimap_2T1w.nii.gz']);
% t1w_nii             = fullfile(ABRIM_MP2RAGE_dir,subj_label,'anat',[subj_label '_run-1_R1map.nii.gz']);
% chi_t1w2atlas_nii   = fullfile(proj_ABRIM_SEPIA_dir,subj_label,'anat','MuSus100_segmentation','MuSus100_intermediate_files','T1w_2_atlas_Step1',[subj_label '_run-1_MEGRE_QSM-LPCNN_Chimap_2atlas.nii.gz']);
% t1w_t1w2atlas_nii   = fullfile(proj_ABRIM_SEPIA_dir,subj_label,'anat','MuSus100_segmentation','MuSus100_intermediate_files','T1w_2_atlas_Step1',[subj_label '_run-1_R1map_2atlas.nii.gz']);
% hybrid_step2_nii    = fullfile(proj_ABRIM_SEPIA_dir,subj_label,'anat','MuSus100_segmentation','MuSus100_intermediate_files','Hybrid_image_Step2.nii.gz');
SpecifyAtlasDirectory;
%% loading data
% subject data
chi_t1w         = load_nii_img_only(chi_t1w_nii);
t1w             = load_nii_img_only(t1w_nii);
chi_t1w2atlas 	= load_nii_img_only(chi_t1w2atlas_nii);
t1w_t1w2atlas  	= load_nii_img_only(t1w_t1w2atlas_nii);

% atlas data
hybrid_atlas	= load_nii_img_only(fullfile(MuSus100_ATLAS_HOME,'atlas','hybrid.nii.gz'));
t1w_atlas       = load_nii_img_only(fullfile(MuSus100_ATLAS_HOME,'atlas','t1.nii.gz'));
chi_atlas       = load_nii_img_only(fullfile(MuSus100_ATLAS_HOME,'atlas','qsm.nii.gz'));
dbn             = uint8(load_nii_img_only(fullfile(MuSus100_ATLAS_HOME,'label','DBN.nii.gz')));
Thalamus     	= uint8(load_nii_img_only(fullfile(MuSus100_ATLAS_HOME,'label','thalamus.nii.gz')));
% mask_ventricle  = load_nii_img_only(fullfile(project_dir, 'atlas', 'MuSus100_mask_ventricle.nii.gz'))>0;
load(fullfile(SEPIA_ANALYSIS_SEGMENTATION_dir,'MuSus100_ventricle_index.mat'));
mask_ventricle = zeros(size(dbn)); mask_ventricle(idx) = 1;
mixed           = uint8(dbn + Thalamus);

% coarse tissue classification on hybrid atlas
mask_mixed              = mixed > 0;
mask_wm                 = t1w_atlas > 1800;
mask_wm(mask_mixed==1)  = 0;
mask_cgm                = and(t1w_atlas > 1000, t1w_atlas <= 1800);
mask_cgm(mask_mixed==1)	= 0;

%% obatin median value per ROI to reduce registration error
% Subcortical measures
label = unique(cat(1,dbn(:),Thalamus(:)));
label = label(label ~= 0);
t1w_atlas_label = zeros(1,length(label));
t1w_subj_label  = zeros(1,length(label));
chi_atlas_label = zeros(1,length(label));
chi_subj_label  = zeros(1,length(label));
hybrid_atlas_label = zeros(1,length(label));
label_size      = zeros(1,length(label));
counter = 0;
for k = 1:length(label)
    if label(k) ~= 0
        counter = counter + 1;
        mask = mixed == label(k);
        t1w_atlas_label(counter)    = median(t1w_atlas(mask>0));
        t1w_subj_label(counter)     = median(t1w_t1w2atlas(mask>0));
        chi_atlas_label(counter)    = median(chi_atlas(mask>0));
        chi_subj_label(counter)     = median(chi_t1w2atlas(mask>0));
        hybrid_atlas_label(counter)	= median(hybrid_atlas(mask>0));
        label_size(counter)         = numel(mask(mask>0));
    end
end

% ventricle
counter = counter+1;
t1w_atlas_label(counter)    = median(t1w_atlas(mask_ventricle>0));
t1w_subj_label(counter)     = median(t1w_t1w2atlas(mask_ventricle>0));
chi_atlas_label(counter)    = median(chi_atlas(mask_ventricle>0));
chi_subj_label(counter)     = median(chi_t1w2atlas(mask_ventricle>0));
hybrid_atlas_label(counter)	= median(hybrid_atlas(mask_ventricle>0));
label_size(counter)         = numel(mask_ventricle(mask_ventricle>0));
% global WM
counter = counter+1;
mask = mask_wm;
t1w_atlas_label(counter)    = median(t1w_atlas(mask>0));
t1w_subj_label(counter)     = median(t1w_t1w2atlas(mask>0));
chi_atlas_label(counter)    = median(chi_atlas(mask>0));
chi_subj_label(counter)     = median(chi_t1w2atlas(mask>0));
hybrid_atlas_label(counter)	= median(hybrid_atlas(mask>0));
label_size(counter)         = numel(mask(mask>0));
% global cortical GM
counter = counter+1;
mask = mask_cgm;
t1w_atlas_label(counter)    = median(t1w_atlas(mask>0));
t1w_subj_label(counter)     = median(t1w_t1w2atlas(mask>0));
chi_atlas_label(counter)    = median(chi_atlas(mask>0));
chi_subj_label(counter)     = median(chi_t1w2atlas(mask>0));
hybrid_atlas_label(counter)	= median(hybrid_atlas(mask>0));
label_size(counter)         = numel(mask(mask>0));

%% weighted by the ROI size (larger the size more relaible)
weights = label_size ./ max(label_size(1:end-3));   % exclude ventricle,WM and cGM as they have much larger size
weights(weights>1) = 1; % ventricle, WM and cGM = 1

%% Matching T1w, weighted 2th order polynomial
[xData, yData, weights_1] = prepareCurveData( t1w_subj_label, t1w_atlas_label, weights );

% Set up fittype and options (2nd order Polyfit)
ft              = fittype( 'poly2' );
opts            = fitoptions( 'Method', 'LinearLeastSquares' );
opts.Weights    = weights_1;

% Fit model to data.
[fitresult_t1w, ~] = fit( xData, yData, ft, opts );

% matching T1w
t1w_subj_match = feval(fitresult_t1w,t1w_t1w2atlas(:));
t1w_subj_match = reshape(t1w_subj_match,size(t1w_t1w2atlas));

t1w_subj_fit = feval(fitresult_t1w,xData);

%% Matching Chimap, weighted 1st order polynomial
[xData, yData, weights_1] = prepareCurveData( chi_subj_label, chi_atlas_label, weights );

% Set up fittype and options. (1st order Polyfit)
ft              = fittype( 'poly1' );
opts            = fitoptions( 'Method', 'LinearLeastSquares' );
opts.Weights    = weights_1;

% Fit model to data.
[fitresult_chi, ~] = fit( xData, yData, ft, opts );

% matching Chi
chi_subj_match = feval(fitresult_chi,chi_t1w2atlas(:));
chi_subj_match = reshape(chi_subj_match,size(chi_t1w2atlas));

chi_subj_fit    = feval(fitresult_chi,xData);

%% Matching Hybrid, weighted 1st order polynomial

[xData, yData, zData, weights_1] = prepareSurfaceData( chi_subj_fit(:), t1w_subj_fit(:), hybrid_atlas_label(:), weights(:) );

% Set up fittype and options.
ft = fittype( 'poly11' );
opts = fitoptions( 'Method', 'LinearLeastSquares' );
opts.Weights = weights_1;

% Fit model to data.
[fitresult_hybrid, ~] = fit( [xData, yData], zData, ft, opts );

%% Apply polynomial coefficients on subject T1w space
% Chimap
chi_t1w_match = feval(fitresult_chi,chi_t1w(:));
chi_t1w_match = reshape(chi_t1w_match,size(chi_t1w));
chi_t1w_match(chi_t1w_match < min(chi_atlas(:))) = min(chi_atlas(:));
chi_t1w_match(chi_t1w_match > max(chi_atlas(:))) = max(chi_atlas(:));
% T1w
t1w_t1w_match = feval(fitresult_t1w,t1w(:));
t1w_t1w_match = reshape(t1w_t1w_match,size(t1w));
t1w_t1w_match(t1w_t1w_match < min(t1w_atlas(:))) = min(t1w_atlas(:));
t1w_t1w_match(t1w_t1w_match > max(t1w_atlas(:))) = max(t1w_atlas(:));
% new Hybrid
hybrid_t1w_match = feval(fitresult_hybrid,[chi_t1w_match(:)  t1w_t1w_match(:)]);
hybrid_t1w_match = reshape(hybrid_t1w_match,size(chi_t1w));
% create weighting on (ventricle) CSF using T1w, make sure they have clear boundary
weight_t1w_csf = t1w_t1w_match / 1000;
weight_t1w_csf(weight_t1w_csf>1) = 1;
weight_t1w_csf(weight_t1w_csf<0) = 0;
% final hybrid image
hybrid_t1w_match = hybrid_t1w_match .* weight_t1w_csf;

% export new hybrid image
% save_nii_img_only(t1w_nii, hybrid_step2_nii, hybrid_t1w_match);
end