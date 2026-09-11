%% Creating imrpoved hybrid image (T1w+QSM) by matching DGM contrast accurately to the hybrid atlas
% Dependencies: (1) SEPIA v1.2.2.4
% Data: (1) ABRIM dataset (https://doi.org/10.34973/7q0a-vj19) (2) MuSus-100 (https://doi.org/10.1007/s00429-022-02547-1)
% Created by Kwok-Shing Chan
% Date created: 16 June 2025
%
function hybrid_t1w_match = compute_MuSus100_hybrid_given_fit(chi_t1w_nii,t1w_nii,fitresult)
%% Main
% current master branch is SEPIA v1.2.2.4
% addpath('/home/common/matlab/sepia/sepia_1.2.2.4/'); % addpath('/project/3015046.06/QSM/tools/sepia');
% sepia_addpath;
% 
sepia_universal_variables;

% input filename
SpecifyAtlasDirectory;
%% loading data
% subject data
chi_t1w         = load_nii_img_only(chi_t1w_nii);
t1w             = load_nii_img_only(t1w_nii);
t1w_atlas       = load_nii_img_only(fullfile(MuSus100_ATLAS_HOME,'atlas','t1.nii.gz'));
chi_atlas       = load_nii_img_only(fullfile(MuSus100_ATLAS_HOME,'atlas','qsm.nii.gz'));

%% Apply polynomial coefficients on subject T1w space
% Chimap
chi_t1w_match = feval(fitresult.chi,chi_t1w(:));
chi_t1w_match = reshape(chi_t1w_match,size(chi_t1w));
chi_t1w_match(chi_t1w_match < min(chi_atlas(:))) = min(chi_atlas(:));
chi_t1w_match(chi_t1w_match > max(chi_atlas(:))) = max(chi_atlas(:));
% T1w
t1w_t1w_match = feval(fitresult.t1w,t1w(:));
t1w_t1w_match = reshape(t1w_t1w_match,size(t1w));
t1w_t1w_match(t1w_t1w_match < min(t1w_atlas(:))) = min(t1w_atlas(:));
t1w_t1w_match(t1w_t1w_match > max(t1w_atlas(:))) = max(t1w_atlas(:));
% new Hybrid
hybrid_t1w_match = feval(fitresult.hybrid,[chi_t1w_match(:)  t1w_t1w_match(:)]);
hybrid_t1w_match = reshape(hybrid_t1w_match,size(chi_t1w));
% create weighting on (ventricle) CSF using T1w, make sure they have clear boundary
weight_t1w_csf = t1w_t1w_match / 1000;
weight_t1w_csf(weight_t1w_csf>1) = 1;
weight_t1w_csf(weight_t1w_csf<0) = 0;
% final hybrid image
hybrid_t1w_match = hybrid_t1w_match .* weight_t1w_csf;

end