
mode = 1;

if mode ==1
input = struct();
input(1).name   = '/project/3015069.01/pilot/qsmhubTestdata/SEPIA_QA/MGRE_original/sub-001_echo-1_part-mag_MEGRE.nii.gz';
input(2).name    = '/project/3015069.01/pilot/qsmhubTestdata/SEPIA_QA/derivatives/HD-BET/sub-001_echo-1_part-mag_MEGRE_brain_mask.nii.gz';
input(3).name   = '/project/3015069.01/pilot/qsmhubTestdata/SEPIA_QA/derivatives/FSL/fsl_anat/sub-001_syntheticMPRAGE.anat/T1_biascorr.nii.gz';
input(4).name    = '/project/3015069.01/pilot/qsmhubTestdata/SEPIA_QA/derivatives/FSL/fsl_anat/sub-001_syntheticMPRAGE.anat/T1_biascorr_brain_mask.nii.gz';
elseif mode == 2
input = struct();
% GRE to T1w, Rigid-body
input(1).name = '/project/3015069.01/pilot/qsmhubTestdata/SEPIA_QA/derivatives/SEPIA/segmentation/GRE_2_T1w_0GenericAffine.mat';
% T1w to MNI 2009c, Affine
input(2).name = '/project/3015069.01/pilot/qsmhubTestdata/SEPIA_QA/derivatives/SEPIA/segmentation/T1w_2_mni09cAsym_0GenericAffine.mat';
% T1w to MNI 2009c, Warp
input(3).name = '/project/3015069.01/pilot/qsmhubTestdata/SEPIA_QA/derivatives/SEPIA/segmentation/T1w_2_mni09cAsym_1InverseWarp.nii.gz';
input(4).name = '/project/3015069.01/pilot/qsmhubTestdata/SEPIA_QA/derivatives/FSL/fsl_anat/sub-001_syntheticMPRAGE.anat/T1_biascorr.nii.gz';
end
output_dir      = '/project/3015069.01/pilot/qsmhubTestdata/SEPIA_QA/derivatives/SEPIA/segmentation/';
algorParam.isBiasFieldCorr	= 1;
algorParam.mode             = mode;

get_CIT168_reinf_learn_labels_in_native_space(input,output_dir,algorParam);

%%

mode = 1;

if mode ==1
input = struct();
input(1).name   = '/project/3015069.01/pilot/qsmhubTestdata/SEPIA_QA/MGRE_original/sub-001_echo-1_part-mag_MEGRE.nii.gz';
input(2).name    = '/project/3015069.01/pilot/qsmhubTestdata/SEPIA_QA/derivatives/HD-BET/sub-001_echo-1_part-mag_MEGRE_brain_mask.nii.gz';
input(3).name   = '/project/3015069.01/pilot/qsmhubTestdata/SEPIA_QA/derivatives/FSL/fsl_anat/sub-001_syntheticMPRAGE.anat/T1_biascorr.nii.gz';
input(4).name    = '/project/3015069.01/pilot/qsmhubTestdata/SEPIA_QA/derivatives/FSL/fsl_anat/sub-001_syntheticMPRAGE.anat/T1_biascorr_brain_mask.nii.gz';
input(5).name    = '/project/3015069.01/pilot/qsmhubTestdata/SEPIA_QA/derivatives/SEPIA/v1p1p2/Sepia_Chimap.nii.gz';
elseif mode == 2
input = struct();
% GRE to T1w, Rigid-body
input(1).name = '/project/3015069.01/pilot/qsmhubTestdata/SEPIA_QA/derivatives/SEPIA/segmentation/GRE_2_T1w_0GenericAffine.mat';
% T1w to MNI 2009c, Affine
input(2).name = '/project/3015069.01/pilot/qsmhubTestdata/SEPIA_QA/derivatives/SEPIA/segmentation/T1w_2_mni09cAsym_0GenericAffine.mat';
% T1w to MNI 2009c, Warp
input(3).name = '/project/3015069.01/pilot/qsmhubTestdata/SEPIA_QA/derivatives/SEPIA/segmentation/T1w_2_mni09cAsym_1InverseWarp.nii.gz';
input(4).name = '/project/3015069.01/pilot/qsmhubTestdata/SEPIA_QA/derivatives/FSL/fsl_anat/sub-001_syntheticMPRAGE.anat/T1_biascorr.nii.gz';
end
output_dir      = '/project/3015069.01/pilot/qsmhubTestdata/SEPIA_QA/derivatives/SEPIA/segmentation/';
algorParam.isBiasFieldCorr	= 1;
algorParam.mode             = mode;

get_AHEAD_atlas_labels_in_native_space(input,output_dir,algorParam);

