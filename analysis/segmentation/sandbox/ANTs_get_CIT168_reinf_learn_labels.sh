#! /bin/bash
#
# This is a shell script to register CIT168 subcortical labels to GRE space
#
# Dependencies: (1)ANTs
#
# Creator: Kwok-shing Chan @DCCN
# kwokshing.chan@donders.ru.nl
# Date created: 6 October 2022
# Date edit:
############################################################

script_dir=`readlink -f "$0"`
SEPIA_HOME=`dirname "$script_dir"`
SEPIA_HOME=${SEPIA_HOME}/../../
SEPIA_ATLAS_dir=${SEPIA_HOME}/atlas/
SEPIA_ANALYSIS_SEGMENTATION_dir=${SEPIA_HOME}analysis/segmentation/
CIT168_reinf_learn_dir=${SEPIA_ATLAS_dir}CIT168_Reinf_Learn_v1.1.0/

# MNI 2009c asym T1w
CIT168_t1_template=CIT168toMNI152-2009c_T1w_brain
# mni09c_nii=${CIT168_reinf_learn_dir}MNI152-Nonlin-Asym-2009c/CIT168toMNI152-2009c_T1w_brain.nii.gz


### user input
mode=$1
if [ $mode -eq 1 ]  # mode 1: registration is required
then
output_dir=$2
gre_nii=$3
gre_mask_nii=$4
t1w_nii=$5
t1w_mask_nii=$6
isBiasCorr=$7

t1_2_mni2009c_inverseWrap_nii=${output_dir}T1w_2_${CIT168_t1_template}_1InverseWarp.nii.gz
t1_2_mni2009c_mat=${output_dir}T1w_2_${CIT168_t1_template}_0GenericAffine.mat
gre_2_t1w_mat=${output_dir}GRE_2_T1w_0GenericAffine.mat

elif [ $mode -eq 2 ] # mode 2: transformation is available
then
output_dir=$2
gre_2_t1w_mat=$3
t1_2_mni2009c_mat=$4
t1_2_mni2009c_inverseWrap_nii=$5
t1w_nii=$5
fi

# derived output
gre_biascorr_nii=${output_dir}GRE_biascorr.nii.gz
gre_biasField_nii=${output_dir}GRE_biasField.nii.gz
gre_brain_nii=${output_dir}GRE_brain.nii.gz
t1w_biascorr_nii=${output_dir}T1w_biascorr.nii.gz
t1w_biasField_nii=${output_dir}T1w_biasField.nii.gz
t1w_brain_nii=${output_dir}T1w_brain.nii.gz

mkdir -p $output_dir

###############################################################################
# Display input
echo "Output directory: ${output_dir}"
echo "GRE image: ${gre_nii}"
echo "GRE mask: ${gre_mask_nii}"
echo "T1 image: ${t1w_nii}"
echo "T1 mask: ${t1w_mask_nii}"
echo "Bias field correction (0:false; 1:true): ${isBiasCorr}"

# ###############################################################################
# # Step 1. GRE to MNI 2009c
# if [ $mode -eq 1 ]  # mode 1: registration is required
# then
# sh ${SEPIA_ANALYSIS_SEGMENTATION_dir}ANTs_gre_2_CIT168RLMNI2009c.sh ${output_dir} ${gre_nii} ${gre_mask_nii} ${t1w_nii} ${t1w_mask_nii} ${isBiasCorr}
# fi

###############################################################################
# Step 2: Apply transforms to label
# 2.1: probabilistic labels
in_dir=${CIT168_reinf_learn_dir}MNI152-Nonlin-Asym-2009c/

in_vol=CIT168toMNI152-2009c_prob
in_nii=${in_dir}${in_vol}.nii.gz
out_nii=${output_dir}${in_vol}_2gre.nii.gz
antsApplyTransforms --interpolation linear --verbose 1 \
        -d 3 -e 3 -i ${in_nii} \
        -r ${gre_nii} \
        -t [${gre_2_t1w_mat},1] \
        -t [${t1_2_mni2009c_mat},1] \
        -t ${t1_2_mni2009c_inverseWrap_nii} \
        -o ${out_nii}

# 2.2: deterministic labels
in_vol=CIT168toMNI152-2009c_det
in_nii=${in_dir}${in_vol}.nii.gz
out_nii=${output_dir}${in_vol}_2gre.nii.gz
antsApplyTransforms --interpolation GenericLabel --verbose 1 \
        -d 3 -e 0 -i ${in_nii} \
        -r ${gre_nii} \
        -t [${gre_2_t1w_mat},1] \
        -t [${t1_2_mni2009c_mat},1] \
        -t ${t1_2_mni2009c_inverseWrap_nii} \
        -o ${out_nii}

cp ${CIT168_reinf_learn_dir}labels.txt ${output_dir}labels.txt
