#! /bin/bash
#
# This is a shell script to register GRE image to MNI 2009c
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
SEPIA_ATLAS_dir=${SEPIA_HOME}atlas/
SEPIA_ANALYSIS_SEGMENTATION_dir=${SEPIA_HOME}analysis/segmentation/
CIT168_reinf_learn_dir=${SEPIA_ATLAS_dir}CIT168_Reinf_Learn_v1.1.0/

### user input
output_dir=$1
gre_nii=$2
gre_mask_nii=$3
t1w_nii=$4
t1w_mask_nii=$5
isBiasCorr=$6

# MNI 2009c asym T1w
mni09c_nii=${CIT168_reinf_learn_dir}MNI152-Nonlin-Asym-2009c/CIT168toMNI152-2009c_T1w_brain.nii.gz

# derived output
gre_biascorr_nii=${output_dir}GRE_biascorr.nii.gz
gre_biasField_nii=${output_dir}GRE_biasField.nii.gz
gre_brain_nii=${output_dir}GRE_brain.nii.gz
t1w_biascorr_nii=${output_dir}T1w_biascorr.nii.gz
t1w_biasField_nii=${output_dir}T1w_biasField.nii.gz
t1w_brain_nii=${output_dir}T1w_brain.nii.gz

mkdir -p $output_dir

###############################################################################
# Step 1. GRE to T1w
sh ${SEPIA_ANALYSIS_SEGMENTATION_dir}ANTs_gre_2_t1w.sh ${output_dir} ${gre_nii} ${gre_mask_nii} ${t1w_nii} ${t1w_mask_nii} ${isBiasCorr}

###############################################################################
# Step 2: T1w to MNI 2009c Asym
if [ $isBiasCorr -eq 1 ]
then
t1w_nii=${t1w_biascorr_nii}
isBiasCorr=0 #already corrected in Step 1
fi

sh ${SEPIA_ANALYSIS_SEGMENTATION_dir}ANTs_t1w_2_CIT168RLMNI2009c.sh ${output_dir} ${t1w_nii} ${t1w_mask_nii} ${isBiasCorr}
