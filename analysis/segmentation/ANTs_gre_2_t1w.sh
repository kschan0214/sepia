#! /bin/bash
#
# This is a shell script to register GRE image to T1w image
#
# Dependencies: (1)ANTs
#
# Creator: Kwok-shing Chan @DCCN
# kwokshing.chan@donders.ru.nl
# Date created: 6 October 2022
# Date edit: 15 June 2025
############################################################

# export ITK_GLOBAL_DEFAULT_NUMBER_OF_THREADS=8

script_dir=`readlink -f "$0"`
SEPIA_HOME=`dirname "$script_dir"`
SEPIA_HOME=${SEPIA_HOME}/../../
SEPIA_ATLAS_dir=${SEPIA_HOME}/atlas/

### user input
output_dir=$1
gre_nii=$2
gre_mask_nii=$3
t1w_nii=$4
t1w_mask_nii=$5
isBiasCorr=$6

# derived output
gre_biascorr_nii=${output_dir}GRE_biascorr.nii.gz
gre_biasField_nii=${output_dir}GRE_biasField.nii.gz
gre_brain_nii=${output_dir}GRE_brain.nii.gz
t1w_biascorr_nii=${output_dir}T1w_biascorr.nii.gz
t1w_biasField_nii=${output_dir}T1w_biasField.nii.gz
t1w_brain_nii=${output_dir}T1w_brain.nii.gz

gre_2_t1w_vol=${output_dir}GRE_2_T1w_
gre_2_t1w_nii=${output_dir}GRE_2_T1w_dof6.nii.gz

mkdir -p $output_dir

###############################################################################
# Step 1. bias field correction
if [ $isBiasCorr -eq 1 ]
then
echo "Bias field correction..."

# GRE
N4BiasFieldCorrection -v -d 3 -i ${gre_nii} \
                        -x ${gre_mask_nii} \
                        -o [${gre_biascorr_nii},${gre_biasField_nii}]

# T1w
N4BiasFieldCorrection -v -d 3 -i ${t1w_nii} \
                        -x ${t1w_mask_nii} \
                        -o [${t1w_biascorr_nii},${t1w_biasField_nii}]

export gre_nii=${gre_biascorr_nii}
export t1w_nii=${t1w_biascorr_nii}
fi

###############################################################################
# Step 2: brain extraction
echo "Applying brain masks..."

MultiplyImages 3 ${gre_nii} ${gre_mask_nii} ${gre_brain_nii}
MultiplyImages 3 ${t1w_nii} ${t1w_mask_nii} ${t1w_brain_nii}

################################################################################
## 3. GRE to T1w space, rigid body transform
echo "Computing GRE to T1w transformation..."

ref_vol=${t1w_brain_nii}
in_vol=${gre_brain_nii}
antsRegistration \
        --dimensionality 3 --float 0 \
        --output [${gre_2_t1w_vol},${gre_2_t1w_nii}] \
        --interpolation Linear \
        --winsorize-image-intensities [0.005,0.995] \
        --use-histogram-matching 0 \
        --initial-moving-transform [${ref_vol},${in_vol},1] \
        --transform Rigid[0.1] \
        --metric MI[${ref_vol},${in_vol},1,32,Regular,0.25] \
        --convergence [1000x500x250x100,1e-6,10] \
        --shrink-factors 4x3x2x1 \
        --smoothing-sigmas 3x2x1x0vox \
        --verbose 1 
