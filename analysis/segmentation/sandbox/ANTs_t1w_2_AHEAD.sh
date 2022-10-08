#! /bin/bash
#
# This is a shell script to register T1w image to MNI 2009c
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

### user input
output_dir=$1
t1w_nii=$2
t1w_mask_nii=$3
ahead_nii=$4
isBiasCorr=$5

# derived output
t1w_brain_nii=${output_dir}Hybrid_image_brain.nii.gz

t1_2_AHEAD_vol=${output_dir}T1w_2_AHEAD_
t1_2_AHEAD_nii=${output_dir}T1w_2_AHEAD_SyN.nii.gz

mkdir -p $output_dir

###############################################################################
# Step 1. bias field correction
if [ $isBiasCorr -eq 1 ]
then
echo "Bias field correction..."

# T1w
N4BiasFieldCorrection -v -d 3 -i ${t1w_nii} \
                        -x ${t1w_mask_nii} \
                        -o [${t1w_biascorr_nii},${t1w_biasField_nii}]

export t1w_nii=${t1w_biascorr_nii}
fi
###############################################################################
# Step 2: brain extraction
MultiplyImages 3 ${t1w_nii} ${t1w_mask_nii} ${t1w_brain_nii}

################################################################################
## 3. T1w to MNI2009casym space, nonlinear

export ref_vol=${ahead_nii}
export in_vol=${t1w_brain_nii}

antsRegistration \
        --dimensionality 3 --float 0 \
        --output [${t1_2_AHEAD_vol},${t1_2_AHEAD_nii}] \
        --interpolation Linear \
        --winsorize-image-intensities [0.005,0.995] \
        --use-histogram-matching 0 \
        --initial-moving-transform [${ref_vol},${in_vol},1] \
        --transform Rigid[0.1] \
        --metric MI[${ref_vol},${in_vol},1,32,Regular,0.1] \
        --convergence [1000x500x250x100,1e-6,10] \
        --shrink-factors 4x3x2x1 \
        --smoothing-sigmas 3x2x1x0vox \
        --transform Affine[0.1] \
        --metric MI[${ref_vol},${in_vol},1,32,Regular,0.1] \
        --convergence [500x250,1e-6,10] \
        --shrink-factors 2x1 \
        --smoothing-sigmas 1x0vox \
        --transform SyN[0.1,3,0] \
        --metric CC[${ref_vol},${in_vol},1,2] \
        --convergence [500x500x250,1e-6,10] \
        --shrink-factors 4x2x1 \
        --smoothing-sigmas 2x1x0vox \
        --verbose 1 
