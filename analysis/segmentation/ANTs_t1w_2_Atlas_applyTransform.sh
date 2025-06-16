#! /bin/bash
#
# This is a shell script to register GRE image to T1w image
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
e=$2
atlas_nii=$3
t1w_nii=$4
t1_2_t1wTemplate_mat=$5
t1_2_t1wTemplate_Wrap_nii=$6

mkdir -p $output_dir

###############################################################################
# define output name
out_nii=$(basename -- "$t1w_nii")
out_nii="${out_nii%%.*}"
out_nii=${output_dir}${out_nii}_2atlas.nii.gz

# this should work for 3D but not 4D
antsApplyTransforms --interpolation linear --verbose 1 \
        -d 3 -e $e -i ${t1w_nii} \
        -r ${atlas_nii} \
        -t ${t1_2_t1wTemplate_Wrap_nii} \
        -t [${t1_2_t1wTemplate_mat}] \
        -o ${out_nii}
