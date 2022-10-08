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
mode=$2
e=$3
label_nii=$4
gre_nii=$5
gre_2_t1w_mat=$6
t1_2_t1wTemplate_mat=$7
t1_2_t1wTemplate_inverseWrap_nii=$8

mkdir -p $output_dir

###############################################################################
# define output name
out_nii=$(basename -- "$label_nii")
out_nii="${out_nii%%.*}"
out_nii=${output_dir}${out_nii}_2gre.nii.gz

# this should work for 3D but not 4D
if [ $mode -eq 1 ]
then
antsApplyTransforms --interpolation GenericLabel --verbose 1 \
        -d 3 -e $e -i ${label_nii} \
        -r ${gre_nii} \
        -t [${gre_2_t1w_mat},1] \
        -t [${t1_2_t1wTemplate_mat},1] \
        -t ${t1_2_t1wTemplate_inverseWrap_nii} \
        -o ${out_nii}
elif [ $mode -eq 2 ]
then
antsApplyTransforms --interpolation linear --verbose 1 \
        -d 3 -e $e -i ${label_nii} \
        -r ${gre_nii} \
        -t [${gre_2_t1w_mat},1] \
        -t [${t1_2_t1wTemplate_mat},1] \
        -t ${t1_2_t1wTemplate_inverseWrap_nii} \
        -o ${out_nii}
fi

