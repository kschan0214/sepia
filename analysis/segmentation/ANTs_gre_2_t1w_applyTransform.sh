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
gre_nii=$2
t1w_nii=$3

# derived output
gre_2_t1w_mat=${output_dir}GRE_2_T1w_0GenericAffine.mat

mkdir -p $output_dir

###############################################################################

in_nii=${gre_nii}
out_nii=$(basename -- "$in_nii")
out_nii="${out_nii%%.*}"
out_nii=${output_dir}${out_nii}_2T1w.nii.gz
antsApplyTransforms --interpolation linear --verbose 1 \
        -d 3 -e 0 -i ${gre_nii} \
        -r ${t1w_nii} \
        -t [${gre_2_t1w_mat},0] \
        -o ${out_nii}
