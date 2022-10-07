#! /bin/bash
#
# This is a shell script to download atlas from online source
#
# Creator: Kwok-shing Chan @DCCN
# kwokshing.chan@donders.ru.nl
# Date created: 6 October 2022
# Date last edit:
############################################################

script_dir=`readlink -f "$0"`
SEPIA_HOME=`dirname "$script_dir"`
# echo $SEPIA_HOME

SEPIA_ATLAS_dir=${SEPIA_HOME}/atlas/

############ CIT168 reinf learning atlas ############
CIT168_reinf_learn_dir=${SEPIA_ATLAS_dir}CIT168_Reinf_Learn_v1.1.0/
CIT168_reinf_learn_file='CIT168_Reinf_Learn_v1.1.0.zip'
CIT168_reinf_learn_link='https://files.osf.io/v1/resources/jkzwp/providers/osfstorage/5b11f8d6f1f288000d6343aa/?zip='

# create new directory for atlas
mkdir -p ${CIT168_reinf_learn_dir}

# download atlas from online source
wget -O ${SEPIA_ATLAS_dir}${CIT168_reinf_learn_file} --no-check-certificate ${CIT168_reinf_learn_link}
# unzip the file
unzip ${SEPIA_ATLAS_dir}${CIT168_reinf_learn_file} -d ${CIT168_reinf_learn_dir}
# delete the zip file
rm ${SEPIA_ATLAS_dir}${CIT168_reinf_learn_file}
############################################################

######################## AHEAD atlas ########################
AHEAD_atlas_dir=${SEPIA_ATLAS_dir}AHEAD_atlas/
AHEAD_atlas_file1='structure_mni09b.tar.gz'
AHEAD_atlas_file2='Templates_mni09b.tar.gz'
AHEAD_atlas_link1='https://uvaauas.figshare.com/ndownloader/files/21209229'
AHEAD_atlas_link2='https://uvaauas.figshare.com/ndownloader/files/21209235'

# create new directory for atlas
mkdir -p ${AHEAD_atlas_dir}

# download atlas from online source
wget -O ${SEPIA_ATLAS_dir}${AHEAD_atlas_file1} --no-check-certificate ${AHEAD_atlas_link1}
wget -O ${SEPIA_ATLAS_dir}${AHEAD_atlas_file2} --no-check-certificate ${AHEAD_atlas_link2}
# untar the file
tar -xvzf ${SEPIA_ATLAS_dir}${AHEAD_atlas_file1} -C ${AHEAD_atlas_dir}
tar -xvzf ${SEPIA_ATLAS_dir}${AHEAD_atlas_file2} -C ${AHEAD_atlas_dir}
# delete the tar balls
rm ${SEPIA_ATLAS_dir}${AHEAD_atlas_file1}
rm ${SEPIA_ATLAS_dir}${AHEAD_atlas_file2}

############################################################
