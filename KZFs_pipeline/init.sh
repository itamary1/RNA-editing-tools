#!/bin/bash


# ##################################################################################################
#
#                             INIT KZFs PIPELINE
#
# ##################################################################################################
#  
#    run this script  once before the first time you run KZFs


SCRIPT_DIR="$(dirname "$(readlink -f "$0")")"
PARENT_DIR="$(dirname ${SCRIPT_DIR})"
EDITING_PIPELINE_DIR="${PARENT_DIR}/RNA_Editing_Pipeline"
${SCRIPT_DIR}/${MAIN_DIR}/${DOWNLOADED_DIR}/SNPs

# set resources dirs names
DOWNLOADED_DIR="Downloaded"
MAIN_DIR="Heavy_Resources"


# install argparse if dont exist
if [ ! -d "~/argparse-bash" ]; then 
 git clone https://github.com/nhoffman/argparse-bash.git $HOME/argparse-bash
fi

# run the the RNA editing pipline init.sh if it wasn't run before
if [ ! -d "${EDITING_PIPELINE_DIR}/Heavy_Resources" ]; then 
 echo "runing RNA_Editing_Pipeline/init.sh it will take few hours"
 bash "${EDITING_PIPELINE_DIR}/init.sh"
fi

# fix project path in all folders
find ${SCRIPT_DIR} \( -type d -name .git -prune \) -o -type f -print0 | xargs -0 sed -i "s|Replace_parentPath|${PARENT_DIR}|g"

# fix SNPs path 
find ${SCRIPT_DIR} \( -type d -name .git -prune \) -o -type f -print0 | xargs -0 sed -i "s|ReplaceSNPPATH|${EDITING_PIPELINE_DIR}/${MAIN_DIR}/${DOWNLOADED_DIR}/SNPs|g"
