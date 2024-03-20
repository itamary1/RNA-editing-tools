#!/bin/bash
#
SCRIPT_DIR="$(dirname "$(readlink -f "$0")")"
DOWNLOADED_DIR=${3:-"Downloaded"};
STAR_INDEX_DIR=${4:-"STAR2.7.10b"};
SALMON_INDEX_DIR=${5:-"Salmon1.10.2"};
MAIN_DIR=${6:-"Heavy_Resources"};


#---------------------------------------------------------------------------
# DOWNLOAD RESOURCES
#---------------------------------------------------------------------------
cd $SCRIPT_DIR
# download resources
sh ${SCRIPT_DIR}/resources/download_resources.sh $MAIN_DIR/$DOWNLOADED_DIR
# index genomes for STAR
sh ${SCRIPT_DIR}/resources/STAR-genome-generate.sh $MAIN_DIR/$DOWNLOADED_DIR $MAIN_DIR/$STAR_INDEX_DIR
# index transcriptome for Salmon
sh ${SCRIPT_DIR}/resources/Salmon-genome-generate.sh $MAIN_DIR/$DOWNLOADED_DIR $MAIN_DIR/$SALMON_INDEX_DIR

# -------------------------------------------------------------------------------
# UPDATE CONFIGS' PATHS
# -------------------------------------------------------------------------------

# update full path for the basic config
sed -i "s|REPLACE_PATH|${SCRIPT_DIR}|g" ${SCRIPT_DIR}/Configs/Docker/nextflow.base.dockers.config

# salmon
# update full path for the basic config
sed -i "s|ReplaceSalmonPath|${SCRIPT_DIR}/${MAIN_DIR}/${SALMON_INDEX_DIR}|g" ${SCRIPT_DIR}/Configs/Docker/SubP_configs/salmon.nf.docker.config
sed -i "s|ReplaceGenemapPath|${SCRIPT_DIR}/${MAIN_DIR}/${DOWNLOADED_DIR}/RefSeqAnnotations|g" ${SCRIPT_DIR}/Configs/Docker/SubP_configs/salmon.nf.docker.config
# star
sed -i "s|ReplaceSTARGenomePATH|${SCRIPT_DIR}/${MAIN_DIR}/$STAR_INDEX_DIR|g" ${SCRIPT_DIR}/Configs/Docker/SubP_configs/star.nf.docker.config
# EI
sed -i "s|RelpaceGenomePath|${SCRIPT_DIR}/${MAIN_DIR}/${DOWNLOADED_DIR}/Genomes|g" ${SCRIPT_DIR}/Configs/Docker/SubP_configs/editing_index.nf.docker.config
sed -i "s|ReplaceRefseqFPath|${SCRIPT_DIR}/${MAIN_DIR}/${DOWNLOADED_DIR}/RefSeqAnnotations|g" ${SCRIPT_DIR}/Configs/Docker/SubP_configs/editing_index.nf.docker.config
sed -i "s|ReplaceExpressionPath|${SCRIPT_DIR}/${MAIN_DIR}/${DOWNLOADED_DIR}/GenesExpression|g" ${SCRIPT_DIR}/Configs/Docker/SubP_configs/editing_index.nf.docker.config
sed -i "s|ReplaceSNPPATH|${SCRIPT_DIR}/${MAIN_DIR}/${DOWNLOADED_DIR}/SNPs|g" ${SCRIPT_DIR}/Configs/Docker/SubP_configs/editing_index.nf.docker.config
sed -i "s|ReplaceRegionsPath|${SCRIPT_DIR}/${MAIN_DIR}/${DOWNLOADED_DIR}/Regions|g" ${SCRIPT_DIR}/Configs/Docker/SubP_configs/editing_index.nf.docker.config
# CMpileup
sed -i "s|RelpaceGenomePath|${SCRIPT_DIR}/${MAIN_DIR}/${DOWNLOADED_DIR}/Genomes|g" ${SCRIPT_DIR}/Configs/Docker/SubP_configs/cmpileup.nf.docker.config