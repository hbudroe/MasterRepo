#!/bin/sh
#SBATCH --job-name="convert-to-qza"
#SBATCH --partition=batch
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem-per-cpu=12G
#SBATCH --time=3-01:00:00
#SBATCH --mail-user=hmb25721@uga.edu
#SBATCH --mail-type=BEGIN,END,FAIL
#SBATCH -e 01-qiime-import.err-%N-%j
#SBATCH -o 01-qiime-import.out-%N-%j

#This script is to convert a csv file of the meioDB taxonomic strings into a .qza qiime artifact file so it can be used for taxonomic assignment in the qiime2 pipeline
module load QIIME2/2025.4-amplicon

INPUT=/home/hmb25721/MasterRepo/MeioDB/meiodb_v5_updatedtaxstringsmeiofauna_09222025_intermediatev6.tsv
OUTPUT=/home/hmb25721/MasterRepo/MeioDB/meiodb_v5_updatedtaxstringsmeiofauna_09222025_intermediatev6.qza

qiime tools import \
  --type FeatureData[Taxonomy] \
  --input-path ${INPUT} \
  --output-path ${OUTPUT} \
  --input-format TSVTaxonomyFormat 

#types ? 'FeatureData[SILVATaxonomy]' tried & error  No transformation from <class 'q2_types.feature_data._formats.TSVTaxonomyFormat'> to <class 'rescript.types._format.SILVATaxonomyDirectoryFormat'>
