#!/bin/bash

#SBATCH --job-name=augustus_comparaison
#SBATCH --partition=lab
#SBATCH --cpus-per-task=1
#SBATCH --mem=2G
#SBATCH --array=1-10

ID=$SLURM_ARRAY_TASK_ID

FICHIER_INPUT_AUGUSTUS="/gstock/user/bieber/resultat_augustus/predictions_paquet_${ID}.gff3"

FICHIER_EXON_MAP_ENSEMBL="/gstock/user/bieber/exonmap_final.json"

FICHIER_OUTPUT_DATASET="/gstock/user/bieber/dataset/dataset_partie_${ID}"

python3 /home/bieber/stage/Projet_Stage_Kylian/Comparateur/comparateur_v2.py $FICHIER_INPUT_AUGUSTUS $FICHIER_EXON_MAP_ENSEMBL $FICHIER_OUTPUT_DATASET