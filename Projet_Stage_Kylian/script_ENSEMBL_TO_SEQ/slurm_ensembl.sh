#!/bin/bash
#SBATCH --job-name=ensembl
#SBATCH --partition=lab

conda activate stage_python39
python /home/bieber/stage/Projet_Stage_Kylian/script_ENSEMBL_TO_SEQ/Retrieve_v2.py /home/bieber/stage/Projet_Stage_Kylian/SwissProt/ENSEMBLuniq.tab /gstock/user/bieber/cds_uniq.fa /gstock/user/bieber/sequence_genomique_uniq.fa /gstock/user/bieber/exon_map_uniq.json