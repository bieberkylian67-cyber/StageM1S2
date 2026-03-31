#!/bin/bash
#SBATCH --job-name=augustus_predictions
#SBATCH --partition=lab
#SBATCH --cpus-per-task=1
#SBATCH --mem=4G
#SBATCH --array=1-10

module load annot/augustus/3.5.0
python3 /home/bieber/stage/Projet_Stage_Kylian/Augustus/script_augustus_v2.py /gstock/user/bieber/paquets_sequences/paquet_${SLURM_ARRAY_TASK_ID}.fa /gstock/user/bieber/resultat_augustus/predictions_paquet_${SLURM_ARRAY_TASK_ID}.gff3
