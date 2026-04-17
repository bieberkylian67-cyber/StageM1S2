#!/bin/bash
#SBATCH --job-name=inference
#SBATCH --partition=gpu
#SBATCH --gres=gpu:1

module load python/uv
uv run python infer_matrice_confusion.py \
  --correct-fasta /home/bieber/stage/Projet_Stage_Kylian/EntrainementIA/esmclassifier/g3po_bbs_correct.fasta \
  --error-fasta /home/bieber/stage/Projet_Stage_Kylian/EntrainementIA/esmclassifier/g3po_bbs_error.fasta \
  --model-path /gstock/user/bieber/entrainement/model_partial_epoch_esm600m.pth \
  --esmc-model esmc_600m \
  --num-layers 2 \
  --max-length 2000
  