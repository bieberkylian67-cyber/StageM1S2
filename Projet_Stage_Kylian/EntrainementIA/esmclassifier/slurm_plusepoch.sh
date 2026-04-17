#!/bin/bash
#SBATCH --job-name=entrainement_IA
#SBATCH --partition=gpu
#SBATCH --gres=gpu:1

module load python/uv
uv run python train.py \
  --sequence-length 40000 \
  --correct-dataset /gstock/user/bieber/dataset_proteine_pasdoublons.fa \
  --error-dataset /gstock/user/bieber/dataset/dataset_partial/dataset_sans_doublons.fa \
  --output-model /gstock/user/bieber/entrainement/model_partial_20epochs.pth \
  --epochs 20 \
  --esmc-model esmc_300m \
  --num-layers 1