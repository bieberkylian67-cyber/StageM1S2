#!/bin/bash
#SBATCH --job-name=entrainement_IA
#SBATCH --partition=gpu
#SBATCH --gres=gpu:1

module load python/uv
uv run python train.py \
  --sequence-length 40000 \
  --correct-dataset /gstock/user/bieber/dataset_proteine_pasdoublons.fa \
  --error-dataset /gstock/user/bieber/dataset/dataset_partial/dataset_sans_doublons.fa \
  --output-model /gstock/user/bieber/entrainement/model_partial_600m.pth \
  --epochs 5 \
  --esmc-model esmc_600m \
  --num-layers 1