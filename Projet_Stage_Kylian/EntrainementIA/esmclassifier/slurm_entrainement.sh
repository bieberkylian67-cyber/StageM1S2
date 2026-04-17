#!/bin/bash
#SBATCH --job-name=entrainement_IA
#SBATCH --partition=gpu
#SBATCH --gres=gpu:1

module load python/uv
uv run python train.py \
  --sequence-length 40000 \
  --correct-dataset /gstock/user/bieber/dataset_proteine.fa \
  --error-dataset /gstock/user/bieber/dataset/dataset_partial/dataset_partial_complet.fa \
  --output-model /gstock/user/bieber/entrainement/model_partial.pth \
  --epochs 5 \
  --esmc-model esmc_300m \
  --num-layers 1