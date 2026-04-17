#!/bin/bash
#SBATCH --job-name=entrainV2
#SBATCH --partition=gpu
#SBATCH --gres=gpu:1

module load python/uv
uv run python train.py \
  --sequence-length 40000 \
  --correct-dataset /gstock/user/bieber/dataset_proteine.fa \
  --error-dataset /gstock/user/bieber/dataset/dataset_exactlyone/dataset_exactlyone_complet.fa \
  --output-model /gstock/user/bieber/entrainement/model_exactlyone.pth \
  --epochs 5 \
  --esmc-model esmc_300m \
  --num-layers 1