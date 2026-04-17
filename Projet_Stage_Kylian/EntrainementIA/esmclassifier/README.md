# Protein Sequence Classification with ESM-C

This project uses [EvolutionaryScale's ESM-C models](https://github.com/evolutionaryscale/esm) to classify protein sequences as "correct" or "error".

## Installation

### Prerequisites

Load **uv** using:

```bash
module load python/uv
```

### Setup

1. Clone or download the project.
2. Create a virtual environment and install dependencies:, then activate it:

   ```bash
   uv sync
   . venv/bin/activate
   ```

## Usage

### Training

```bash
python train.py \
  --sequence-length 512 \ # Maximum protein sequence length
  --correct-dataset data/correct.fasta \
  --error-dataset data/error.fasta \
  --output-model model.pth \
  --epochs 5 \
  --esmc-model esmc_300m \
  --num-layers 1
```

The script saves the best weights based on validation accuracy.

### Inference

```bash
python infer.py \
  --input-fasta data/test.fasta \
  --model-path model.pth \
  --esmc-model esmc_300m \
  --num-layers 1
```
