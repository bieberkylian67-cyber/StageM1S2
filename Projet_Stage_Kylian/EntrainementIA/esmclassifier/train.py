#!/usr/bin/env python3

import argparse
import random
import os
import numpy as np
import torch
from torch.utils.data import Dataset, DataLoader
from typing import List, Tuple
from sklearn.model_selection import train_test_split
from model import ESMClassifier


AMINO_ACIDS = set("ACDEFGHIKLMNPQRSTVWYBXZJUO")
AVAILABLE_MODELS = ESMClassifier.AVAILABLE_MODELS


def read_fasta(path: str, max_length: int) -> List[Tuple[str, str]]:
    records = []
    name = None
    seq_lines = []
    with open(path) as f:
        for line in f:
            line = line.strip()
            if not line:
                continue
            if line.startswith(">"):
                if name is not None:
                    seq = "".join(seq_lines).upper()
                    seq = "".join([c if c in AMINO_ACIDS else "X" for c in seq])
                    if len(seq) > max_length:
                        continue # skip sequences too long
                    records.append((name, seq))
                name = line[1:].strip()
                seq_lines = []
            else:
                seq_lines.append(line)
    if name is not None:
        seq = "".join(seq_lines).upper()
        seq = "".join([c if c in AMINO_ACIDS else "X" for c in seq])
        if len(seq) > max_length:
            seq = seq[:max_length]
        records.append((name, seq))
    return records


class SeqDataset(Dataset):
    def __init__(self, items: List[Tuple[str, str, int]]):
        self.items = items

    def __len__(self):
        return len(self.items)

    def __getitem__(self, idx):
        return self.items[idx]


def main():
    print("Starting training script", flush=True)
    parser = argparse.ArgumentParser()
    parser.add_argument("--sequence-length", type=int, required=True)
    parser.add_argument("--correct-dataset", type=str, required=True)
    parser.add_argument("--error-dataset", type=str, required=True)
    parser.add_argument("--seed", type=int, default=0)
    parser.add_argument("--output-model", type=str, required=True)
    parser.add_argument("--batch-size", type=int, default=32)
    parser.add_argument("--epochs", type=int, default=5)
    parser.add_argument(
        "--num-layers",
        type=int,
        default=1,
        help="number of hidden layers in classifier head",
    )
    parser.add_argument(
        "--esmc-model",
        type=str,
        default="esmc_300m",
        choices=AVAILABLE_MODELS,
        help=f"ESMC pretrained model to use ({', '.join(AVAILABLE_MODELS)})",
    )
    args = parser.parse_args()

    SEQUENCE_LENGTH = args.sequence_length
    CORRECT_DATASET = args.correct_dataset
    ERROR_DATASET = args.error_dataset
    SEED = args.seed
    OUTPUT_MODEL = args.output_model
    BATCH_SIZE = args.batch_size
    EPOCHS = args.epochs
    ESMC_MODEL = args.esmc_model

    # Check parameters
    assert SEQUENCE_LENGTH > 0, "Sequence length must be strictly positive"
    assert os.access(CORRECT_DATASET, os.R_OK), (
        f"Cannot read correct dataset file: {CORRECT_DATASET}"
    )
    assert os.access(ERROR_DATASET, os.R_OK), (
        f"Cannot read error dataset file: {ERROR_DATASET}"
    )
    assert EPOCHS > 0, "Number of epochs must be strictly positive"
    assert ESMC_MODEL in AVAILABLE_MODELS, (
        f"ESMC model must be one of {AVAILABLE_MODELS}"
    )

    DEVICE = torch.device("cuda" if torch.cuda.is_available() else "cpu")

    # Set seeds
    random.seed(SEED)
    np.random.seed(SEED)
    torch.manual_seed(SEED)
    if torch.cuda.is_available():
        torch.cuda.manual_seed_all(SEED)

    # Read datasets
    correct = read_fasta(CORRECT_DATASET, SEQUENCE_LENGTH)
    error = read_fasta(ERROR_DATASET, SEQUENCE_LENGTH)
    print(f"Read {len(correct)} correct and {len(error)} error sequences", flush=True)

    items = []
    for name, seq in correct:
        items.append((name, seq, 0))
    for name, seq in error:
        items.append((name, seq, 1))

    train_items, temp_items = train_test_split(items, test_size=0.2, random_state=SEED)
    val_items, test_items = train_test_split(
        temp_items, test_size=0.5, random_state=SEED
    )
    print(
        f"Split sizes: train={len(train_items)}, val={len(val_items)}, test={len(test_items)}",
        flush=True,
    )

    # Create datasets and data loaders
    train_dataset = SeqDataset(train_items)
    val_dataset = SeqDataset(val_items)
    test_dataset = SeqDataset(test_items)
    train_loader = DataLoader(
        train_dataset, batch_size=BATCH_SIZE, shuffle=True, collate_fn=lambda x: x
    )
    val_loader = DataLoader(
        val_dataset, batch_size=BATCH_SIZE, shuffle=False, collate_fn=lambda x: x
    )
    test_loader = DataLoader(
        test_dataset, batch_size=BATCH_SIZE, shuffle=False, collate_fn=lambda x: x
    )

    # Load model
    model = ESMClassifier(
        esmc_model=ESMC_MODEL,
        num_layers=args.num_layers,
    )
    model.to(DEVICE)
    print(f"Using model {ESMC_MODEL} on device {DEVICE}", flush=True)

    # Loss and optimizer
    loss_fn = torch.nn.CrossEntropyLoss()
    optimizer = torch.optim.Adam(model.parameters(), lr=1e-3)

    # Training loop
    best_val_acc = 0.0
    for epoch in range(EPOCHS):
        train_loss, train_acc = model.train_epoch(
            train_loader, optimizer, loss_fn, DEVICE
        )
        val_loss, val_acc = model.evaluate_epoch(val_loader, loss_fn, DEVICE)
        print(
            f"Epoch {epoch + 1}/{EPOCHS}: Train Loss={train_loss:.4f}, Train Acc={train_acc:.4f}, Val Loss={val_loss:.4f}, Val Acc={val_acc:.4f}",
            flush=True,
        )

        if val_acc > best_val_acc:
            best_val_acc = val_acc
            torch.save(model.state_dict(), OUTPUT_MODEL)

    # Final evaluation on test set
    test_loss, test_acc = model.evaluate_epoch(test_loader, loss_fn, DEVICE)
    print(f"Test Loss={test_loss:.4f}, Test Acc={test_acc:.4f}", flush=True)
    print("Finished training", flush=True)
    print(
        f"Split sizes: train={len(train_items)}, val={len(val_items)}, test={len(test_items)}",
        flush=True,
    )


if __name__ == "__main__":
    main()
