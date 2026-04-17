#!/usr/bin/env python3

import argparse
import torch
from model import ESMClassifier
from Bio import SeqIO
from sklearn.metrics import confusion_matrix 

def main():
    y_true = []
    y_pred = []

    parser = argparse.ArgumentParser()
    parser.add_argument("--correct-fasta", required=True)
    parser.add_argument("--error-fasta", required=True)
    parser.add_argument("--model-path", default="out.model", help="Path to saved model")
    parser.add_argument("--esmc-model", default="esmc_300m", help="ESMC model type")
    parser.add_argument(
        "--num-layers", type=int, default=1, help="Number of layers in classifier"
    )
    parser.add_argument(
        "--max-length", type=int, default=512, help="Maximum sequence length"
    )
    args = parser.parse_args()

    # Load the model
    DEVICE = "cuda" if torch.cuda.is_available() else "cpu"

    model = ESMClassifier(
        esmc_model=args.esmc_model,
        num_layers=args.num_layers,
    )
    state_dict = torch.load(args.model_path, map_location="cpu")
    model.load_state_dict(state_dict)
    model.to(DEVICE)
    model.eval()

    for record in SeqIO.parse(args.correct_fasta, "fasta"):
        name = record.id
        sequence = str(record.seq).upper()
        if len(sequence) > args.max_length:
            print(f"Skipping {name} (length {len(sequence)} exceeds max length)")
            continue
        # Get embeddings
        pred_classes, prob_errors = model.predict([sequence], DEVICE)
        y_true.append(0)
        y_pred.append(pred_classes[0])

        label = "CORRECT" if pred_classes[0] == 0 else "ERROR" # class 0 = correct
        print(
            f"{name}: Predicted {label}, probability : {prob_errors[0]:.4f}"
        )

    for record in SeqIO.parse(args.error_fasta, "fasta"):
        name = record.id
        sequence = str(record.seq).upper()
        if len(sequence) > args.max_length:
            print(f"Skipping {name} (length {len(sequence)} exceeds max length)")
            continue
        # Get embeddings
        pred_classes, prob_errors = model.predict([sequence], DEVICE)
        y_true.append(1)
        y_pred.append(pred_classes[0])

        label = "CORRECT" if pred_classes[0] == 0 else "ERROR" # class 0 = correct
        print(
            f"{name}: Predicted {label}, probability : {prob_errors[0]:.4f}"
        )

    matrice_confusion = confusion_matrix(y_true, y_pred)
    print(matrice_confusion)
    TN = matrice_confusion[0,0]
    TP = matrice_confusion[1,1]
    FN = matrice_confusion[1,0]
    FP = matrice_confusion[0,1]
    
    recall = TP / (TP + FN)
    print(f"Recall = {recall}")

    specificity = TN / (TN + FP)
    print(f"Specificity = {specificity}")    

if __name__ == "__main__":
    main()
