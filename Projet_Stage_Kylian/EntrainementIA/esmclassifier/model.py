from esm.models.esmc import ESMC
from esm.sdk.api import ESMProtein, LogitsConfig
import torch
import torch.nn as nn


class ESMClassifier(nn.Module):
    AVAILABLE_MODELS = ["esmc_300m", "esmc_600m"]
    EMBEDDING_DIMS = {"esmc_300m": 960, "esmc_600m": 1152}

    def __init__(
        self,
        esmc_model: str,
        freeze_backbone: bool = True,
        hidden_dim: int = 256,
        num_layers: int = 1,
        dropout: float = 0.1,
    ):
        super().__init__()
        self.esmc_model = esmc_model
        self.freeze_backbone = freeze_backbone
        self.embedding_dim = self.EMBEDDING_DIMS.get(esmc_model, 960)
        self.hidden_dim = hidden_dim
        self.num_layers = int(num_layers)
        self.dropout = float(dropout)

        if self.num_layers <= 0:
            self.head = nn.Linear(self.embedding_dim, 2)
        else:
            layers = []
            # first layer (embedding_dim -> hidden_dim)
            layers.append(nn.Linear(self.embedding_dim, self.hidden_dim))
            layers.append(nn.ReLU())
            layers.append(nn.Dropout(self.dropout))
            # hidden layers (hidden_dim -> hidden_dim)
            for _ in range(self.num_layers - 1):
                layers.append(nn.Linear(self.hidden_dim, self.hidden_dim))
                layers.append(nn.ReLU())
                layers.append(nn.Dropout(self.dropout))
            # final classifier
            layers.append(nn.Linear(self.hidden_dim, 2))
            self.head = nn.Sequential(*layers)

        # Load backbone
        self.backbone = self._load_pretrained_backbone(self.esmc_model)

        if self.freeze_backbone:
            for param in self.backbone.parameters():
                param.requires_grad = False

    def _load_pretrained_backbone(self, model_name):
        try:
            res = ESMC.from_pretrained(model_name)
            if isinstance(res, (tuple, list)) and len(res) == 2:
                model_backbone, alphabet = res
            else:
                model_backbone = res
            return model_backbone
        except Exception as e:
            raise RuntimeError(f"Failed to load ESMC backbone '{self.esmc_model}'")

    def forward(self, embeddings):
        # embeddings: (batch, embedding_dim)
        return self.head(embeddings)

    def get_embeddings(self, seqs, device):
        proteins = [ESMProtein(sequence=seq) for seq in seqs]
        protein_tensors = [self.backbone.encode(protein) for protein in proteins]
        logits_outputs = [
            self.backbone.logits(tensor, LogitsConfig(return_embeddings=True))
            for tensor in protein_tensors
        ]
        embeddings = [output.embeddings for output in logits_outputs]
        # Assuming embeddings are mean-pooled per sequence
        pooled_embeddings = torch.stack(
            [emb.mean(dim=1).squeeze(0) for emb in embeddings]
        )
        return pooled_embeddings.to(device)
    
    def predict(self, seqs, device):
        self.eval()
        with torch.no_grad():
            embeddings = self.get_embeddings(seqs, device)
            logits = self.head(embeddings)
            probs = torch.softmax(logits, dim=1)
            pred_classes = logits.argmax(dim=1).cpu().numpy()
            prob_errors = probs[:, 1].cpu().numpy()
        return pred_classes, prob_errors

    def train_epoch(self, data_loader, optimizer, loss_fn, device):
        self.train()
        total_loss = 0.0
        correct = 0
        total = 0
        for batch in data_loader:
            seqs = [b[1] for b in batch]
            labels = torch.tensor(
                [b[2] for b in batch], dtype=torch.long, device=device
            )
            embeddings = self.get_embeddings(seqs, device)
            logits = self(embeddings)
            loss = loss_fn(logits, labels)
            optimizer.zero_grad()
            loss.backward()
            optimizer.step()
            total_loss += loss.item() * labels.size(0)
            preds = logits.argmax(dim=1)
            correct += (preds == labels).sum().item()
            total += labels.size(0)
        acc = correct / total if total > 0 else 0.0
        return total_loss / total if total > 0 else 0.0, acc

    def evaluate_epoch(self, data_loader, loss_fn, device):
        self.eval()
        total_loss = 0.0
        correct = 0
        total = 0
        with torch.no_grad():
            for batch in data_loader:
                seqs = [b[1] for b in batch]
                labels = torch.tensor(
                    [b[2] for b in batch], dtype=torch.long, device=device
                )
                embeddings = self.get_embeddings(seqs, device)
                logits = self(embeddings)
                loss = loss_fn(logits, labels)
                total_loss += loss.item() * labels.size(0)
                preds = logits.argmax(dim=1)
                correct += (preds == labels).sum().item()
                total += labels.size(0)
        if total == 0:
            return 0.0, 0.0
        return total_loss / total, correct / total
