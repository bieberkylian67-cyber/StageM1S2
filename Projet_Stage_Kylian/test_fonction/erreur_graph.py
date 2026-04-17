import sys
import re
import matplotlib.pyplot as plt
from collections import Counter

input_file = sys.argv[1]

error_counts = Counter()

with open(input_file) as f:
    for line in f:
        if line.startswith(">"):
            # Cherche N1, N2, C1, C2, 1, 2
            codes = re.findall(r'\b(N1|N2|C1|C2|(?<![NC])1(?!\d)|(?<![NC])2(?!\d))\b', line)
            for code in codes:
                error_counts[code] += 1

# Noms lisibles pour le graphique
labels_map = {
    "N1": "N-term deletion",
    "N2": "N-term insertion", 
    "C1": "C-term deletion",
    "C2": "C-term insertion",
    "1": "deletion",
    "2": "insertion"
}

print("Distribution des erreurs :")
total = sum(error_counts.values())
for error, count in sorted(error_counts.items()):
    label = labels_map.get(error, error)
    print(f"  {label}: {count} ({count/total*100:.1f}%)")

labels = [labels_map.get(k, k) for k in error_counts.keys()]
values = list(error_counts.values())

plt.figure(figsize=(10, 6))
bars = plt.bar(labels, values, color="steelblue")
plt.title("Distribution of prediction errors in the dataset")
plt.xlabel("Error type")
plt.ylabel("Number of occurrences")
plt.xticks(rotation=20)
for bar, val in zip(bars, values):
    plt.text(bar.get_x() + bar.get_width()/2, bar.get_height() + 50,
             str(val), ha='center', fontsize=10)
plt.tight_layout()
plt.savefig("error_distribution.png", dpi=300)
print("Graphique sauvegardé : error_distribution.png")