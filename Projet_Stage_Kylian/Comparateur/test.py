#Ce script est un test 
from Bio import SeqIO
import sys

# Vérification des arguments
if len(sys.argv) < 3:
    print("Usage: python script.py fichier.fa ID_CIBLE")
    sys.exit(1)

fichier_input = sys.argv[1]
target_id = sys.argv[2]
fichier_sortie = f"extract_{target_id}.fa"

def extraire_et_sauvegarder(input_fa, id_cible, output_fa):
    found = False
    with open(output_fa, "w") as f_out:
        for record in SeqIO.parse(input_fa, "fasta"):
            # On vérifie si l'ID est exactement celui recherché
            if record.id == id_cible:
                # SeqIO.write s'occupe de mettre le >header et la séquence proprement
                SeqIO.write(record, f_out, "fasta")
                found = True
                break # On arrête dès qu'on l'a trouvé pour gagner du temps
    
    if found:
        print(f"Séquence {id_cible} extraite avec succès dans {output_fa}")
    else:
        print(f" ID {id_cible} non trouvé dans le fichier.")

# Lancement script
extraire_et_sauvegarder(fichier_input, target_id, fichier_sortie)
