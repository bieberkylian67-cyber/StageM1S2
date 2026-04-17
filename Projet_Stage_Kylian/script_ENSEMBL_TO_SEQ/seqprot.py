import sys
from Bio import SeqIO
from Bio.SeqRecord import SeqRecord


fichier_entree = sys.argv[1]
fichier_sortie = sys.argv[2]

donnees_proteines = []


for record in SeqIO.parse(fichier_entree, "fasta"):
    
    #(multiples de 3)
    if len(record.seq) % 3 == 0:
        
        sequence_traduite = record.seq.translate(to_stop=True)
        description_maj = record.description.replace("CDS", "protein")
        
        enregistrement = SeqRecord(
            sequence_traduite,
            id=record.id,
            description=description_maj
        )
        
        donnees_proteines.append(enregistrement)

# Écriture du résultat
SeqIO.write(donnees_proteines, fichier_sortie, "fasta")