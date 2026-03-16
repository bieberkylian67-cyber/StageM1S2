#Ce script permet de découper en plusieurs fichiers un gros fichier qui contient beaucoup de séquences fasta
import sys

# Input du script
nombre_paquets = int(sys.argv[1])
nombre_seq = int(sys.argv[2])
fichier_seq = sys.argv[3]

seq_par_paquet = nombre_seq // nombre_paquets
compteur = 0
num_paquet = 1

# Ouvrir le premier fichier de sortie
fichier_sortie = open(f"paquet_{num_paquet}.fa", "w")

# Lire le fichier genomic
with open(fichier_seq, "r") as f:
    for ligne in f:
        if ligne.startswith(">"):
            compteur += 1
            if compteur % 1000 == 0:  # afficher tous les 1000 séquences
                print(f"Traité : {compteur}/{nombre_seq} séquences")
            if compteur > seq_par_paquet * num_paquet and num_paquet < nombre_paquets:
                fichier_sortie.close()
                num_paquet += 1
                fichier_sortie = open(f"paquet_{num_paquet}.fa", "w")
                print(f"Paquet {num_paquet} créé !")
        fichier_sortie.write(ligne)

# Ferme le dernier fichier
fichier_sortie.close()

print(f"Script terminé ! {nombre_paquets} paquets créés")


