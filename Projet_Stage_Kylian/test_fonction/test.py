# Ce script permet de récupérer les IDs Ensembl de plusieurs IDs swissprot
# Input : 
# argument 1 : fichier qui contient tous les ids swissprot
# argument 2 : le fichier de sortie 

# Output : fichier en format table avec en 1er colonne l'id swissprot (AC) et en 2ème colonne l'id ENSEMBL unique

import sys
import requests
import time 

# Fonction qui permet d'extraire le fichier uniprotKB via l'API Batch
def extracteur_fichierJSON(liste_ids):
    url = "https://rest.uniprot.org/uniprotkb/accessions"
    params = {"accessions": ",".join(liste_ids), "format": "json"}
    reponse = requests.get(url, params=params)
    if reponse.status_code != 200:
        print(f"Erreur API: {reponse.status_code}")
        return []
    return reponse.json().get("results", [])

# Fonction qui permet de trouver l'ID de l'isoforme affichée (canonique)
def ID_ENSEMBL_Canonique(resultat_fichierJSON):
    for comment in resultat_fichierJSON.get("comments", []):
        if comment.get("commentType") == "ALTERNATIVE PRODUCTS":
            for isoform in comment.get("isoforms", []):
                if isoform.get('isoformSequenceStatus') == "Displayed":
                    return isoform["isoformIds"][0]
    return None

# Fonction qui permet de trouver les IDs Ensembl Transcript correspondants
def Id_ENSEMBL_Transcript(resultat_fichierJSON, resultat_id):
    liste_ID_Transcript = []
    # Parfois 'uniProtKBCrossReferences' peut manquer
    for cross_ref in resultat_fichierJSON.get("uniProtKBCrossReferences", []):
        if cross_ref["database"] == "Ensembl":
            # On vérifie si l'isoforme correspond ou s'il n'y a pas d'isoforme spécifique listé
            if cross_ref.get('isoformId') == resultat_id or resultat_id is None:
                juste_id = cross_ref["id"].split(".") 
                liste_ID_Transcript.append(juste_id[0])
    return liste_ID_Transcript

if __name__ == "__main__":
    if len(sys.argv) < 3:
        print("Usage: python script.py <input_swissprot> <output_table>")
        sys.exit(1)

    fichier_swissprot = sys.argv[1]
    chemin_fichier_sortie = sys.argv[2]

    # Lecture des IDs d'entrée
    with open(fichier_swissprot, "r") as fichier:
        liste_id_swissprot = [line.strip() for line in fichier if line.strip()]

    with open(chemin_fichier_sortie, "w") as fichier_sortie:
        fichier_sortie.write("AC\tENSEMBL\n")
        
        # Requêtes par groupes de 200 pour l'API UniProt
        for i in range(0, len(liste_id_swissprot), 200):
            groupe = liste_id_swissprot[i:i+200]
            print(f"Traitement du groupe {i//200 + 1}/{(len(liste_id_swissprot)-1)//200 + 1}...")
            
            try:
                resultats = extracteur_fichierJSON(groupe)
                for proteine in resultats:
                    ac = proteine["primaryAccession"]
                    canonique = ID_ENSEMBL_Canonique(proteine)
                    transcripts = Id_ENSEMBL_Transcript(proteine, canonique)
                    
                    # --- MODIFICATION MAJEURE ICI ---
                    # Au lieu de boucler sur tous les transcrits, on n'en prend qu'un seul.
                    if transcripts:
                        t_unique = transcripts[0] 
                        fichier_sortie.write(f"{ac}\t{t_unique}\n")
                    # --------------------------------
            except Exception as e:
                print(f"Erreur lors du traitement d'un groupe : {e}")
            
            # Pause pour respecter les limites de l'API
            time.sleep(1)

    print(f"Extraction terminée. Fichier sauvegardé sous : {chemin_fichier_sortie}")