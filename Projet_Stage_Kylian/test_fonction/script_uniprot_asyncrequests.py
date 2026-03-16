#Ce script permet de récupérer les IDs Ensembl de plusieurs IDs swissprot
#Input : 
# argument 1 : fichier qui contient tous les ids swissprot
# argument 2 : le fichier de sortie 

#Output : fichier en format table avec en 1er colonne l'id swissprot (AC) et en 2ème colonne l'id ENSEMBL

#Importation librairie nécessaire pour le script
import sys
import requests
import time 

#Fonction qui permet d'extraire le fichier uniprotKB mais on utilise API batch pour faire des requêtes groupées
def extracteur_fichierJSON(liste_ids):
    url = "https://rest.uniprot.org/uniprotkb/accessions"
    params = {"accessions": ",".join(liste_ids), "format": "json"}
    reponse = requests.get(url, params=params)
    print(reponse.status_code)
    return reponse.json()["results"]


    #Fonction qui permet de trouver les id ENSEMBL qui est la canonique : 
def ID_ENSEMBL_Canonique(resultat_fichierJSON):
    for comment in resultat_fichierJSON.get("comments",[]):
        if comment["commentType"] == "ALTERNATIVE PRODUCTS":
            for isoform in comment["isoforms"]:
                if isoform['isoformSequenceStatus'] == "Displayed":
                    resultat_id = isoform["isoformIds"][0]
                    return resultat_id


    #Fonction qui permet de trouver les id ENSEMBL Transcript : 
def Id_ENSEMBL_Transcript(resultat_fichierJSON, resultat_id):
    liste_ID_Transcript = []
    for cross_ref in resultat_fichierJSON["uniProtKBCrossReferences"]:
        if cross_ref["database"] == "Ensembl":
            if cross_ref.get('isoformId') == resultat_id: #get() évite une KeyError si 'isoformId' est absent dans l'entrée
                juste_id = cross_ref["id"].split(".") # permet d'enlever les numéro de versions sur les id ENSEMBL
                liste_ID_Transcript.append(juste_id[0])
    return liste_ID_Transcript


if __name__ == "__main__":
    fichier_swissprot = sys.argv[1]
    chemin_fichier_sortie = sys.argv[2]
    with open(fichier_swissprot, "r") as fichier:
        contenu = fichier.readlines()
        
    liste_id_swissprot = []
        
    for id in contenu:
        a = id.strip("\n")
        liste_id_swissprot.append(a)
        
    with open(chemin_fichier_sortie, "w") as fichier_sortie:
        fichier_sortie.write("AC\tENSEMBL\n")
        
        for i in range(0, len(liste_id_swissprot), 200):
            groupe = liste_id_swissprot[i:i+200]
            print(f"Groupe {i//200 + 1}/{len(liste_id_swissprot)//200 + 1}")
            resultats = extracteur_fichierJSON(groupe)
            for proteine in resultats:
                ac = proteine["primaryAccession"]
                canonique = ID_ENSEMBL_Canonique(proteine)
                transcripts = Id_ENSEMBL_Transcript(proteine, canonique)
                for t in transcripts:
                    fichier_sortie.write(f"{ac}\t{t}\n")
            time.sleep(1)
