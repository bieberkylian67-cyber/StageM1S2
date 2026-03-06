
import sys
import requests

#Fonction qui permet d'extraire le fichier uniprotKB : 
def extracteur_fichierJSON(id_swissprot):
    url_uniprotkb = f"https://rest.uniprot.org/uniprotkb/{id_swissprot}.json"
    fichier_uniprotkb = requests.get(url_uniprotkb)
    fichier_uniprotkb_JSON = fichier_uniprotkb.json() # on utilise le format .json car c'est un dictionnaire donc plus simple de l'utilsier que .txt
    return fichier_uniprotkb_JSON



#Fonction qui permet de trouver les id ENSEMBL qui est la canonique : 
def ID_ENSEMBL_Canonique(resultat_fichierJSON):
    for comment in resultat_fichierJSON["comments"]:
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
            if cross_ref['isoformId'] == resultat_id:
                liste_ID_Transcript.append(cross_ref["id"])
    return liste_ID_Transcript

if __name__ == "__main__":
    id_swissprot = sys.argv[1]
    resultat_fichierJSON = extracteur_fichierJSON(id_swissprot)
    resultat_ID = ID_ENSEMBL_Canonique(resultat_fichierJSON)
    resultat_LISTE = Id_ENSEMBL_Transcript(resultat_fichierJSON, resultat_ID)
    print(resultat_LISTE)

    