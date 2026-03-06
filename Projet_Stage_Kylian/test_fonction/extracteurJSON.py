import sys
import requests

id_swissprot = sys.argv[1]

#Permet d'extraire le fichier UniprotKB au format JSON 
def extracteur_fichierJSON(id_swissprot):
    url_uniprotkb = f"https://rest.uniprot.org/uniprotkb/{id_swissprot}.json"
    fichier_uniprotkb = requests.get(url_uniprotkb)
    fichier_uniprotkb_JSON = fichier_uniprotkb.json() # on utilise le format .json car c'est un dictionnaire donc plus simple de l'utilsier que .txt
    return fichier_uniprotkb_JSON

resultat_fichierJSON = extracteur_fichierJSON(id_swissprot)

if __name__ == "__main__":
    print(resultat_fichierJSON)



