#Ce script permet de récupérer l'id ENSEMBL mais fonctionne que pour un id swissprot

#Input 
#arg1 : id_swissprot
#arg2 : chemin fichier de sortie

 #Output
 #Sort en sortie un fichier au format tsv avec en 1er colonne le AC et en deuxièmme colonne l'ID
 

#Importation librairie nécessaire pour le script
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
            if cross_ref.get('isoformId') == resultat_id: #get() évite une KeyError si 'isoformId' est absent dans l'entrée
                liste_ID_Transcript.append(cross_ref["id"])
    return liste_ID_Transcript


#Fonction qui permet d'écrire le résultat dans un fichier :
def ecrire_fichier_resultat(id_swissprot, resultat_LISTE, chemin_fichier_sortie):
    with open (chemin_fichier_sortie, "w") as fichier:
        fichier.write("AC\tENSEMBL\n") #écrit le header du fichier 
        for i in resultat_LISTE:
            fichier.write(f"{id_swissprot}\t{i}\n")


#Fonction finale qui utilise toutes les fonctions précédentes : 
def fonction_finale(id_swissprot, chemin_fichier_sortie):
    resultat_fichierJSON = extracteur_fichierJSON(id_swissprot)
    resultat_ID = ID_ENSEMBL_Canonique(resultat_fichierJSON)
    resultat_LISTE = Id_ENSEMBL_Transcript(resultat_fichierJSON, resultat_ID)
    ecrire_fichier_resultat(id_swissprot, resultat_LISTE, chemin_fichier_sortie)


if __name__ == "__main__":
    id_swissprot = sys.argv[1]
    chemin_fichier_sortie = sys.argv[2]
    fonction_finale(id_swissprot, chemin_fichier_sortie)