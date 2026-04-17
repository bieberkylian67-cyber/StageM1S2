#Projet ESM-C
Ce projet contient les scripts python permettant d'obtenir le dataset contenant les erreurs et aussi les 
scripts pour faire l'entrainement sur ESM classifier 

##Prérequis
Le script dans le répertoire SCRIPT_ENSEMBL_TO_SEQ utilise la version de python 3.9 il faut installer cet environnement --> environnement_stage.yml

Tous les autres scripts utilisent la version 3.12 de python --> environnement de base 

Le Gene predictor utilisé est Augustus 

Il faut adapter les chemins des fichiers input et ouput pour pouvoir utiliser les scripts.

##Dataset 
Dataset incorrect = " dataset_partial_complet.fa "
Chemin : /gstock/user/bieber/dataset/dataset_partial/dataset_partial_complet.fa

Dataset correct version 1 = " dataset_proteine.fa "  --> version avec les isoformes
Dataset correct version 2 = "dataset_proteine_pasdoublons.fa " --> version sans les isoformes
Chemin : /gstock/user/bieber
