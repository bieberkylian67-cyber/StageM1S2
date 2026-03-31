#Ce script permet de comparer les positions des exons prédits par Augustus et les positions des exons de ENSEMBL

import sys 
import json
import Bio
#Input 
# Argument 1 : fichier qui contient les prédictions de Augustus
# Argument 2 : fichier ENSEMBL qui contient les bonnes positions des exons 

#Output 

fichier_input_augustus = sys.argv[1]
fichier_input_ensembl = sys.argv[2]
dataset_bad_predictions = sys.argv[3]
fichier_annotation = sys.argv[4]

def convertir_gff3_to_dict(fichier_input_augustus):
    with open (fichier_input_augustus, "r") as fichier:
        dictionnaire_augustus = {}
        id_memoire = None 
        lecture_proteine = False
        gene_memoire = None 

        for ligne in fichier:
            if ligne.startswith("#") and "name =" not in ligne and "start gene" not in ligne and "protein sequence =" not in ligne and lecture_proteine == False:
                continue
            elif "name = " in ligne:
                id_memoire = ligne.split("name =")[1].split(")")[0].strip()
                dictionnaire_augustus[id_memoire] = {}
            elif "start gene" in ligne:
                gene_memoire = ligne.split()[3]
                dictionnaire_augustus[id_memoire][gene_memoire] = {"exons": [], "sequence": "","strand":None}
            elif "CDS" in ligne and not ligne.startswith("#"):
                liste_ligne = ligne.split("\t")
                debut = int(liste_ligne[3])
                fin = int(liste_ligne[4])
                strand_symbole = liste_ligne[6]
                strand_num = 1 if strand_symbole == "+" else -1
                tuple_position = (debut, fin)
                dictionnaire_augustus[id_memoire][gene_memoire]["exons"].append((tuple_position))
                dictionnaire_augustus[id_memoire][gene_memoire]["strand"] = strand_num
            elif "protein sequence =" in ligne: 
                lecture_proteine = True
                sequence_prot = ligne.split("[")[1]
                
                if "]" in sequence_prot:
                    dictionnaire_augustus[id_memoire][gene_memoire]["sequence"] += sequence_prot.split("]")[0].strip()
                    lecture_proteine = False 
                else:
                    dictionnaire_augustus[id_memoire][gene_memoire]["sequence"] += sequence_prot.strip()
            elif lecture_proteine == True and ligne.startswith("#"):
                suite_seq = ligne.replace("#","").strip()
                if "]" in suite_seq:
                    dictionnaire_augustus[id_memoire][gene_memoire]["sequence"] += suite_seq.split("]")[0].strip()
                    lecture_proteine = False
                else:
                    dictionnaire_augustus[id_memoire][gene_memoire]["sequence"] += suite_seq

    return dictionnaire_augustus

def dictionnaire_ensembl(fichier_input_ensembl):
    dictionnaire = {}
    with open (fichier_input_ensembl, "r") as fichier:
        for ligne in fichier:
            convertir = json.loads(ligne)
            dictionnaire.update(convertir)
    return dictionnaire

def convertir_positions(fichier_input_augustus, fichier_input_ensembl):
    dico_conversion = {}
    dictionnaire_augustus = convertir_gff3_to_dict(fichier_input_augustus)
    start = dictionnaire_ensembl(fichier_input_ensembl)
    for transcript_id in dictionnaire_augustus :
        dico_conversion[transcript_id] = {}
        start_ensembl = start[transcript_id]["start"]
        for gene_number, infos_gene in dictionnaire_augustus[transcript_id].items():
            exon_convertirt = []
            strand_aug = infos_gene["strand"]
            for (debut, fin) in infos_gene["exons"]:
                debut_convert = (start_ensembl + debut - 1)
                fin_convert = (start_ensembl + fin - 1)
                exon_convertirt.append((debut_convert, fin_convert))
            dico_conversion[transcript_id][gene_number] = {"exons": exon_convertirt, "sequence":infos_gene["sequence"],"strand": strand_aug}
    return dico_conversion

def filtre_mauvaise_prediction(dico_conversion, dico_ensembl,dataset_bad_predictions):

    with open(dataset_bad_predictions, "w") as fichier:
        for id_transcript, total_genes in dico_conversion.items():
            nombre_gene_total = len(total_genes)
            compteur = 1 
            exon_map_ensembl = dico_ensembl.get(id_transcript, {}).get("Exon", [])
            nombre_exons_ensembl = len(exon_map_ensembl)

            for genes, exons in total_genes.items():
                if nombre_gene_total > 1:
                    fin_header = f"_g{compteur}"
                else:
                    fin_header = ""

                header_final = f">{id_transcript}{fin_header}\n"
                liste_exon_augustus = exons["exons"]
                nombre_exons_augustus = len(liste_exon_augustus)
                compteur += 1 

                if nombre_exons_augustus != nombre_exons_ensembl:
                    fichier.write(header_final)
                    fichier.write(f"{exons['sequence']}\n")

                else:
                    trie_exons_augustus = sorted(liste_exon_augustus)
                    conversion_tuple_ensembl = []

                    for exon in exon_map_ensembl:
                        debut = exon["start"]
                        fin = exon["end"]
                        stock_temporaire = (debut,fin)
                        conversion_tuple_ensembl.append(stock_temporaire)
                    trie_exons_ensembl = sorted(conversion_tuple_ensembl)

                    if trie_exons_augustus == trie_exons_ensembl:
                        continue

                    elif trie_exons_augustus != trie_exons_ensembl:
                        fichier.write(header_final)
                        fichier.write(f"{exons['sequence']}\n")
    return dataset_bad_predictions

def annotations(dico_conversion, dico_ensembl, fichier_annotation):
    with open(fichier_annotation,"w") as fichier:
        header = f"ID\tStart_Augustus\tEnd_Augustus\tType_Erreur\tDecalage_Start\tDecalage_End\n"
        fichier.write(header)
        for id_transcript, gene_augustus in dico_conversion.items():
            rassembler_exons = []
            for g, cle in gene_augustus.items():
                rassembler_exons.extend(cle["exons"])
            rassembler_exons.sort()
            exons_ensembl = dico_ensembl[id_transcript]["Exon"]
            for start, end in rassembler_exons:
                start_augustus = int(start)
                end_augustus = int(end)

                type_erreur = "Insertion"
                decalage_start = start_augustus
                decalage_end = end_augustus
                for exon in exons_ensembl:
                    start_ensembl = int(exon["start"])
                    end_ensembl = int(exon["end"])
                    if start_augustus == start_ensembl and end_augustus == end_ensembl:
                        type_erreur = "Pas d'erreur"
                        decalage_start = 0
                        decalage_end = 0
                        break
                    elif (start_ensembl <= start_augustus <= end_ensembl) or (start_ensembl <= end_augustus <= end_ensembl):
                        type_erreur = "Mismatch"
                        decalage_start = start_augustus - start_ensembl
                        decalage_end = end_augustus - end_ensembl
                        break
                    elif (start_augustus < start_ensembl) and (end_augustus < start_ensembl) or (start_augustus > end_ensembl) and (end_augustus > end_ensembl):
                        type_erreur = "Délétion"
                        decalage_start = start_augustus
                        decalage_end = end_augustus
                        break
                ligne = f"{id_transcript}\t{start_augustus}\t{end_augustus}\t{type_erreur}\t{decalage_start}\t{decalage_end}\n"
                fichier.write(ligne)

dico_conversion =  convertir_positions(fichier_input_augustus, fichier_input_ensembl)
dico_ensembl = dictionnaire_ensembl(fichier_input_ensembl)

#resultat = annotations(dico_conversion, dico_ensembl, fichier_annotation)

print(dico_conversion)
