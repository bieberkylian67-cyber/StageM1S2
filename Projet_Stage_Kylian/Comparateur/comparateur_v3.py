#Ce script permet de comparer les positions des exons prédits par Augustus et les positions des exons de ENSEMBL

import sys 
import json
#Input 
# Argument 1 : fichier qui contient les prédictions de Augustus
# Argument 2 : fichier ENSEMBL qui contient les bonnes positions des exons 

#Output un fichier dataset

#Code pour les erreurs : 
#N1 N-terminal extension
#N2 N-terminal deletion
#C1 C-terminal extension
#C2 C-terminal deletion
#1  insertion
#2  deletion
#3  inconsistent segment


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

def filtre_mauvaise_prediction(dico_conversion, dico_ensembl, dico_erreurs, dataset_bad_predictions):

    with open(dataset_bad_predictions, "w") as fichier:
        for id_transcript, total_genes in dico_conversion.items():
            nombre_gene_total = len(total_genes)
            compteur = 1 
            exon_map_ensembl = dico_ensembl.get(id_transcript, {}).get("Exon", [])
            nombre_exons_ensembl = len(exon_map_ensembl)
            erreurs_header = " ".join(dico_erreurs.get(id_transcript, []))

            for genes, exons in total_genes.items():
                if nombre_gene_total > 1:
                    fin_header = f"_g{compteur}"
                else:
                    fin_header = ""

                header_final = f">{id_transcript}{fin_header} {erreurs_header}\n"
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

def annotations(dico_conversion, dico_ensembl):
    dico_erreurs = {}
    for id_transcript, gene_augustus in dico_conversion.items():
        rassembler_exons = []
        liste_exon_ensembl = []
        erreurs = []
        
        # 1. On prépare Augustus
        for g, cle in gene_augustus.items():
            rassembler_exons.extend(cle["exons"])
        rassembler_exons.sort()

        # 2. On prépare Ensembl
        exons_ensembl = dico_ensembl[id_transcript]["Exon"]
        for exon in exons_ensembl:
            start_ensembl = int(exon["start"])
            end_ensembl = int(exon["end"])
            tuple_ensembl = (start_ensembl, end_ensembl)
            liste_exon_ensembl.append(tuple_ensembl)

        # MEMOIRE EXONS  (Pour les Délétions)
        coches_ensembl = [False] * len(liste_exon_ensembl)

        # 3. On compare chaque exon Augustus
        for start, end in rassembler_exons:
            start_augustus = int(start)
            end_augustus = int(end)
            trouve = False 

            # On fouille la liste Ensembl pour cet exon précis
            for i, (s_ens, e_ens) in enumerate(liste_exon_ensembl):
                
                # CAS MATCH
                if start_augustus == s_ens and end_augustus == e_ens:
                    trouve = True 
                    coches_ensembl[i] = True # On coche !
                    break
                
                # CAS MISMATCH (Chevauchement)

                elif start_augustus <= e_ens and end_augustus >= s_ens:
                    coches_ensembl[i] = True
                    trouve = True


                    if start_augustus < s_ens and end_augustus > e_ens:
                        # Augustus trop long des 2 côtés
                        erreurs.append(f"N1({start_augustus}, {s_ens})")
                        erreurs.append(f"C1({e_ens}, {end_augustus})")

                    elif start_augustus > s_ens and end_augustus < e_ens:
                        # Augustus trop court des 2 côtés
                        erreurs.append(f"N2({s_ens}, {start_augustus})")
                        erreurs.append(f"C2({end_augustus}, {e_ens})")

                    elif start_augustus == s_ens and end_augustus > e_ens:
                        # Start identique, end trop long
                        erreurs.append(f"C1({e_ens}, {end_augustus})")
                        

                    elif start_augustus < s_ens and end_augustus == e_ens:
                        # Start trop long, end identique
                        erreurs.append(f"N1({start_augustus}, {s_ens})")
                        
                    elif start_augustus > s_ens and end_augustus == e_ens:
                        # Start trop court, end identique
                        erreurs.append(f"N2({s_ens}, {start_augustus})")

                    elif start_augustus == s_ens and end_augustus < e_ens:
                        # Start identique, end trop court
                        erreurs.append(f"C2({end_augustus}, {e_ens})")

                    elif start_augustus < s_ens and end_augustus < e_ens:
                        # Start trop long, end trop court 
                        erreurs.append(f"N1({start_augustus}, {s_ens})")
                        erreurs.append(f"C2({end_augustus}, {e_ens})")

                    elif start_augustus > s_ens and end_augustus > e_ens:
                        # Start trop court, end trop long
                        erreurs.append(f"N2({s_ens}, {start_augustus})")
                        erreurs.append(f"C1({e_ens}, {end_augustus})")

                    break
            
            if not trouve:
                erreurs.append(f"1({start_augustus}, {end_augustus})")
            
        #4. DERNIÈRE ÉTAPE : LES DÉLÉTIONS
        # On regarde quels exons Ensembl n'ont jamais été "cochés"
        for i, (s_ens, e_ens) in enumerate(liste_exon_ensembl):
            if coches_ensembl[i] == False:
                erreurs.append(f"2({s_ens}, {e_ens})")
            
        dico_erreurs[id_transcript] = erreurs
    return dico_erreurs 

if __name__ == "__main__":
    fichier_input_augustus = sys.argv[1]
    fichier_input_ensembl = sys.argv[2]
    dataset_bad_predictions = sys.argv[3]
        
    # 1. On convertit les positions d'Augustus
    dico_conv = convertir_positions(fichier_input_augustus, fichier_input_ensembl)

    # 2. On récupère le dictionnaire Ensembl
    dico_ensembl = dictionnaire_ensembl(fichier_input_ensembl)

    # 3. On génère le fichier de statistiques/annotations
    dico_erreurs = annotations(dico_conv, dico_ensembl)

    # 4. (Optionnel) Ton filtre de mauvaises prédictions
    filtre_mauvaise_prediction(dico_conv, dico_ensembl, dico_erreurs, dataset_bad_predictions)

    print("Script terminé")

