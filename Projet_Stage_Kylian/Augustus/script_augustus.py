#Ce script permet de lancer augustus sur chaque séquence d'un fichier fasta 
import sys
import os 
# 1. On récupère le paquet envoyé par Slurm (ex: paquet_1.fa)
nom_paquet = sys.argv[1]
nom_resultat = sys.argv[2]
# On crée un nom de fichier temporaire unique pour ce paquet
nom_temp = "/home/bieber/stage/Projet_Stage_Kylian/Augustus/temporaire/temp_" + os.path.basename(nom_paquet)

dictionnaire = {
    # === human ===
    "ENST00": "human",    # Homo sapiens
    "ENSMUS": "human",    # Mus musculus → pas de mouse dans Augustus
    "ENSRNO": "human",    # Rattus norvegicus → pas de rat dans Augustus
    "ENSBTA": "human",    # Bos taurus
    "ENSSSC": "human",    # Sus scrofa
    "ENSCAF": "human",    # Canis lupus familiaris
    "ENSECA": "human",    # Equus caballus
    "ENSEAS": "human",    # Equus asinus
    "ENSFCA": "human",    # Felis catus
    "ENSOAR": "human",    # Ovis aries
    "ENSBGR": "human",    # Bos mutus
    "ENSBIX": "human",    # Bison bison
    "ENSCDR": "human",    # Camelus dromedarius
    "ENSVPA": "human",    # Vicugna pacos
    "ENSLAF": "human",    # Loxodonta africana
    "ENSDNO": "human",    # Dasypus novemcinctus
    "ENSMOD": "human",    # Monodelphis domestica
    "ENSMEU": "human",    # Notamacropus eugenii
    "ENSOAN": "human",    # Ornithorhynchus anatinus
    "ENSUMA": "human",    # Ursus maritimus
    "ENSDLE": "human",    # Delphinapterus leucas
    "ENSTTR": "human",    # Tursiops truncatus
    "ENSNVI": "human",    # Neovison vison
    "ENSPTI": "human",    # Panthera tigris
    "ENSVVU": "human",    # Vulpes vulpes
    "ENSRFE": "human",    # Rhinolophus ferrumequinum
    # Primates
    "ENSPPY": "human",    # Pongo abelii
    "ENSPTR": "human",    # Pan troglodytes
    "ENSPPA": "human",    # Pan paniscus
    "ENSMFA": "human",    # Macaca fascicularis
    "ENSMMU": "human",    # Macaca mulatta
    "ENSMNE": "human",    # Macaca nemestrina
    "ENSGGO": "human",    # Gorilla gorilla
    "ENSPAN": "human",    # Papio anubis
    "ENSCJA": "human",    # Callithrix jacchus
    "ENSSBO": "human",    # Saimiri boliviensis
    "ENSTGE": "human",    # Theropithecus gelada
    "ENSNLE": "human",    # Nomascus leucogenys
    "ENSCSA": "human",    # Cercocebus atys
    "ENSCMM": "human",    # Macaque sp.
    "ENSCJP": "human",    # Aotus nancymaae
    "ENSMIC": "human",    # Microcebus murinus
    "ENSOGA": "human",    # Otolemur garnettii
    "ENSTSY": "human",    # Carlito syrichta
    # Rongeurs → human (souris absent)
    "ENSCGR": "human",    # Cricetulus griseus
    "ENSOCU": "human",    # Oryctolagus cuniculus
    "ENSCPO": "human",    # Cavia porcellus
    "ENSMAU": "human",    # Mesocricetus auratus
    "ENSCHI": "human",    # Chinchilla lanigera
    "ENSMUG": "human",    # Meriones unguiculatus
    "ENSHGL": "human",    # Heterocephalus glaber
    "ENSSTO": "human",    # Ictidomys tridecemlineatus
    "ENSODE": "human",    # Octodon degus
    "ENSMSI": "human",    # Mus spicilegus
    "ENSJJA": "human",    # Jaculus jaculus
    "ENSMOC": "human",    # Microtus ochrogaster
    "ENSNGA": "human",    # Nannospalax galili
    "ENSMLU": "human",    # Myotis lucifugus
    # Reptiles/Amphibiens → human
    "ENSXET": "human",    # Xenopus tropicalis
    "ENSACA": "human",    # Anolis carolinensis
    "ENSPCA": "human",    # Pelodiscus sinensis
    "ENSCPR": "human",    # Crocodylus porosus

    #    chicken 
    "ENSGAL": "chicken",  # Gallus gallus
    "ENSTGU": "chicken",  # Taeniopygia guttata
    "ENSMGA": "chicken",  # Meleagris gallopavo
    "ENSAPL": "chicken",  # Anas platyrhynchos
    "ENSANA": "chicken",  # Anas sp.
    "ENSAME": "chicken",  # Dromaius novaehollandiae

    #    zebrafish 
    "ENSDAR": "zebrafish",  # Danio rerio
    "ENSSSA": "zebrafish",  # Salmo salar
    "ENSTRU": "zebrafish",  # Takifugu rubripes
    "ENSTNI": "zebrafish",  # Tetraodon nigroviridis
    "ENSORL": "zebrafish",  # Oryzias latipes
    "ENSCAR": "zebrafish",  # Carassius auratus
    "ENSOMY": "zebrafish",  # Oncorhynchus mykiss
    "ENSDLA": "zebrafish",  # Dicentrarchus labrax
    "ENSSAU": "zebrafish",  # Sparus aurata
    "ENSONI": "zebrafish",  # Oreochromis niloticus
    "ENSSTU": "zebrafish",  # Salmo trutta
    "ENSHCO": "zebrafish",  # Hippocampus comes
    "ENSCCR": "zebrafish",  # Cyprinus carpio
    "ENSIPU": "zebrafish",  # Ictalurus punctatus
    "ENSELU": "zebrafish",  # Electrophorus electricus

    #    Espèces spécifiques Augustus 
    "ENSLAC": "elephant_shark",   # Latimeria chalumnae → requin le + proche
    "ENSCMI": "elephant_shark",   # Callorhinchus milii
    "ENSCIN": "ciona",            # Ciona intestinalis
    "ENSXMA": "Xiphophorus_maculatus",  # Xiphophorus maculatus
    "ENSPMA": "sealamprey",       # Petromyzon marinus
}

# 3. On lit le paquet
with open(nom_paquet, "r") as f:
    header = ""
    sequence = ""

    for ligne in f:
        ligne = ligne.strip()
        if not ligne: 
            continue

        if ligne.startswith(">"):
            # Si on a une séquence prête, on la traite
            if header != "":
                # Trouver l'espèce (on prend les 6 caractères après le >)
                prefixe = header[1:7]
                espece = dictionnaire.get(prefixe, "human")

                # Créer le fichier temporaire
                with open(nom_temp, "w") as f_tmp:
                    f_tmp.write(header + "\n" + sequence + "\n")

                # Lancer Augustus
                commande = f"augustus --species={espece} {nom_temp} >> {nom_resultat}"
                os.system(commande)

            # Reset pour la séquence suivante
            header = ligne
            sequence = ""
        else:
            sequence = sequence + ligne

    # Traiter la toute dernière séquence du fichier
    if header != "":
        prefixe = header[1:7]
        espece = dictionnaire.get(prefixe, "human")
        with open(nom_temp, "w") as f_tmp:
            f_tmp.write(header + "\n" + sequence + "\n")
        os.system("augustus --species=" + espece + " " + nom_temp + " >> " + nom_resultat)

# 4. On nettoie le fichier temporaire à la fin
if os.path.exists(nom_temp):
    os.remove(nom_temp)
