

import sys
import subprocess
from Bio import SeqIO  # importe le module SeqIO de Biopython pour lire les fichiers FASTA

fichier_input = sys.argv[1]
fichier_sortie = sys.argv[2]

dictionnaire = {
    #    Humain 
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
    # Rongeurs 
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
    # Reptiles/Amphibiens
    "ENSXET": "human",    # Xenopus tropicalis
    "ENSACA": "human",    # Anolis carolinensis
    "ENSPCA": "human",    # Pelodiscus sinensis
    "ENSCPR": "human",    # Crocodylus porosus

    #    Chicken
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

with open(fichier_sortie, "w") as f:
    for sequences in SeqIO.parse(fichier_input, "fasta"):
        id = sequences.id
        prefixe = id[0:6]
        species = dictionnaire.get(prefixe, "human")
        with open("test_tmp.fa", "w") as tmp:
            tmp.write(sequences.format("fasta"))

        commande = ["augustus", f"--species={species}", "test_tmp.fa"]
        resultat = subprocess.run(commande, text=True, capture_output=True)
        resultat_texte = resultat.stdout.replace("test_tmp", sequences.id)
        f.write(resultat_texte)