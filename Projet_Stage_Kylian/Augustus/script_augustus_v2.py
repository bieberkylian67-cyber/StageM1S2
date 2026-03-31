#Ce script permet de lancer augustus sur chaque séquence d'un fichier fasta 
import sys
import os 
from Bio import SeqIO # on importe dans le module biopython SeqIo pour pouvoir extraire plus facilement les séquences et ids
# Input
#Argument 1 = nom du fichier entrée qui contient les séquences au format fasta
#Argument 2 = nom du fichier de sortie qui contient les prédictions de augustus au format .gff3
nom_paquet = sys.argv[1]
nom_resultat = sys.argv[2]
# On crée un nom de fichier temporaire unique pour ce paquet
nom_temp = "/gstock/user/bieber/temporaire/temp_" + os.path.basename(nom_paquet) #fichier temporaire pour stocker

dictionnaire = {
    # 
    "ENST00": "human",    # Homo sapiens
    "ENSMUS": "human",    # Mus musculus --> pas de modèle souris dans Augustus
    "ENSRNO": "human",    # Rattus norvegicus --> pas de modèle rat dans Augustus
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

    #    oiseau
    "ENSGAL": "chicken",  # Gallus gallus
    "ENSTGU": "chicken",  # Taeniopygia guttata
    "ENSMGA": "chicken",  # Meleagris gallopavo
    "ENSAPL": "chicken",  # Anas platyrhynchos
    "ENSANA": "chicken",  # Anas sp.
    "ENSAME": "chicken",  # Dromaius novaehollandiae

    #    poisson 
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
    "ENSLAC": "elephant_shark",   # Latimeria chalumnae c'est le requin le plus proche
    "ENSCMI": "elephant_shark",   # Callorhinchus milii
    "ENSCIN": "ciona",            # Ciona intestinalis
    "ENSXMA": "Xiphophorus_maculatus",  # Xiphophorus maculatus
    "ENSPMA": "sealamprey",       # Petromyzon marinus
}

# Permet de lire le fichier d'entrée et de récupérer les séquences et les ids 
for sequences in SeqIO.parse(nom_paquet, "fasta"):
    id = sequences.id # .id permet de récupérer l'ID facilement 
    prefixe = id[0:6]
    species = dictionnaire.get(prefixe, "human") # prend le species qui correspond dans le dictionnaire , .get pour éviter une erreur de clé et que le script s'arrête et prend human par défaut 
    with open (nom_temp, "w") as temp:
        SeqIO.write(sequences, temp, "fasta") # 3 arguments : sequences = la séquence à écrire, temp = le fichier dans lequel écrire, "fasta" = le format 
    commande = f"augustus --species={species} --gff3=on {nom_temp} >> {nom_resultat}"
    os.system(commande)
        
#Efface le fichier temporaire 
if os.path.exists(nom_temp):
    os.remove(nom_temp)
