import pandas as pd
import numpy as np
import random
import os
from scipy.io.wavfile import write
from scipy.interpolate import interp1d
from mygene import MyGeneInfo
import requests
from itertools import chain


def build_random_signal(data_file:str, output_folder:str):
    """Build signals from data file, assign random order and random interval to genes
    Create one signal file per patient

    Args:
        - data_file (str) : path to the data_file
        - output_folder (str) : path to the output folder
    
    """

    # params
    min_dist = 2
    max_dist = 5

    # laod data
    df = pd.read_csv(data_file)

    # define a random order of genes
    gene_list = list(df.keys())[1:]
    gene_list = np.random.permutation(gene_list)

    # define random interval between genes
    distances = [random.randint(min_dist, max_dist) for _ in range(len(gene_list) - 1)]

    # Build signal - x part
    x_positions = [0] 
    for d in distances:
        x_positions.append(x_positions[-1] + d)

    # Build signal - y part
    id_to_y = {}
    for index, row in df.iterrows():
        y = []
        for gene in gene_list:
            y.append(row[gene])
        id_to_y[row['ID']] = y

    # craft & save signal
    if not os.path.isdir(output_folder):
        os.mkdir(output_folder)

    # loop over patients
    for i in id_to_y:

        # vars
        cmpt = 0
        y = id_to_y[i]

        # save signal in file
        output_file = open(f"{output_folder}/{i}_signal.csv", "w")
        output_file.write("x,y\n")
        for x in x_positions:
            output_file.write(f"{x},{y[cmpt]}\n")
            cmpt+=1
        output_file.close()


def build_signal_from_computed_positions(data_file:str, output_folder:str, gene_to_pos:dict):
    """Build signal from pre-computed positions for each genes
    Create one signal file per patient

    Args:
        - data_file (str) : path to the data_file
        - output_folder (str) : path to the output folder
        - gene_to_pos (dict) : gene to position
    
    """

    # laod data
    df = pd.read_csv(data_file)
    gene_list = list(gene_to_pos.keys())
    
    # Build signal
    id_to_x = {}
    id_to_y = {}
    for index, row in df.iterrows():
        y = []
        x = []
        for gene in gene_list:
            y.append(row[gene])
            x.append(gene_to_pos[gene])
        id_to_y[row['ID']] = y
        id_to_x[row['ID']] = x

    # craft & save signal
    if not os.path.isdir(output_folder):
        os.mkdir(output_folder)

    # loop over patients
    for i in id_to_y:

        # vars
        cmpt = 0
        y = id_to_y[i]
        x_positions = id_to_x[i]

        # save signal in file
        output_file = open(f"{output_folder}/{i}_signal.csv", "w")
        output_file.write("x,y\n")
        for x in x_positions:
            output_file.write(f"{x},{y[cmpt]}\n")
            cmpt+=1
        output_file.close()



def ensembl_to_uniprot(ensembl_ids:list) -> dict:
    """Get Prot ID from genes ensembl_ids

    Args:
        - ensembl_ids (list) : list of ensembl id genes

    Returns:
        (dict) : gene ensembl id to protein id [!] values can be list
    
    """

    # drop part after dots in id
    ensembl_ids = [i.split('.')[0] for i in ensembl_ids]

    # run request
    mg = MyGeneInfo()
    results = mg.querymany(ensembl_ids, scopes='ensembl.gene', fields='uniprot', species='human')
    gene_to_uniprot = {}

    # extract results
    for res in results:
        if 'uniprot' in res:
            up = res['uniprot']
            if isinstance(up, dict):
                gene_to_uniprot[res['query']] = up.get('Swiss-Prot') or up.get('TrEMBL')
            else:
                gene_to_uniprot[res['query']] = up
    
    return gene_to_uniprot


def get_string_ids(uniprot_ids:list) -> dict:
    """Use uniprot ids to get string ids

    Args:
        - uniprot_ids (list) : list of uniprot_ids

    Returns:
        (dict) : uniprot ids to string id 
    
    """

    # params
    species=9606    
    identifiers = '\r'.join(uniprot_ids)
    params = {
        'identifiers': identifiers,
        'species': species,
        'caller_identity': 'gene_in_music'
    }

    # call API
    response = requests.post("https://string-db.org/api/json/get_string_ids", data=params)
    data = response.json()

    # extract data
    uniprot_to_string = {item['queryItem']: item['stringId'] for item in data}

    return uniprot_to_string


def get_string_interactions(string_ids:list):
    """Get interactions between strings ids

    Args:
        - string_ids (list) : list of string ids

    Returns:
        - (json) : list of dict containing interaction results
    
    """

    # params
    species=9606
    identifiers = '\r'.join(string_ids)
    params = {
        'identifiers': identifiers,
        'species': species,
        'caller_identity': 'gene_in_music'
    }
    
    # call API
    response = requests.post("https://string-db.org/api/json/network", data=params)
    
    # return results
    return response.json()


def build_similarity_matrix(interactions) -> pd.DataFrame:
    """Build a similarity matrix between the identified proteins from interaction extracted
    with strings API

    Args:
        - interaction (list) : list of dict, kind of json, obtained from strings API

    Returns:
        - (pd.DataFrame) : proximity matrix
    
    """

    # forge matrix
    proteins = set()
    for interaction in interactions:
        proteins.add(interaction['preferredName_A'])
        proteins.add(interaction['preferredName_B'])
    proteins = list(proteins)
    matrix = pd.DataFrame(index=proteins, columns=proteins, data=0.0)

    # fill matrix
    for interaction in interactions:
        p1 = interaction['preferredName_A']
        p2 = interaction['preferredName_B']
        score = interaction['score']
        matrix.at[p1, p2] = score
        matrix.at[p2, p1] = score

    # fill diagonale
    for p in proteins:
        matrix.at[p, p] = 1.0

    return matrix



def build_gene_similarity_matrix(genes, gene_to_proteins, protein_similarity_matrix, output_file):
    """ [!] Vibe coding sucks
    Construit une matrice de similarité entre gènes à partir d'une matrice de similarité entre protéines,
    en prenant le score maximal entre toutes les paires de protéines associées aux deux gènes.

    Params:
    - genes : liste de gènes Ensembl (sans version)
    - gene_to_proteins : dict {gene: [prot1, prot2, ...]}  (même une seule prot = liste d’un élément)
    - protein_similarity_matrix : DataFrame (protéines en lignes & colonnes)

    Return:
    - gene_similarity_matrix : DataFrame (gènes en lignes & colonnes) avec similarité max entre protéines
    """
    
    # drop part after dots in id
    genes = [i.split('.')[0] for i in genes]

    # compute matrix
    gene_matrix = pd.DataFrame(index=genes, columns=genes, data=0.0)
    for gene1 in genes:
        for gene2 in genes:
            prots1 = gene_to_proteins.get(gene1, [])
            prots2 = gene_to_proteins.get(gene2, [])
            max_score = 0.0

            for p1 in prots1:
                for p2 in prots2:
                    if p1 in protein_similarity_matrix.index and p2 in protein_similarity_matrix.columns:
                        score = protein_similarity_matrix.at[p1, p2]
                        max_score = max(max_score, score)

            gene_matrix.at[gene1, gene2] = max_score

    # save matrix
    gene_matrix.to_csv(output_file)



def turn_signal_into_audio(signal_file:str, target_duration:float) -> None:
    """Turn a signal extracted from data file to an audio signal and save it in
    a wave file
    
    Args:
        - signal_file (str) : path to a signal file
        - target_duration (float) : duration of the audio fole
    
    """

    # load signal
    df = pd.read_csv(signal_file)
    x = np.array(list(df['x']))
    y = np.array(list(df['y']))

    # Redimensionner x pour fit la target_duration
    x = (x - x[0]) / (x[-1] - x[0]) * target_duration

    # audio stuff
    sample_rate = 44100  # 44.1 kHz
    duration = x[-1]
    t = np.linspace(0, duration, int(sample_rate * duration))

    # Interpolation pour obtenir des y réguliers
    interpolateur = interp1d(x, y, kind='linear')
    y_interp = interpolateur(t)

    # Normaliser pour correspondre à une plage 16-bit
    y_norm = np.int16(y_interp / np.max(np.abs(y_interp)) * 32767)

    # Sauvegarder dans un fichier WAV
    write(signal_file.replace(".csv", ".wav"), sample_rate, y_norm)


if __name__ == "__main__":

    # build_random_signal("data/toy_data.csv")
    
    # turn_signal_into_audio("signals/42_signal.csv", 4.0)

    # get gene list
    df = pd.read_csv("data/gene_reads_artery_aorta.csv")
    gene_list = list(df.keys())[5:500]
    
    # get gene to uniprot
    gene_to_uniprot = ensembl_to_uniprot(gene_list)
    gene_to_uniprot_format = {}
    for gene in gene_to_uniprot:
        uniprot = gene_to_uniprot[gene]
        if not isinstance(uniprot, list):
            gene_to_uniprot_format[gene] = [uniprot]
        else:
            gene_to_uniprot_format[gene] = uniprot

    # get proteine to proximity matrix
    uniprot_flat = list(chain.from_iterable(list(gene_to_uniprot.values())))
    uniprot_to_strings = get_string_ids(uniprot_flat)
    interactions = get_string_interactions(list(uniprot_to_strings.values()))
    uniprot_proximity_matrix = build_similarity_matrix(interactions)

    # compute gene proximity matrix
    build_gene_similarity_matrix(gene_list, gene_to_uniprot_format, uniprot_proximity_matrix, "data/test.csv")
