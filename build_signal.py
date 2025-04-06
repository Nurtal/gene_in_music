import pandas as pd
import numpy as np
import random
import os
from scipy.io.wavfile import write
from scipy.interpolate import interp1d
from mygene import MyGeneInfo


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





def ensembl_to_uniprot(ensembl_ids:list) -> dict:
    """Get Prot ID from genes ensembl_ids

    Args:
        - ensembl_ids (list) : list of ensembl id genes

    Returns:
        (dict) : gene ensembl id to protein id 
    
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
    
    gene_list = ['ENSG00000100982.12','ENSG00000226580.1','ENSG00000230978.2','ENSG00000160190.14','ENSG00000283051.1','ENSG00000223461.1','ENSG00000184113.10']
    machin = ensembl_to_uniprot(gene_list)
    print(machin)
