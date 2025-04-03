import pandas as pd
import numpy as np
import random
import os


def build_random_signal(data_file:str):
    """Build signals from data file, assign random order and random interval to genes
    Create one signal file per patient

    Args:
        - data_file (str) : path to the data_file
    
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
    if not os.path.isdir("signals"):
        os.mkdir("signals")

    # loop over patients
    for i in id_to_y:

        # vars
        cmpt = 0
        y = id_to_y[i]

        # save signal in file
        output_file = open(f"signals/{i}_signal.csv", "w")
        output_file.write("x,y\n")
        for x in x_positions:
            output_file.write(f"{x},{y[cmpt]}\n")
            cmpt+=1
        output_file.close()
    


if __name__ == "__main__":

    build_random_signal("data/toy_data.csv")

    
