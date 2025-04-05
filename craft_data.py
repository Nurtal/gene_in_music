import pandas as pd
from functools import reduce
import random


def read_gct(filepath:str) -> pd.DataFrame:
    """Load content of a gct file into dataframe and return it

    Args:
        - filepath (str) : path to the gct file
    
    """

    # vibe coding
    with open(filepath, 'r') as f:
        version = f.readline().strip()
        dims = f.readline().strip().split()
        n_rows, n_cols = map(int, dims)

    # Now read the rest of the file as a regular table
    df = pd.read_csv(filepath, sep='\t', skiprows=2)

    return df


def pick_random_genes(gct_file_list:list, n_pick:int) -> list:
    """Pick a random set of genes among the genes shared bythe gct files in gct_file_list

    Args:
        - gct_file_list (list) : list of path for gct files
        - n_pick (int) : number of gene to randomly pick

    Returns:
        - (list) : list of picked genes
    
    """

    # look for genes in data
    gene_list_list = []
    for gct_file in gct_file_list:
        df = read_gct(gct_file)
        gene_list_list.append(list(df['Name']))

    # compute intersection
    gene_list = list(reduce(lambda a, b: set(a) & set(b), gene_list_list))
    
    # random pick
    gene_list = random.sample(gene_list, n_pick)

    return gene_list



if __name__ == "__main__":


    pick_random_genes(["data/gene_reads_artery_aorta.gct", "data/gene_reads_artery_coronary.gct"], 5)
