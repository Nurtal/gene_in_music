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


def craft_reduce_datasets(gct_file_list:list, n_genes:int) -> None:
    """Create a parsable csv file for each gct file gct_file_list, use only a random selection
    of genes to craft dataset

    Args:
        - gct_file_list (list) : list of path for gct files
        - n_pick (int) : number of gene to randomly pick
    """

    # pick a random list of genes
    gene_list = pick_random_genes(gct_file_list, n_genes)

    # craft a dataset for each gct file
    for gct_file in gct_file_list:

        # load
        df = read_gct(gct_file)

        # reformat
        df = df[df['Name'].isin(gene_list)]
        df = df.drop(columns=['Description'])
        df = df.rename(columns={'Name':'ID'})
        df = df.set_index('ID')
        df = df.T

        # save
        df.to_csv(gct_file.replace(".gct", ".csv"))



def craft_datasets(gct_file_list:list) -> None:
    """Create a parsable csv file for each gct file gct_file_list, use all possible genes

    Args:
        - gct_file_list (list) : list of path for gct files
    """

    # look for genes in data
    gene_list_list = []
    for gct_file in gct_file_list:
        df = read_gct(gct_file)
        gene_list_list.append(list(df['Name']))

    # compute intersection
    gene_list = list(reduce(lambda a, b: set(a) & set(b), gene_list_list))

    # craft a dataset for each gct file
    for gct_file in gct_file_list:

        # load
        df = read_gct(gct_file)

        # reformat
        df = df[df['Name'].isin(gene_list)]
        df = df.drop(columns=['Description'])
        df = df.rename(columns={'Name':'ID'})
        df = df.set_index('ID')
        df = df.T

        # save
        df.to_csv(gct_file.replace(".gct", ".csv"))

    

if __name__ == "__main__":


    # craft_reduce_datasets(["data/gene_reads_artery_aorta.gct", "data/gene_reads_artery_coronary.gct"], 5)
    craft_datasets(["data/gene_reads_artery_aorta.gct", "data/gene_reads_artery_coronary.gct"])
