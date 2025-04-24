import pandas as pd
from functools import reduce
import random
import mygene
import os
import shutil


def entrez_to_ensembl(entrez_gene_list:list) -> list:
    """Use a manually downloaded ressource file to convert a list of entre gene into a list
    of ensembl gene

    Args:
        - entrez_gene_list (list) : list of entrez gene

    Returns:
        - (list) : list of ensembl genes
    
    """
    
    # param
    data_file = "data/mart_export2.txt"
    df = pd.read_csv(data_file).dropna()
    df['NCBI gene (formerly Entrezgene) ID'] = df['NCBI gene (formerly Entrezgene) ID'].astype(int).astype(str)

    # get gene to scan
    to_scan = []
    for g in entrez_gene_list:
        if g in list(df['NCBI gene (formerly Entrezgene) ID']) and g not in to_scan:
            to_scan.append(g)
    df = df[df['NCBI gene (formerly Entrezgene) ID'].isin(to_scan)]

    # get ensembl genes
    ensembl_gene_list = list(df['Gene stable ID'])

    return ensembl_gene_list

    
            


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
        df['ID'] = df.index
        cols = df.columns.tolist()
        new_order = [cols[-1]] + cols[:-1]
        df = df[new_order]

        # save
        df.to_csv(gct_file.replace(".gct", ".csv"), index=False)



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
        df['ID'] = df.index
        cols = df.columns.tolist()
        new_order = [cols[-1]] + cols[:-1]
        df = df[new_order]

        # save
        df.to_csv(gct_file.replace(".gct", ".csv"), index=False)



def craft_gsea_dataset(gct_file_list:list, gmt_file:str, output_folder:str) -> None:
    """craft gsea datasets, one file for each gene set extracted from gmt_file / for each of the data files in gct_file_list
    
    Args:
        - gct_file_list (list) : list of path for gct files
        - gmt_file (str) : path to the gmt file containing the gene sets
        - output_folder (str) : path to the folder use to save generated datasets
    
    """

    # prepare result folder
    if os.path.exists(output_folder) and os.path.isdir(output_folder):
        shutil.rmtree(output_folder)
    os.mkdir(output_folder)

    # look for genes in data
    gene_list_list = []
    for gct_file in gct_file_list:
        df = read_gct(gct_file)
        gene_list_list.append(list(df['Name']))

    # compute intersection
    gene_intersection_list = list(reduce(lambda a, b: set(a) & set(b), gene_list_list))

    # get gene set to gene list
    gene_sets = {}
    with open(gmt_file, "r") as f:
        for line in f:
            parts = line.strip().split("\t")
            name = parts[0]
            genes = parts[2:]  # on ignore la description
            gene_sets[name] = genes

    # init mygene stuff
    mg = mygene.MyGeneInfo()
    
    # craft a dataset for each gct file
    for gct_file in gct_file_list:

        # load
        df = read_gct(gct_file)

        # reformat
        df = df[df['Name'].isin(gene_intersection_list)]
        df = df.drop(columns=['Description'])
        df = df.rename(columns={'Name':'ID'})
        df = df.set_index('ID')
        df = df.T
        df['ID'] = df.index
        cols = df.columns.tolist()
        new_order = [cols[-1]] + cols[:-1]
        df = df[new_order]
        df.columns = df.columns.str.replace(r"\.\d+$", "", regex=True)

        # split to gene set
        for gene_set in gene_sets:
            gene_list = gene_sets[gene_set]

            # convert entrez gene from gmt data to ensembl gene to match gct files
            print(gene_list)
            ensembl_gene_list = entrez_to_ensembl(gene_list) 
                        
            # select ensembl gene found in data
            ensembl_gene_list_to_keep = ['ID']
            for ensembl_gene in ensembl_gene_list:
                if ensembl_gene in list(df.keys()):
                    ensembl_gene_list_to_keep.append(ensembl_gene)

            # create subset
            if len(ensembl_gene_list_to_keep) > 1:
                data_file_name = gct_file.split("/")[-1].replace(".gct", f"_{gene_set}.csv")
                df_spe = df[ensembl_gene_list_to_keep]
                df_spe.to_csv(f"{output_folder}/{data_file_name}", index=False)


      

if __name__ == "__main__":


    # craft_reduce_datasets(["data/gene_reads_artery_aorta.gct", "data/gene_reads_artery_coronary.gct"], 5)
    # craft_datasets(["data/gene_reads_artery_aorta.gct", "data/gene_reads_artery_coronary.gct"])
    # craft_gsea_dataset(["data/gene_reads_artery_aorta.gct", "data/gene_reads_artery_coronary.gct"], "data/h.all.v2024.1.Hs.entrez.gmt", "/tmp/zog")

    entrez_to_ensembl(['AGAP12P-203', 'OR4C45-202'])
