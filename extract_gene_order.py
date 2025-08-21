import pandas as pd
import numpy as np
from sklearn.manifold import MDS
import mygene
import craft_data
import random


def get_proximity_from_data(data_file_list:list, matrix_save_file:str) -> None:
    """Compute proximity beween genes as the absolute correlation of genes expression within the merged datasets

    Args:
        - data_file_list (list) : list of path of datasets to use
        - matrix_save_file (str) : path to save the matrix as a csv file
    
    """

    # load and merge data
    dfs = [pd.read_csv(f, index_col=0) for f in data_file_list]
    df = pd.concat(dfs, axis=0)

    # compute absolute correlation matrix
    corr_matrix = df.corr(method='pearson')  # corrélation entre colonnes
    abs_corr = corr_matrix.abs()
    abs_corr = abs_corr.fillna(0)
    np.fill_diagonal(abs_corr.values, 1.0)

    # save matrix
    abs_corr.to_csv(matrix_save_file)


def build_random_gene_order_from_data(data_file:str) -> dict:
    """Build random gene order from a data file (trash dev purpose)"""

    # extract gene in random order
    df = pd.read_csv(data_file).drop(columns=['ID', 'LABEL'])
    gene_list = list(df.keys())
    random.shuffle(gene_list)
    
    # craft gene to position
    gene_to_pos = {}
    position = 0
    for gene in gene_list:
        gene_to_pos[gene] = position
        position+=1

    return gene_to_pos



def build_order_from_proximity(prox_matrix_file:str) -> dict:
    """Build gene positions from a proximity matrix

    Args:
        - prox_matrix_file (str) : path to the proximity matrix file

    Returns:
        - (dict) : gene to computed position
    
    """

    # load data 
    similarity_matrix = pd.read_csv(prox_matrix_file).rename(columns={'Unnamed: 0':'ID'})
    similarity_matrix = similarity_matrix.set_index('ID')

    # Convert similarity to distance
    distance_matrix = 1 - similarity_matrix

    # Apply MDS
    mds = MDS(n_components=1, dissimilarity='precomputed', random_state=42)
    coords_1d = mds.fit_transform(distance_matrix).flatten()

    # Translate vers le haut pour que le minimum soit 0
    coords_1d = coords_1d - coords_1d.min()
    
    # craft gene to position
    gene_to_position = dict(zip(similarity_matrix.index, coords_1d))

    return gene_to_position


def reorder_cols_from_gmt(data_file, gmt_file, output_file):
    """Use a gmt file to reorder cols to, supposed to work, need to test it sometimes, whatever"""

    # load data & extract genes
    df = pd.read_csv(data_file)
    present_gene_list = []
    for v in list(df.keys()):
        if v not in ['ID', 'LABEL'] and v not in present_gene_list:
            present_gene_list.append(v)

    # get gene set to gene list
    gene_sets = {}
    with open(gmt_file, "r") as f:
        for line in f:
            parts = line.strip().split("\t")
            name = parts[0]
            genes = parts[2:]  # on ignore la description
            gene_sets[name] = genes

    
    # build gene order
    gene_order = []
    mg = mygene.MyGeneInfo()
    for gene_set in gene_sets:
        gene_list = gene_sets[gene_set]
        ensembl_gene_list = craft_data.entrez_to_ensembl(gene_list) 
        for g in ensembl_gene_list:
            if g in present_gene_list and g not in gene_order:
                gene_order.append(g)

    # add gene that has not been spoted as part of a pathway
    for g in present_gene_list:
        if g not in gene_order:
            gene_order.append(g)

    # reorder cols
    col_order = ['ID']
    for g in gene_order:
        col_order.append(g)
    col_order.append('LABEL')
    df = df[col_order]
    df.to_csv(output_file, index=False)

        


def extract_order_from_graph_distances(distance_matrix_file:str) -> dict:
    """Compute gene order from proximity matrix builded from STringDB graph

    Args:
        - distance_matrix_file (str) : path to the distance matric file

    Returns:
        - (dict) : gene to position
    
    """

    # load distance matrix
    dist_matrix = pd.read_csv(distance_matrix_file, index_col=0)

    # find a linear order for genes
    nodes = list(dist_matrix.index)
    placed = {}
    
    # 1. Trouver la paire la plus proche
    min_dist = float("inf")
    pair = None
    for i in nodes:
        for j in nodes:
            if i != j and dist_matrix.loc[i, j] < min_dist:
                min_dist = dist_matrix.loc[i, j]
                pair = (i, j)
    
    a, b = pair
    placed[a] = 0.0
    placed[b] = min_dist
    used = {a, b}
    
    # 2. Construire l’ordre
    current = b
    while len(used) < len(nodes):
        # trouver le plus proche voisin de "current" qui n'est pas encore utilisé
        dists = dist_matrix.loc[current].drop(labels=used)
        next_gene = dists.idxmin()
        next_dist = dists.min()
        
        placed[next_gene] = placed[current] + next_dist
        used.add(next_gene)
        current = next_gene
    
    # 3. Retourner l’ordre trié par coordonnée
    ordered_genes = sorted(placed.items(), key=lambda x: x[1])

    # convert to dict
    gene_to_pos = {}
    for elt in ordered_genes:
        gene_to_pos[elt[0]] = float(elt[1])

    return gene_to_pos
    




if __name__ == "__main__":

    # get_proximity_from_data(['data/gene_reads_artery_aorta.csv', 'data/gene_reads_artery_coronary.csv'], "data/prox_matrix.csv")
    # build_order_from_proximity("data/prox_matrix.csv")

    # reorder_cols_from_gmt("data/small_rnaseq.csv", "data/h.all.v2024.1.Hs.entrez.gmt", "/tmp/zog.csv")
    
    # build_random_gene_order_from_data("data/small_rnaseq.csv")

    extract_order_from_graph_distances("/tmp/dist.csv")
