import pandas as pd
import numpy as np
from sklearn.manifold import MDS


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


if __name__ == "__main__":

    get_proximity_from_data(['data/gene_reads_artery_aorta.csv', 'data/gene_reads_artery_coronary.csv'], "data/prox_matrix.csv")
    build_order_from_proximity("data/prox_matrix.csv")
    
