import os
import pandas as pd
import umap.umap_ as umap
from sklearn.preprocessing import StandardScaler
from sklearn.decomposition import PCA


def get_pca(df:pd.DataFrame) -> pd.DataFrame:
    """Compute PCA from a pandas dataframe, return coordinates in the PCA space""" 

    # Extract
    ids = df['ID']
    labels = df['LABEL']
    X = df.drop(columns=['ID', 'LABEL'])

    # Normalize
    scaler = StandardScaler()
    X_scaled = scaler.fit_transform(X)

    # Compute PCA
    pca = PCA(n_components=2)
    X_pca = pca.fit_transform(X_scaled)

    # craft & return df
    df = pd.DataFrame({
        'ID': ids,
        'LABEL': labels,
        'PC1': X_pca[:, 0],
        'PC2': X_pca[:, 1]
    })
    return df


def get_umap(df:pd.DataFrame) -> pd.DataFrame:
    """Compute UMAP from a pandas dataframe, return coordinates in the UMAP space"""

    # Extract
    ids = df['ID']
    labels = df['LABEL']
    X = df.drop(columns=['ID', 'LABEL'])

    # Normalize
    scaler = StandardScaler()
    X_scaled = scaler.fit_transform(X)

    # Compute UMAP
    umap_model = umap.UMAP(n_components=2, random_state=42)
    X_umap = umap_model.fit_transform(X_scaled)

    # craft & return df
    umap_df = pd.DataFrame({
        'ID': ids,
        'LABEL': labels,
        'UMAP1': X_umap[:, 0],
        'UMAP2': X_umap[:, 1]
    })
    return umap_df

def run(data_file, output_dir):
    """ """

    # init output folder
    if not os.path.isdir(output_dir):
        os.mkdir(output_dir)
        os.mkdir(f"{output_dir}/data")
        
    # load data
    df = pd.read_csv(data_file)

    # save raw data
    df.to_csv(f"{output_dir}/data/raw.csv", index=False)

    # generate & save PCA
    df_pca = get_pca(df)
    df_pca.to_csv(f"{output_dir}/data/pca.csv", index=False)

    # save umap data
    df_umap = get_umap(df)
    df_umap.to_csv(f"{output_dir}/data/umap.csv", index=False)

    # save wav data

    







if __name__ == "__main__":

    run("data/small_rnaseq.csv", "/tmp/zooog")
