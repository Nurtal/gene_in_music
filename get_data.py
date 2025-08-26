import kagglehub
import pandas as pd
import os
import requests
import gzip
import shutil


def get_data_from_kaggle():
    """Download a real dataset directly from kaggle"""

    # Download latest version of the dataset
    path = kagglehub.dataset_download("albertozorzetto/rnaseq-aging-dementia-and-tbi")

    # reformat data
    df = pd.read_csv(f"{path}/fpkm_table_normalized.csv")
    gene_to_symbol = {}
    df_gene = pd.read_csv(f"{path}/rows-genes.csv")
    for index, row in df_gene.iterrows():
        gene = row['gene_id']
        symbol = row['gene_symbol']
        gene_to_symbol[gene] = symbol
    df["gene_id \ rnaseq_profile_id"] = df["gene_id \ rnaseq_profile_id"].replace(gene_to_symbol)
    df = df.rename(columns={"gene_id \ rnaseq_profile_id":"ID"})
    df = df.set_index('ID')
    df = df.T
    df['ID'] = df.index

    id_to_group = {}
    df_sample = pd.read_csv(f"{path}/columns-samples.csv")
    for index, row in df_sample.iterrows():
        id = row['rnaseq_profile_id']
        group = row['structure_acronym']
        id_to_group[str(id)] = group
    df['GROUP'] = df['ID'].replace(id_to_group)

    # reorder variables
    vars = ['ID']
    for v in list(df.keys()):
        if v != 'ID':
            vars.append(v)
    df = df[vars]

    # save data
    if not os.path.isdir("data"):
        os.mkdir("data")
    df.to_csv("data/kaggle_dementia.csv", index=False)



def get_data_from_stringdb(link_save_path:str, info_save_path:str) -> None:
    """Download ressources file for gene order computation from stringdb
    
    Args:
        - link_save_path (str) : path to link save file (usually a txt file, gz decompression is handle within the function)
        - info_save_path (str) : path to info save file (usually a txt file, gz decompression is handle within the function)
    
    """

    # params
    link_file_url = "https://stringdb-downloads.org/download/protein.links.v12.0/9606.protein.links.v12.0.txt.gz"
    info_file_url = "https://stringdb-downloads.org/download/protein.info.v12.0/9606.protein.info.v12.0.txt.gz"
    link_tmp_path = "/tmp/link_protein.gz"
    info_tmp_path = "/tmp/info_protein.gz"

    # download protein link file and extract gz content
    response = requests.get(link_file_url)
    with open(link_tmp_path, 'wb') as f:
        f.write(response.content)
    with gzip.open(link_tmp_path, "rb") as f_in:
        with open(link_save_path, "wb") as f_out:
            shutil.copyfileobj(f_in, f_out)
    os.remove(link_tmp_path)

    # download protein info file and extract gz content
    response = requests.get(info_file_url)
    with open(info_tmp_path, 'wb') as f:
        f.write(response.content)
    with gzip.open(info_tmp_path, "rb") as f_in:
        with open(info_save_path, "wb") as f_out:
            shutil.copyfileobj(f_in, f_out)
    os.remove(info_tmp_path)


if __name__ == "__main__":

    # get_data_from_kaggle()
    get_data_from_stringdb("/tmp/soubidou.txt", "/tmp/mashcidnefff.txt")
    
