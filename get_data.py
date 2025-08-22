import kagglehub
import pandas as pd
import os



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
    


if __name__ == "__main__":

    get_data_from_kaggle()
    
