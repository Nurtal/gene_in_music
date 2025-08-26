import stringdb
import pandas as pd


def split_dataset_into_pathways(data_file):
    """ """

    # load dataset
    df = pd.read_csv(data_file)
    genes = []
    for k in list(df.keys()):
        if k not in ['ID', 'GROUP', 'LABEL']:
            genes.append(k)

    string_ids = stringdb.get_string_ids(genes)
    enrichment_df = stringdb.get_enrichment(string_ids.queryItem)
    print(enrichment_df)


if __name__ == "__main__":

    split_dataset_into_pathways("data/kaggle_dementia.csv")
