import pandas as pd


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



if __name__ == "__main__":


    df = read_gct("data/gene_reads_artery_aorta.gct")

    print(df)
