import os


def get_pca(df):
    """ """ 

def run(data_file, output_dir):
    """ """

    # init output folder
    if not os.path.isdir(output_dir):
        os.mkdir(output_dir)
        os.mkdir(f"{output_dir}/data")
        
    # load data
    df = pd.read_csv(data_file)

    # save raw data
    df.to_csv(f"{output_dir}/data/raw.csv")

    # save PCA
    get_pca(df)

    # save umap daya

    # save wav data

    

    print("Tardis")






if __name__ == "__main__":

    run("data/small_rnaseq.csv", "/tmp/zooog")
