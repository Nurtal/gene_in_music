import pandas as pd
import numpy as np
from sklearn.manifold import MDS
import mygene
import craft_data
import random
import time

from mygene import MyGeneInfo



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
    



def extract_order_from_protein_distances(data_file:str, protein_link_file:str, protein_info_file:str, log_file:str, position_file:str):
    """Extract gene order using proximiy between associated proteins, from a local data file.

    Args:
        - data_file (str) : path to input data file, should contain gene in 'preferred_name' format (e.g EGFR, IFN ...)
        - protein_link_file (str) : path to ressource file downloaded from stringdb
        - protein_info_file (str) : path to ressource file downloaded from stringdb, used for ids converstion
        - log_file (str) : path to generated log file, contain list of ids that failed conversion
        - position_file (str) : path to generated position file, where results are saved 
       
    """

    # load data and extract genes to consider
    df = pd.read_csv(data_file)
    genes = []
    for k in list(df.keys()):
        if k not in ['ID', 'GROUP', 'LABEL']:
            genes.append(k)

    # get gene to gene id
    id_to_gene = {}
    df = pd.read_csv(protein_info_file, sep="\t")
    df = df[df['preferred_name'].isin(genes)]
    for index, row in df.iterrows():
        id_to_gene[row['#string_protein_id']] = row['preferred_name']

    # check missing genes
    log_data = open(log_file, "w")
    log_data.write("MISSING\n")
    for g in genes:
        if g not in list(id_to_gene.values()):
            log_data.write(f"{g}\n")
    log_data.close()

    # load and filter ressource file on genes
    df = pd.read_csv(protein_link_file, sep=" ")
    df = df[df['protein1'].isin(list(id_to_gene.keys()))]
    df = df[df['protein2'].isin(list(id_to_gene.keys()))]

    # compute gene order
    id_to_pos = {}
    pos = 0
    
    # init - find the closes entry in the dataset
    row_max = df.loc[df["combined_score"].idxmax()]
    id_to_pos[row_max['protein1']] = pos
    pos += float(1/row_max['combined_score'])
    id_to_pos[row_max['protein2']] = pos
    pivot = row_max['protein2']
    root = row_max['protein1']
    df = df[df['protein1']!=root]
    df = df[df['protein2']!=root]

    # its been a while since i used one of these
    iteration = 0
    start = time.time()
    while len(list(id_to_pos.keys())) < len(list(id_to_gene.keys())):

        # create sub df focus on pivot
        df_sub = df[df['protein1']==pivot]

        # check that there is something in the sub df
        if df_sub.shape[0] > 0:

            # find closest entry
            row_max = df_sub.loc[df_sub["combined_score"].idxmax()]

            # assign pos to closes entry
            pos += float(1/row_max['combined_score'])
            id_to_pos[row_max['protein2']] = pos

            # update, closes entry become next pivot, remove old pivot from df
            pivot = row_max['protein2']
            root = row_max['protein1']
            df = df[df['protein1']!=root]
            df = df[df['protein2']!=root]

            # display pregress
            iteration +=1
            progress = (len(list(id_to_pos.keys())) / len(list(id_to_gene.keys()))) * 100.0
            current_time = time.time()
            duration = current_time - start
            print(f"[ORDERING GENES][{iteration}][DURATION:{duration}] => {progress} ({len(list(id_to_pos.keys()))} / {len(list(id_to_gene.keys()))})")
        else:
            print("[!]PREMATURE STOP")
            break
        
    # translate ids
    id_to_pos_clean = {}
    for i in id_to_pos:
        id_to_pos_clean[id_to_gene[i]] = id_to_pos[i]
    id_to_pos = id_to_pos_clean

    # save in a result file
    result_data = open(position_file, "w")
    result_data.write("GENE,POS\n")
    for i in id_to_pos:
        result_data.write(f"{i},{id_to_pos[i]}\n")
    result_data.close()
            

def get_ensembl_genes(genes:list) -> dict:
    """Extract a dictionarry gene to ensembl ids, where input genes are suppused to be gene symbol (e.g EGFR)

    Args:
        - genes (list) : list of genes symbol

    Returns:
        - (dict) : keys are symbol and values are associated list of ensembl ids
    
    """

    # init mygene object
    mg = MyGeneInfo()
    results = mg.querymany(
        genes,
        scopes="symbol",     # type d'input
        fields="ensembl.gene",  # ce qu'on veut récupérer
        species="human"      # Homo sapiens
    )

    # Convertir en DataFrame
    df = pd.DataFrame(results)
    df_simple = df[["query", "ensembl"]]

    # Normaliser le champ "ensembl" (car parfois c'est une liste de dicts)
    def extract_ensg(x):
        if isinstance(x, list):
            return [d["gene"] for d in x]  # plusieurs ENSG possibles
        elif isinstance(x, dict):
            return x.get("gene")
        return None

    # compute gene to ensembl
    gene_to_ensembl = {}
    df_simple["ENSG"] = df_simple["ensembl"].apply(extract_ensg)
    for index, row in df_simple.iterrows():

        k = row['query']
        ensembl = row['ENSG']

        if isinstance(ensembl, list):
            gene_to_ensembl[k] = ensembl
        elif ensembl:
            gene_to_ensembl[k] = [ensembl]

    return gene_to_ensembl


def compute_gene_to_gene_distances(data_file:str, link_file:str, alias_file:str, output_file:str) -> None:
    """Use protein links and alias to craft two gene distance file (all distances and filtered distances, keep only the closest distance between symbols)

    Args:
        - data_file (str) : path to input data file, supposed to contains gene symbols as column, a first id column and a last label / group column
        - link_file (str) : path to stringdb links file
        - alias_file (str) : path to stringdb alias file
        - output_file (str) : path to save file (generate also a _filtered file containing only the best symbols distances) 
    
    """

    # load gene list
    df = pd.read_csv(data_file)
    genes = list(df.keys())[1:-1]

    # convert to ENSG id
    symbol_to_id = get_ensembl_genes(genes)
    id_to_symbol = {}
    ensg_list = []
    for symbol in symbol_to_id:
        sub_id_list = symbol_to_id[symbol]
        for i in sub_id_list:
            if i not in ensg_list:
                ensg_list.append(i)
                id_to_symbol[i] = symbol

    # Mapp ENSG to ENSP
    aliases = pd.read_csv(alias_file, sep="\t", skiprows=1, names=["protein", "alias", "source"])
    links = pd.read_csv(link_file, sep=" ")
    ensg_map = aliases[aliases["alias"].isin(ensg_list)][["alias", "protein"]]

    # extract score for all ENSP
    proteins = ensg_map["protein"].unique().tolist()
    sub_links = links[
        (links["protein1"].isin(proteins)) & 
        (links["protein2"].isin(proteins))
    ].copy()

    # add ENSG
    prot_to_ensg = dict(zip(ensg_map["protein"], ensg_map["alias"]))
    sub_links["gene1"] = sub_links["protein1"].map(prot_to_ensg)
    sub_links["gene2"] = sub_links["protein2"].map(prot_to_ensg)

    # Agréger au max par paire de gènes
    gene_links = (
        sub_links.groupby(["gene1", "gene2"])["combined_score"]
        .max()
        .reset_index()
    )

    # Convertir en distances
    gene_links["symbol1"] = gene_links["gene1"].replace(id_to_symbol)
    gene_links["symbol2"] = gene_links["gene2"].replace(id_to_symbol)
    gene_links["distance"] = 1 - (gene_links["combined_score"] / 1000)

    # keep only best values
    gene_links_filtered = gene_links.loc[gene_links.groupby(["symbol1", "symbol2"])["distance"].idxmin()].reset_index(drop=True)

    # save
    gene_links.to_csv(output_file, index=False)
    gene_links_filtered.to_csv(output_file.replace(".csv", "_filtered.csv"), index=False)


    
def extract_order_from_gene_distances(data_file:str, gene_distance_file:str, position_file:str, log_file:str) -> dict:
    """Extract gene order using proximiy between associated proteins, from a local data file.

    Args:
        - data_file (str) : path to input data file, should contain gene in 'preferred_name' format (e.g EGFR, IFN ...)
        - gene_distance_file (str) : file containing computed distances between genes
        - position_file (str) : path to generated position file, where results are saved 
        - log_file (str) : path to log file, store unavailable genes

    Returns:
        - (dict) : symbol to position
       
    """

    # load data and extract genes to consider
    df = pd.read_csv(data_file)
    genes = list(df.keys())[1:-1]

    # check genes available in distance file
    df = pd.read_csv(gene_distance_file)
    available_genes = []
    for x in list(df['symbol1']):
        if x not in available_genes:
            available_genes.append(x)
    for x in list(df['symbol2']):
        if x not in available_genes:
            available_genes.append(x)

    # check available genes
    log_data = open(log_file, "w")
    for g in genes:
        if g not in available_genes:
            log_data.write(f"MISSING DISTANCE INFORMATION FOR GENE {g}\n")
    log_data.close()

    # Trouver la paire avec la distance minimale
    idxmin = df["distance"].idxmin()
    root = df.loc[idxmin, "symbol1"]
    second = df.loc[idxmin, "symbol2"]
    d_root_second = df.loc[idxmin, "distance"]

    # Positions initiales
    positions = {root: 0, second: d_root_second}
    visited = {root, second}
    current = second

    # Boucle jusqu'à visiter tous les gènes
    while len(visited) < len(set(df["symbol1"]) | set(df["symbol2"])):
        # Trouver la plus petite distance depuis le "current" vers un gène non visité
        candidates = df[
            ((df["symbol1"] == current) & (~df["symbol2"].isin(visited)))
            | ((df["symbol2"] == current) & (~df["symbol1"].isin(visited)))
        ]

        if candidates.empty:
            break  # bloqué (graphe non connexe)

        next_idx = candidates["distance"].idxmin()
        row = candidates.loc[next_idx]

        if row["symbol1"] == current:
            nxt = row["symbol2"]
        else:
            nxt = row["symbol1"]

        # Position = position du current + distance
        positions[nxt] = positions[current] + row["distance"]

        # Avancer
        visited.add(nxt)
        current = nxt

    # save data
    data = []
    for symbol in positions:
        data.append({"GENE":symbol, 'POSITION':positions[symbol]})
    df = pd.DataFrame.from_dict(data)
    df.to_csv(position_file, index=False)

    return positions
        

    
    



if __name__ == "__main__":

    # get_proximity_from_data(['data/gene_reads_artery_aorta.csv', 'data/gene_reads_artery_coronary.csv'], "data/prox_matrix.csv")
    # build_order_from_proximity("data/prox_matrix.csv")

    # reorder_cols_from_gmt("data/small_rnaseq.csv", "data/h.all.v2024.1.Hs.entrez.gmt", "/tmp/zog.csv")
    
    # build_random_gene_order_from_data("data/small_rnaseq.csv")

    # extract_order_from_graph_distances("/tmp/dist.csv")
    
    # extract_order_from_protein_distances("data/fake_gene_data.csv", "data/9606.protein.links.v12.0.txt", "data/9606.protein.info.v12.0.txt", "/tmp/zog.log", "data/computed_positions.csv")

    # compute_gene_to_gene_distances("data/fake_gene_data.csv", "data/9606.protein.links.v12.0.txt", "data/stringdb_alias.txt", "data/gene_distances.csv")
    m = extract_order_from_gene_distances("data/fake_gene_data.csv", "data/gene_distances.csv", "/tmp/pos.csv", "/tmp/log.txt")
