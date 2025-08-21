import pickle
import pandas as pd
import networkx as nx


def compute_graph_distance(graph_file_name:str, matrix_file_name:str):
    """Generate distance matrix of nodes within a stringdb graph

    Args:
        - graph_file_name (str) : path to pickle file to load the graph
        - matrix_file_name (str) : path to matrix file to save the distances
    
    """

    # load graph
    G = pickle.load(open(graph_file_name, 'rb'))

    # 1. Transformer les poids STRING en coût (plus le score est grand, plus le coût est petit)
    for u, v, data in G.edges(data=True):
        data['cost'] = 1.0 / data['weight']

    # 2. Calculer toutes les distances pondérées (weighted shortest path)
    distances = dict(nx.all_pairs_dijkstra_path_length(G, weight="cost"))

    # 3. Convertir en DataFrame (matrice de distances)
    nodes = list(G.nodes())
    dist_matrix = pd.DataFrame(index=nodes, columns=nodes, dtype=float)
    for i in nodes:
        for j in nodes:
            dist_matrix.loc[i, j] = distances[i].get(j, float("inf"))  # inf si pas de chemin

    # 4. Sauvegarde CSV
    dist_matrix.to_csv(matrix_file_name)


if __name__ == "__main__":

    compute_graph_distance("/tmp/string_network.pickle", "/tmp/dist.csv")

    
