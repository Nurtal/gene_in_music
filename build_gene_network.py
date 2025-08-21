import requests
import networkx as nx
import matplotlib.pyplot as plt
import pickle


def build_gene_network(gene_list:list, graph_image_file:str, graph_save_file:str, score_threshold:int):
    """Build a gene graph using gene in gene_list as nodes and stringDB as prior knowledge to build the edges

    Args:
        - gene_list (list) : list of genes (nodes) to use
        - graph_image_file (str) : path to save the image of the graph
        - graph_save_file (str) : path to the save file (must be a .pickle)
        - score_treshold (int) : confidence threshold to build an edge between nodes, 700 is considered high confidence
    
    """

    # set params
    species = 9606

    # Construire l'URL de STRING API
    base_url = "https://string-db.org/api"
    output_format = "tsv-no-header"
    method = "network"
    request_url = f"{base_url}/{output_format}/{method}"

    # set parameters
    params = {
        "identifiers": "%0d".join(gene_list),
        "species": species,
        "required_score": score_threshold,
        "network_type": "functional",
    }

    # Send request to STRING
    response = requests.post(request_url, data=params)
    lines = response.text.strip().split("\n")

    # Build graph
    G = nx.Graph()
    for line in lines:
        items = line.split("\t")
        if len(items) >= 6:
            gene1 = items[2]
            gene2 = items[3]
            score = float(items[5])
            G.add_edge(gene1, gene2, weight=score)

    # save graph data
    pickle.dump(G, open(graph_save_file, 'wb'))
    
    # Craft figure
    plt.figure(figsize=(8, 6))
    pos = nx.spring_layout(G, seed=42)
    nx.draw(G, pos, with_labels=True, node_color='skyblue', edge_color='gray', node_size=1500, font_size=12)
    plt.title("STRING Network")
    plt.savefig(graph_image_file)
    plt.close()


if __name__ == "__main__":

    build_gene_network(["IFNA1", "IFNB1", "STAT1", "IRF7", "MTOR", "AKT1", "RPS6KB1"], "/tmp/zog.png", "/tmp/string_network.pickle", 500)
