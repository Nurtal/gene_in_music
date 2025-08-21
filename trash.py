import requests
import networkx as nx
import matplotlib.pyplot as plt

# -----------------------
# Paramètres de base
# -----------------------
genes_of_interest = ["IFNA1", "IFNB1", "STAT1", "IRF7", "MTOR", "AKT1", "RPS6KB1"]  # ton toy set
species = 9606  # Homo sapiens
score_threshold = 700  # seuil de confiance STRING (700 = haut)

# -----------------------
# Construire l'URL de STRING API
# -----------------------
base_url = "https://string-db.org/api"
output_format = "tsv-no-header"
method = "network"

request_url = f"{base_url}/{output_format}/{method}"

# Paramètres
params = {
    "identifiers": "%0d".join(genes_of_interest),  # join avec retour chariot
    "species": species,
    "required_score": score_threshold,
    "network_type": "functional",
}

# -----------------------
# Requête à STRING
# -----------------------
response = requests.post(request_url, data=params)
lines = response.text.strip().split("\n")

# -----------------------
# Construire le graphe avec NetworkX
# -----------------------
G = nx.Graph()

for line in lines:
    items = line.split("\t")
    if len(items) >= 6:
        gene1 = items[2]
        gene2 = items[3]
        score = float(items[5])
        G.add_edge(gene1, gene2, weight=score)

# -----------------------
# Visualisation rapide
# -----------------------
plt.figure(figsize=(8, 6))
pos = nx.spring_layout(G, seed=42)
nx.draw(G, pos, with_labels=True, node_color='skyblue', edge_color='gray', node_size=1500, font_size=12)
plt.title("Toy STRING Network")
plt.savefig("/tmp/graph.png")
plt.close()

