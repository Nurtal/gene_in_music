import mygene

# Initialiser le client
mg = mygene.MyGeneInfo()

# Exemple de liste de gènes à mapper
genes = ["RNA5SP42", "C1orf94", "LOC100420336", "EGFR", "TP53"]

# Interroger mygene
results = mg.querymany(
    genes,
    scopes="symbol",           # on dit que nos IDs sont des symboles (HGNC/LOC/ORF/etc.)
    fields="symbol,name,entrezgene,ensembl.gene",
    species="human"
)

# Afficher résultats
for r in results:
    print(r)
