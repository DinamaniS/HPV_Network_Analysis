import pandas as pd
import networkx as nx
import numpy as np
import os
from networkx.algorithms import community

# ======================================================
# ✅ Define paths
# ======================================================
base_dir = r"H:\Dinamani\Cervical_Cancer_Virus\Cervical_Cancer_Proteome_data\Network_New\Core_One_viral_Protein"
output_dir = os.path.join(base_dir, "Supplementary_Advanced_Topology")
os.makedirs(output_dir, exist_ok=True)

# ======================================================
# ✅ Load host–virus interaction data
# ======================================================
hpv16_df = pd.read_excel(os.path.join(base_dir, "HPV16_Host.xlsx"))
hpv18_df = pd.read_excel(os.path.join(base_dir, "HPV18_Host.xlsx"))

# ======================================================
# ✅ Build bipartite graphs
# ======================================================
def build_graph(df):
    G = nx.Graph()
    for _, row in df.iterrows():
        host = row["Host_Protein"]
        virus = row["Virus_Protein"]
        G.add_edge(host, virus)
    return G

G16 = build_graph(hpv16_df)
G18 = build_graph(hpv18_df)

# ======================================================
# ✅ Compute Advanced Metrics
# ======================================================
def compute_metrics(G, virus_name):
    metrics = {}

    # Degree-related metrics
    degree_centrality = nx.degree_centrality(G)
    betweenness_centrality = nx.betweenness_centrality(G)
    closeness_centrality = nx.closeness_centrality(G)
    try:
        eigenvector_centrality = nx.eigenvector_centrality(G, max_iter=1000)
    except nx.NetworkXError:
        eigenvector_centrality = {n: np.nan for n in G.nodes()}

    clustering_coeff = nx.clustering(G)
    
    # Community detection using modularity
    comm = community.greedy_modularity_communities(G)
    modularity_value = community.quality.modularity(G, comm)

    # Save node-level data
    df = pd.DataFrame({
        "Node": list(G.nodes()),
        "Degree": [G.degree(n) for n in G.nodes()],
        "Degree_Centrality": [degree_centrality[n] for n in G.nodes()],
        "Betweenness_Centrality": [betweenness_centrality[n] for n in G.nodes()],
        "Closeness_Centrality": [closeness_centrality[n] for n in G.nodes()],
        "Eigenvector_Centrality": [eigenvector_centrality[n] for n in G.nodes()],
        "Clustering_Coefficient": [clustering_coeff[n] for n in G.nodes()],
        "Community": [
            next((i+1 for i, c in enumerate(comm) if n in c), None) for n in G.nodes()
        ]
    })

    # Export
    df.to_excel(os.path.join(output_dir, f"{virus_name}_Advanced_Metrics.xlsx"), index=False)

    # Network-level summary
    metrics = {
        "Nodes": G.number_of_nodes(),
        "Edges": G.number_of_edges(),
        "Density": nx.density(G),
        "Average_Clustering": nx.average_clustering(G),
        "Number_of_Communities": len(comm),
        "Modularity": modularity_value,
        "Average_Degree": np.mean([G.degree(n) for n in G.nodes()])
    }

    return metrics

metrics_16 = compute_metrics(G16, "HPV16")
metrics_18 = compute_metrics(G18, "HPV18")

# ======================================================
# ✅ Save Summary
# ======================================================
summary = pd.DataFrame([
    {"Metric": k, "HPV16": v, "HPV18": metrics_18[k]} for k, v in metrics_16.items()
])
summary.to_excel(os.path.join(output_dir, "Supplementary_Advanced_Network_Summary.xlsx"), index=False)

print("\n✅ Supplementary topology analysis complete!")
print(f"Results saved to: {output_dir}")
