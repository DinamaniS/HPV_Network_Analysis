import pandas as pd
import networkx as nx
import numpy as np
import matplotlib.pyplot as plt
import os
from matplotlib import rcParams

# =============================
# ✅ STEP 0 Setup Output Directory
# =============================
output_dir = r"H:\Dinamani\Oncolytic_Viral_Hepatocellular_Carcinoma\Cervical_Cancer_Proteome_data\Network_New\NEW_Degree\New"
os.makedirs(output_dir, exist_ok=True)

# =============================
# ✅ STEP 1 Load Mapped Data
# =============================
hpv16_filtered = pd.read_excel(os.path.join(output_dir, "HPV16_Host.xlsx"))
hpv18_filtered = pd.read_excel(os.path.join(output_dir, "HPV18_Host.xlsx"))

# =============================
# ✅ STEP 2 Build Network Graphs
# =============================
def build_graph(df):
    G = nx.Graph()
    for _, row in df.iterrows():
        G.add_edge(row['Host_Protein'], row['Virus_Protein'])
    return G

G_16 = build_graph(hpv16_filtered)
G_18 = build_graph(hpv18_filtered)

# =============================
# ✅ STEP 3 Compute Topological Metrics
# =============================
def compute_heterogeneity(G):
    degrees = [deg for _, deg in G.degree()]
    return np.std(degrees) / np.mean(degrees)

def compute_centralization(G):
    centralities = list(nx.degree_centrality(G).values())
    return max(centralities) - np.mean(centralities)

hetero_16 = compute_heterogeneity(G_16)
central_16 = compute_centralization(G_16)

hetero_18 = compute_heterogeneity(G_18)
central_18 = compute_centralization(G_18)

# Save metrics
metrics_df = pd.DataFrame({
    "Metric": ["Heterogeneity", "Centralization"],
    "HPV16": [hetero_16, central_16],
    "HPV18": [hetero_18, central_18]
})
metrics_df.to_excel(os.path.join(output_dir, "Network_Topology_Metrics.xlsx"), index=False)

# =============================
# ✅ STEP 4 Generate Random Network Metrics
# =============================
def random_metrics(G):
    nodes = G.number_of_nodes()
    edges = G.number_of_edges()
    R = nx.gnm_random_graph(nodes, edges)
    return compute_heterogeneity(R), compute_centralization(R)

r_hetero_16, r_central_16 = random_metrics(G_16)
r_hetero_18, r_central_18 = random_metrics(G_18)

# Save random vs real metrics
random_df = pd.DataFrame({
    "Metric": ["Heterogeneity", "Centralization"],
    "HPV16_real": [hetero_16, central_16],
    "HPV16_random": [r_hetero_16, r_central_16],
    "HPV18_real": [hetero_18, central_18],
    "HPV18_random": [r_hetero_18, r_central_18]
})
random_df.to_excel(os.path.join(output_dir, "Random_vs_Real_Metrics.xlsx"), index=False)

# =============================
# ✅ STEP 5 Plot and Save Heterogeneity & Centralization Barplots
# =============================

# Set global font to Times New Roman
rcParams['font.family'] = 'Times New Roman'

labels = ['HPV16', 'HPV16_Random', 'HPV18', 'HPV18_Random']
hetero_values = [hetero_16, r_hetero_16, hetero_18, r_hetero_18]
central_values = [central_16, r_central_16, central_18, r_central_18]
colors = ['red', 'grey', 'green', 'grey']

def save_plot(values, ylabel, title, filename):
    plt.figure(figsize=(7, 6))
    plt.bar(labels, values, color=colors)
    plt.xticks(rotation=45, ha='right', rotation_mode='anchor',
               fontsize=18, fontweight='bold')
    plt.yticks(fontsize=14)
    plt.ylabel(ylabel, fontsize=24, fontweight='bold')
    plt.title(title, fontsize=24, fontweight='bold')
    plt.tight_layout()

    # Save in TIFF and JPEG at both 300 and 600 dpi
    for fmt in ['tiff', 'jpeg']:
        for dpi in [300, 600]:
            plt.savefig(os.path.join(output_dir, f"{filename}_{dpi}dpi.{fmt}"),
                        dpi=dpi, format=fmt)
    plt.close()

# --- Save Heterogeneity ---
save_plot(hetero_values, "Heterogeneity Score", "Network Heterogeneity", "Network_Heterogeneity")

# --- Save Centralization ---
save_plot(central_values, "Centralization Score", "Network Centralization", "Network_Centralization")

# =============================
# ✅ STEP 6 Top Viral Proteins by Degree
# =============================
def top_viral_proteins(df, G):
    virus_nodes = df['Virus_Protein'].unique().tolist()
    degrees = {v: G.degree(v) for v in virus_nodes if v in G.nodes()}
    sorted_degrees = sorted(degrees.items(), key=lambda x: x[1], reverse=True)
    return sorted_degrees[:3], sorted_degrees

top_virus_16, all_virus_16 = top_viral_proteins(hpv16_filtered, G_16)
top_virus_18, all_virus_18 = top_viral_proteins(hpv18_filtered, G_18)

# Save full viral protein degrees
df_16_degrees = pd.DataFrame(all_virus_16, columns=['Virus_Protein', 'Degree'])
df_18_degrees = pd.DataFrame(all_virus_18, columns=['Virus_Protein', 'Degree'])

df_16_degrees.to_excel(os.path.join(output_dir, "HPV16_Viral_Protein_Degrees.xlsx"), index=False)
df_18_degrees.to_excel(os.path.join(output_dir, "HPV18_Viral_Protein_Degrees.xlsx"), index=False)

print("Top HPV16 Viral Proteins:", top_virus_16)
print("Top HPV18 Viral Proteins:", top_virus_18)

# =============================
# ✅ STEP 7 Export Host Proteins for Metascape
# =============================
hpv16_genes = sorted(hpv16_filtered['Host_Protein'].dropna().unique().tolist())
hpv18_genes = sorted(hpv18_filtered['Host_Protein'].dropna().unique().tolist())

with open(os.path.join(output_dir, "HPV16_HostProteins.txt"), "w") as f:
    f.write("\n".join(hpv16_genes))

with open(os.path.join(output_dir, "HPV18_HostProteins.txt"), "w") as f:
    f.write("\n".join(hpv18_genes))

print(f"✅ All files have been saved to {output_dir}")
