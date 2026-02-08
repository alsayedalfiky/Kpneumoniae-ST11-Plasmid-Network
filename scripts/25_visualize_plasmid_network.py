#!/usr/bin/env python3
"""
Plasmid Association Network Visualization
This script generates a plasmid-plasmid network based on a bipartite projection 
of ST-Replicon associations.

Source: physical_plasmid_metadata_fixed.tsv
Logic: Bipartite weighted projection onto plasmid nodes.
Output: High-resolution network plot (PNG).
"""

import os
import pandas as pd
import networkx as nx
import matplotlib.pyplot as plt

# ---------------------------------------------------------
# 1. Setup and Directory Management
# ---------------------------------------------------------
# Ensure output directory exists for the final figure
os.makedirs("results/networks", exist_ok=True)

# ---------------------------------------------------------
# 2. Data Loading
# ---------------------------------------------------------
# Load the merged metadata file containing both physical links and STs
df = pd.read_csv("results/merged_physical_plasmid_metadata_fixed.tsv", sep="\t")

# ---------------------------------------------------------
# 3. Bipartite Graph Construction
# ---------------------------------------------------------
# Create a cross-tabulation of STs vs Replicons to build the bipartite structure
st_plasmid_matrix = pd.crosstab(df['ST'], df['Replicon'])

G_bipartite = nx.Graph()
for st in st_plasmid_matrix.index:
    for plasmid in st_plasmid_matrix.columns:
        weight = st_plasmid_matrix.loc[st, plasmid]
        if weight > 0:
            # Connect the ST node to the Plasmid node
            G_bipartite.add_edge(f"ST_{st}", plasmid, weight=weight)

# ---------------------------------------------------------
# 4. Weighted Projection onto Plasmid Nodes
# ---------------------------------------------------------
# Extracts the relationships between plasmids based on shared STs
plasmid_nodes = [n for n in G_bipartite.nodes() if not str(n).startswith("ST_")]
G_plasmid = nx.bipartite.weighted_projected_graph(G_bipartite, plasmid_nodes)

# ---------------------------------------------------------
# 5. Visual Attribute Configuration
# ---------------------------------------------------------
# Calculate degree, node sizes, and edge widths based on connections and weights
degrees = dict(G_plasmid.degree())
node_sizes = [degrees[n] * 30 for n in G_plasmid.nodes()]
edge_widths = [d['weight'] / 5 for _, _, d in G_plasmid.edges(data=True)]

# Color mapping based on degree
node_colors = [degrees[n] for n in G_plasmid.nodes()]

# ---------------------------------------------------------
# 6. Layout and Plotting
# ---------------------------------------------------------
# Using Kamada-Kawai layout for aesthetic node distribution
pos = nx.kamada_kawai_layout(G_plasmid)

plt.figure(figsize=(12, 10))

# Draw the network elements
nodes = nx.draw_networkx_nodes(G_plasmid, pos,
                               node_size=node_sizes,
                               node_color=node_colors,
                               cmap=plt.cm.plasma,
                               edgecolors='black')

nx.draw_networkx_edges(G_plasmid, pos, width=edge_widths, edge_color='gray')

# ---------------------------------------------------------
# 7. Labeling Hubs (Degree > 39)
# ---------------------------------------------------------
# Label plasmids that are highly connected within the population
labels = {node: str(node).split("_")[0] for node, deg in degrees.items() if deg > 39}
nx.draw_networkx_labels(G_plasmid, pos, labels=labels, font_size=10, font_weight='bold')

# ---------------------------------------------------------
# 8. Figure Finalization
# ---------------------------------------------------------
plt.colorbar(nodes, label="Node degree (connections)")
plt.title("Plasmid Association Network in Klebsiella pneumoniae")
plt.axis('off')
plt.tight_layout()

# Save at publication quality
plt.savefig("results/networks/plasmid_network_refined.png", dpi=600)
print("Network visualization successfully saved to results/networks/plasmid_network_refined.png")

plt.show()
