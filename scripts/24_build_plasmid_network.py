"""
Lineage-Associated Plasmid Network Construction
Constructs a network where nodes are plasmid replicons and edges represent 
their co-presence within the same Sequence Type (ST).

Logic:
- Grouping: By Sequence Type (ST).
- Edges: Created between all unique replicons found in the same ST.
- Weight: The number of STs where both replicons are present.
- Interpretation: High-weight edges indicate plasmids that share the same host lineages.
"""

import pandas as pd
import networkx as nx
from itertools import combinations
from collections import Counter
import os

# -----------------------------
# 1. Configuration
# -----------------------------
PHYSICAL_LINKS = "results/physical_plasmid_links.tsv"     # From Step 23 (Filtered ResFinder)
METADATA_FILE = "results/plasmid_with_metadata_fixed.tsv" # From Step 19 (Fixed STs)

OUTPUT_EDGES = "results/plasmid_network_edges_fixed.tsv"
OUTPUT_DEGREES = "results/plasmid_network_degrees_fixed.tsv"
OUTPUT_GRAPHML = "results/plasmid_association_fixed.graphml"

def main():
    if not os.path.exists(PHYSICAL_LINKS) or not os.path.exists(METADATA_FILE):
        print("Error: Required input files not found.")
        return

    print("Loading data for Lineage-Associated Network...")

    # -----------------------------
    # 2. Load and Filter Data
    # -----------------------------
    # Load list of plasmid-borne genes (to define relevant genomes)
    hits = pd.read_csv(PHYSICAL_LINKS, sep="\t", header=None)
    hits["Genome"] = hits[0]
    
    # Load metadata
    meta = pd.read_csv(METADATA_FILE, sep="\t")

    # Merge to link Plasmids to STs, keeping only rows relevant to our physical links
    # Note: We group by ST, so we need the unique set of replicons per ST
    joined = pd.merge(hits[["Genome"]], meta[["Genome", "Replicon", "ST"]], on="Genome", how="inner")

    # -----------------------------
    # 3. Build Network Logic (Population Level)
    # -----------------------------
    print(f"Grouping replicons by Sequence Type (ST)...")
    
    # Create the "Plasmid Pool" for each ST
    # This creates a list of all unique replicons found in ST11, ST258, etc.
    st_to_reps = joined.groupby("ST")["Replicon"].apply(lambda s: sorted(set(s))).to_dict()

    # Create Edges
    edges = []
    for st, reps in st_to_reps.items():
        # Connect every plasmid in this ST to every other plasmid in this ST
        if len(reps) > 1:
            for i in range(len(reps)):
                for j in range(i+1, len(reps)):
                    edges.append((reps[i], reps[j]))

    # -----------------------------
    # 4. Calculate Weights and Metrics
    # -----------------------------
    # Weight = Number of STs that contain this pair of plasmids
    edge_counter = Counter(tuple(sorted(edge)) for edge in edges)
    
    edge_df = pd.DataFrame(
        [(a, b, w) for (a, b), w in edge_counter.items()],
        columns=["Replicon_A", "Replicon_B", "Weight"]
    )

    # Build Graph for Degree Calculation
    G = nx.Graph()
    for _, row in edge_df.iterrows():
        G.add_edge(row["Replicon_A"], row["Replicon_B"], weight=row["Weight"])

    # Calculate Degrees (Connectivity in the population)
    degree_dict = dict(G.degree())
    degree_df = pd.DataFrame(list(degree_dict.items()), columns=["Replicon", "Degree"])

    # -----------------------------
    # 5. Save Outputs
    # -----------------------------
    edge_df.to_csv(OUTPUT_EDGES, sep="\t", index=False)
    degree_df.to_csv(OUTPUT_DEGREES, sep="\t", index=False)
    nx.write_graphml(G, OUTPUT_GRAPHML)

    print("-" * 30)
    print(f"Network Construction Complete.")
    print(f"Nodes (Replicons): {G.number_of_nodes()}")
    print(f"Edges (Co-associations): {G.number_of_edges()}")
    print(f"GraphML saved to: {OUTPUT_GRAPHML}")
    print("-" * 30)

if __name__ == "__main__":
    main()
