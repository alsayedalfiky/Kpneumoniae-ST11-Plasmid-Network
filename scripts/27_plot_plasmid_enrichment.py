#!/usr/bin/env python3
"""
Visualization of ST-Plasmid Evolutionary Partnerships
Generates a bubble plot (dot plot) showcasing statistically significant 
associations between Klebsiella Sequence Types and Plasmid Replicons.

Input: results/significant_plasmid_st_associations.csv
Visual Mapping:
    - X-axis: Plasmid Replicon Type
    - Y-axis: Sequence Type (ST)
    - Bubble Size: Prevalence (% of isolates in ST carrying the plasmid)
    - Bubble Color: Odds Ratio (Strength of statistical enrichment)
"""

import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns

# ---------------------------------------------------------
# 1. Load Enrichment Results
# ---------------------------------------------------------
# Loading the results generated in Step 26
df = pd.read_csv('results/significant_plasmid_st_associations.csv')

# ---------------------------------------------------------
# 2. Data Preparation and Normalization
# ---------------------------------------------------------
# Calculate Prevalence (%) to serve as the bubble size
# This represents the internal carriage rate within each lineage
df['Prevalence'] = (df['Observed'] / df['ST_Total']) * 100

# Handle 'inf' OddsRatio and apply capping for visualization
# Capping at 250 ensures that high-outlier associations (e.g., ST11/IncFII) 
# do not distort the color scale for other significant associations.
cap_val = 250
df['Odds_Capped'] = df['OddsRatio'].replace(np.inf, cap_val)
df['Odds_Capped'] = df['Odds_Capped'].clip(upper=cap_val)

# Sort the layout for improved readability (by ST then Prevalence)
df = df.sort_values(['ST', 'Prevalence'], ascending=[True, False])

# ---------------------------------------------------------
# 3. Main Plot Construction
# ---------------------------------------------------------
plt.figure(figsize=(14, 10))
sns.set_style("whitegrid")

# Create the scatter plot
# s = bubble size (Prevalence scaled by 8 for visibility)
# c = color (Capped Odds Ratio)
scatter = plt.scatter(
    x=df['Replicon'], 
    y=df['ST'].astype(str), 
    s=df['Prevalence'] * 8, 
    c=df['Odds_Capped'], 
    cmap='YlOrRd', 
    edgecolors='black', 
    linewidth=0.5,
    alpha=0.85
)

# ---------------------------------------------------------
# 4. Colorbar and Legend Configuration
# ---------------------------------------------------------
# Add Colorbar to represent the Odds Ratio (Enrichment Strength)
cbar = plt.colorbar(scatter)
cbar.set_label('Odds Ratio (Enrichment Strength)', rotation=270, labelpad=20, fontsize=12)

# Manually create a Size Legend to represent Prevalence percentages
for size in [25, 50, 75, 100]:
    plt.scatter([], [], c='gray', alpha=0.5, s=size*8, label=f'{size}%', edgecolors='black')

plt.legend(title="% Carriage (Prevalence)", labelspacing=1.5, borderpad=1.5,
           loc='upper left', bbox_to_anchor=(1.25, 1), frameon=True)

# ---------------------------------------------------------
# 5. Final Formatting and Export
# ---------------------------------------------------------
plt.title('Plasmid-ST Evolutionary Partnerships\nSize = % of ST Carrying Plasmid | Color = Odds Ratio', 
          fontsize=16, pad=20)
plt.xlabel('Plasmid Replicon Type', fontsize=12)
plt.ylabel('Klebsiella Sequence Type (ST)', fontsize=12)
plt.xticks(rotation=45, ha='right')

plt.tight_layout()

# Save at publication quality (300 DPI)
plt.savefig('results/plasmid_partnership_refined.png', dpi=300)
print("Figure successfully saved to: results/plasmid_partnership_refined.png")

plt.show()
