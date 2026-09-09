import pandas as pd
import matplotlib.pyplot as plt
from matplotlib.patches import Patch

# Load feature data
df = pd.read_csv("results/features.csv")

# Define gene groups
disease_genes = ["ATXN1", "ATXN3", "C9orf72", "FMR1"]

# Everything else in your dataset is treated as a housekeeping/control gene
control_genes = [
    "ACTB", "B2M", "HPRT1", "RPL13A",
    "SDHA", "TBP", "TUBB"
]

# Assign a color to each gene
colors = []

for gene in df["gene"]:
    if gene in disease_genes:
        colors.append("tomato")
    else:
        colors.append("steelblue")


# =====================================================
# GRAPH 1: GC CONTENT SPIKINESS
# =====================================================

plt.figure(figsize=(11, 6))

plt.bar(
    df["gene"],
    df["spikiness"],
    color=colors
)

plt.xlabel("Gene")
plt.ylabel("Spikiness")
plt.title("GC Content Spikiness Across Genes")

plt.xticks(rotation=45)
plt.legend(
    handles=[
        Patch(color="tomato", label="Disease-associated"),
        Patch(color="steelblue", label="Housekeeping/control")
    ]
)

plt.tight_layout()
plt.show()


# =====================================================
# GRAPH 2: MAXIMUM Z-SCORE
# =====================================================

plt.figure(figsize=(11, 6))

plt.bar(
    df["gene"],
    df["max_z_score"],
    color=colors
)

plt.xlabel("Gene")
plt.ylabel("Maximum Z-score")
plt.title("Maximum GC Content Z-score Across Genes")

plt.xticks(rotation=45)
plt.legend(
    handles=[
        Patch(color="tomato", label="Disease-associated"),
        Patch(color="steelblue", label="Housekeeping/control")
    ]
)

plt.tight_layout()
plt.show()