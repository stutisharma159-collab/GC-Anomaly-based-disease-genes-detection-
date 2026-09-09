import pandas as pd
from sklearn.ensemble import IsolationForest

# Load features
df = pd.read_csv("results/features.csv")

# Targeted GC-spike features
features = [
    "spikiness",
    "max_z_score"
]

X = df[features]

# Isolation Forest
model = IsolationForest(
    n_estimators=200,
    contamination=0.20,
    random_state=42
)

df["anomaly"] = model.fit_predict(X)
df["anomaly_score"] = model.decision_function(X)

df["anomaly_label"] = df["anomaly"].map({
    1: "normal",
    -1: "anomaly"
})

# Save separately
output_file = "results/anomaly_results_targeted_n25.csv"
df.to_csv(output_file, index=False)

print("\nTargeted GC-Spike Anomaly Detection")
print("===================================")

print(
    df[
        ["gene", "anomaly_label", "anomaly_score",
         "spikiness", "max_z_score"]
    ]
    .sort_values("anomaly_score")
    .to_string(index=False)
)

print(f"\nSaved to: {output_file}")