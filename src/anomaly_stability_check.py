import pandas as pd
from sklearn.ensemble import IsolationForest
from sklearn.preprocessing import StandardScaler

df = pd.read_csv("results/features.csv")

feature_cols = ["spikiness", "max_z_score"]
X = df[feature_cols]

scaler = StandardScaler()
X_scaled = scaler.fit_transform(X)

contamination = 4/11

results = {}
for seed in [0, 1, 21, 42, 99]:
    model = IsolationForest(contamination=contamination, random_state=seed)
    scores = model.fit_predict(X_scaled)
    results[f"seed_{seed}"] = scores

stability_df = df[["gene"]].copy()
for col, scores in results.items():
    stability_df[col] = scores

stability_df["flagged_count"] = (stability_df.filter(like="seed_") == -1).sum(axis=1)

stability_df.to_csv("results/anomaly_stability_check.csv", index=False)
print(stability_df)