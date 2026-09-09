import pandas as pd
from sklearn.ensemble import IsolationForest
from sklearn.preprocessing import StandardScaler

df = pd.read_csv("results/features.csv")

feature_cols = ["spikiness", "max_z_score"]
X = df[feature_cols]

scaler = StandardScaler()
X_scaled = scaler.fit_transform(X)

model = IsolationForest(contamination=4, random_state=42)
df["anomaly_score_targeted"] = model.fit_predict(X_scaled)

df.to_csv("results/anomaly_results_targeted.csv", index=False)
print(df[["gene", "spikiness", "max_z_score", "anomaly_score_targeted"]].sort_values("max_z_score", ascending=False))