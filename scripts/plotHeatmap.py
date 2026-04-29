import argparse
import numpy as np
import pandas as pd

import seaborn as sb
import matplotlib.pyplot as plt

from pathlib import Path
from scipy.spatial.distance import squareform
from scipy.cluster.hierarchy import leaves_list, linkage

from rates2chip.utilities import import_pandas, export_pandas

def cluster_by_similarity(corr: pd.DataFrame) -> pd.DataFrame:
	distance = (1.0 - corr.abs())

	# Convert the square distance matrix into SciPy's condensed vector format.
	# Because the matrix is symmetric, only the upper triangle is needed.
	condensed_distance = squareform(distance.to_numpy(), checks=False)
	
	# Run average-linkage hierarchical clustering on the condensed distances.
	# This builds the same type of merge tree used by dendrogram-based clustering.
	linkage_matrix = linkage(condensed_distance, method="average")

	# Extract the feature order implied by the dendrogram leaves.
	order_idx = leaves_list(linkage_matrix)

	# Map those integer positions back to the original feature names.
	order = corr.index[order_idx]

	# Reorder both rows and columns so similar features sit next to each other.
	return corr.loc[order, order]

if __name__ == "__main__":
	parser = argparse.ArgumentParser()
	parser.add_argument("--input", help="input table path")
	parser.add_argument("--features", nargs="+", default=[], help="specific feature columns to include")
	parser.add_argument("--exclude", nargs="+", default=[], help="feature columns to exclude")
	parser.add_argument("--corr_out", required=False, help="path for long-form pairwise correlation table")
	parser.add_argument("--heatmap_out", required=False, help="path for correlation heatmap image")
	parser.add_argument("--cluster",action='store_true')
	parser.add_argument("--figsize", nargs=2, type=float, default=[20.0, 16.0], help="heatmap figure size")
	args = parser.parse_args()

	df = import_pandas(Path(args.input))

	# Drop structural columns
	df.drop(columns=["gene", "chromosome", "strand", "start", "stop"],inplace=True)

	# Drop target column
	if args.target and args.target in df.columns:
		df = df.drop(columns=[args.target])

	# Gather features to plot
	if args.features:
		features = [c for c in args.features if c in df.columns]
	else:
		features = [c for c in df.columns if c not in args.exclude]

	# Calculate correlations
	x = df[features].corr(method="spearman")

	# Cluster
	if args.cluster:
		x = cluster_by_similarity(x)

	if args.corr_out:
		export_pandas(x, Path(args.corr_out))

	if args.heatmap_out:
		plt.figure(figsize=(args.figsize[0], args.figsize[1]))
		sb.heatmap(x, cmap="vlag", center=0.0, vmin=-1.0, vmax=1.0, square=True)
		plt.title("Feature-Feature Spearman Correlation")
		plt.xticks(rotation=45, ha='right')
		plt.tight_layout()
		plt.savefig(Path(args.heatmap_out), dpi=300)
		plt.close()
