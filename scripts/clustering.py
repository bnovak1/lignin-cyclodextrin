"""
Run HDBSCAN clustering on collective variables data for a given value of min_samples and save the model and number of clusters to files.
"""

import argparse
from pathlib import Path
import sys

import joblib
import numpy as np
from sklearn.cluster import HDBSCAN
from sklearn.pipeline import Pipeline
from sklearn.preprocessing import StandardScaler

sys.path.append(Path("../scripts"))
from read_colvar_file import read_colvar

# Command line inputs
parser = argparse.ArgumentParser()
parser.add_argument("--colvar", type=str, help="Input file with collective variables data")
parser.add_argument("--min_samples", type=int, help="Minimum number of samples")
parser.add_argument("--model_file", type=str, default=None,
                    help="File name to save the clustering model to")
parser.add_argument("--nclusters_file", type=str, default=None,
                    help="File name to save the number of clusters to")
args = parser.parse_args()

# Read collective variables
colvars = read_colvar(args.colvar)

# Scale and cluster the data
X = colvars.copy().values
pipe = Pipeline(
    [
        ("scaler", StandardScaler()),
        ("clustering", HDBSCAN(min_cluster_size=args.min_samples, n_jobs=-1)),
    ]
)
model = pipe.fit(X)

# Save the model
if args.model_file:
    joblib.dump(model, args.model_file)

# Save the number of clusters
if args.nclusters_file:
    nclusters = model["clustering"].labels_.max() + 1
    np.savetxt(args.nclusters_file, [nclusters], fmt="%d")
