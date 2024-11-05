"""
Use the elbow method on the average distance to the kth nearest neighbor to determine max_eps for OPTICS.
"""

import argparse
from pathlib import Path
import sys

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from sklearn.neighbors import NearestNeighbors
sys.path.append(Path("../scripts"))
from read_colvar_file import read_colvar

# Command line inputs
parser = argparse.ArgumentParser()
parser.add_argument("--lignol", type=str, help="Lignol abbreviation")
parser.add_argument("--colvar", type=str, help="Input file with collective variables data")
parser.add_argument("--outdir", type=str, help="Output directory name")
args = parser.parse_args()

# Read data
X = read_colvar(args.colvar, args.lignol).values
    
# Calculate average distances to the kth nearest neighbors
nbrs = NearestNeighbors(n_neighbors=20).fit(X)
distances, _ = nbrs.kneighbors(X)
distances_mean = np.mean(distances, axis=0)

# Plot the average distances to the kth nearest neighbors
k = np.arange(20)
plt.plot(k[1:], distances_mean[1:])
plt.xlabel("k")
plt.ylabel("Average distance to the kth nearest neighbor")
# Set the major ticks at integers
plt.xticks(np.arange(1, 20, 1))
plt.savefig(Path(args.outdir, "elbow.png"), dpi=300, bbox_inches="tight")

# Estimate the elbow location
k2 = k[2:-1]
deriv2 = np.abs(np.diff(distances_mean[1:], 2))
idx = np.argmax(deriv2) + 2
eps = distances_mean[idx]
max_eps = 10 * eps
np.savetxt(Path(args.outdir, "max_eps.txt"), [max_eps])
