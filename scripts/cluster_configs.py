"""
Read cluster labels and a molecular dynamics trajectory, and save the trajectory frames corresponding to a specified cluster id to a new XTC file.

Usage
-----
python cluster_configs.py --labels <labels_file> --gro <gro_file> --xtc_in <input_xtc_file> --xtc_out <output_xtc_file> --cluster_id <cluster_id>

Parameters
----------
labels : str
    Path to the input file containing cluster labels.
gro : str
    Path to the GRO file containing lignin and cyclodextrin structures.
xtc_in : str
    Path to the input XTC file containing the trajectory of lignin and cyclodextrin.
xtc_out : str
    Path to the output XTC file where the trajectory frames for the specified cluster will be saved.
cluster_id : int
    The ID of the cluster for which the trajectory frames will be saved.

Returns
-------
None
"""
import argparse
from pathlib import Path
import sys

import MDAnalysis as mda
import numpy as np

sys.path.append(Path("../scripts"))
from read_colvar_file import read_colvar


# Parse command line arguments
parser = argparse.ArgumentParser()
parser.add_argument("--labels", type=str, help="Input file with cluster labels")
parser.add_argument("--gro", type=str, help="gro file with lignin & cyclodextrin")
parser.add_argument("--xtc_in", type=str, help="Input xtc file with lignin & cyclodextrin")
parser.add_argument("--xtc_out", type=str, 
                    help="Input xtc file with lignin & cyclodextrin for a cluster")
parser.add_argument("--cluster_id", type=int, help="Cluster id")
args = parser.parse_args()

# Read in cluster labels
labels = np.loadtxt(args.labels, dtype=int)

# Read in the trajectory
u = mda.Universe(args.gro, args.xtc_in)

# Save the trajectory frames for the cluster
idx = np.where(labels == args.cluster_id)[0]
u.atoms.write(args.xtc_out, frames=u.trajectory[idx])

