"""
This script analyzes the number of clusters formed by HDBSCAN clustering 
algorithm over a range of `min_samples` values and identifies the optimal 
`min_samples` value based on the stability of the number of clusters.

Command Line Arguments
----------------------
lignol : str
    Lignol abbreviation.
min_samples_rng : int, nargs=2
    Range of values of min_samples used in clustering.
dir : str
    Base directory name for input and output files.

Workflow
--------
1. Parse command line arguments.
2. Loop over the specified range of `min_samples` values to read the number of clusters from files.
3. Identify the optimal `min_samples` value based on the position of the last increase in the number of clusters.
4. Determine the optimal `min_samples` value based on the last increase in the number of clusters.
5. Save the optimal `min_samples` value to a file.
6. Plot the number of clusters as a function of `min_samples` and save the plot in both HTML and PNG formats.
"""

import argparse
from pathlib import Path

import matplotlib.pyplot as plt
import plotly.express as px
import numpy as np


def main():
    """
    Main function to execute the script.

    Parses command line arguments, reads the number of clusters for each
    `min_samples` value, identifies the optimal `min_samples` value, and
    generates plots.

    Returns
    -------
    None
    """
    # Command line inputs
    parser = argparse.ArgumentParser()
    parser.add_argument("--lignol", type=str, help="Lignol abbreviation")
    parser.add_argument(
        "--min_samples_rng",
        type=int,
        nargs=2,
        help="Range of values of min_samples used in clustering",
    )
    parser.add_argument("--dir", type=str, help="Base directory name for input and output files")
    args = parser.parse_args()

    # Loop over min_samples values to read in the number of clusters
    n_clusters = []
    for min_samples in range(args.min_samples_rng[0], args.min_samples_rng[1] + 1):

        infile = Path(args.dir, "HDBSCAN_model", f"nclusters_{min_samples}.dat")
        n_clusters.append([min_samples, int(np.loadtxt(infile, dtype=int))])

    n_clusters = np.array(n_clusters)

    # Find the last position where the number of clusters increases
    diff = np.diff(n_clusters[:, 1])
    idx_increase = np.where(diff > 0)[0][-1]
                
    # min_samples taken at the next occurrence of n_clusters corresponding to the value at the last increase in n_clusters
    n_clusters_increase = n_clusters[idx_increase, 1]
    idx = np.where(n_clusters[:, 1] == n_clusters_increase)[0]
    idx_diff = idx - idx_increase
    min_samples = n_clusters[idx[np.where(idx_diff > 0)[0][0]], 0]

    # Save the min_samples value to a file
    np.savetxt(Path(args.dir, "min_samples_best.dat"), [min_samples], fmt="%d")

    # Plot the number of clusters as a function of min_samples with plotly
    fig = px.line(
        x=n_clusters[:, 0],
        y=n_clusters[:, 1],
        markers=True,
        labels={"x": "min_samples", "y": "Number of clusters"},
    )

    # Save the plot in HTML format
    fig.write_html(Path(args.dir, "min_samples.html"))

    # Plot the number of clusters as a function of min_samples with matplotlib
    plt.plot(n_clusters[:, 0], n_clusters[:, 1], "k.-")

    # Use a red circle to highlight the min_samples value used
    idx = n_clusters[:, 0] == min_samples
    plt.plot(
        n_clusters[idx, 0],
        n_clusters[idx, 1],
        "ro",
        mfc="none",
        label=f"min_samples = {min_samples} used\n{n_clusters[idx, 1][0]} clusters",
    )

    # Labels and legend
    plt.xlabel("min_samples")
    plt.ylabel("Number of clusters")
    plt.legend()

    # Save the plot in PNG format
    plt.savefig(Path(args.dir, "min_samples.png"), dpi=300, bbox_inches="tight")


if __name__ == "__main__":
    main()
