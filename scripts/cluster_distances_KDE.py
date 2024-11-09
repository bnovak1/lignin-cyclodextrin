"""
* Reads the cyclodextrin COM to lignin dimer (COM, head, tail, ends) normal distances for each frame from a file created by PLUMED (colvar input)
* Reads the cluster labels for each frame from a file (clust_labels input)
* Reads the cluster group information from a specified JSON file (center_bound_clusters input)
* Calculates the normal distance Kernel Density Estimates (KDEs) for each cluster group
* Plots the KDEs in a single figure for each cluster group
* Saves the KDE figure for each cluster group to PNG files in the specified output directory (outdir input)

Functions
---------
main(args):

Usage
-----
To use this script from the command line, provide the path to the input file containing the collective variables and the path to the output file for the KDE plots:
    python cluster_distances_KDE.py --colvar path/to/colvar/file --output path/to/output/plot
    
Example JSON file content
-------------------------
{
    "head-sec": [
        4,
        5,
        6,
        7,
        8
    ],
    "center-sec": [
        2
    ],
    "tail-sec": [
        3
    ]
}
"""

import argparse
import json
from pathlib import Path
import sys

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import seaborn as sns

sys.path.append("../scripts")
from read_colvar_file import read_colvar


def main(arguments):
    """
    Read collective variables, calculate their Kernel Density Estimates (KDEs), and plot the KDEs.

    Parameters
    ----------
    args : argparse.Namespace
        Command-line arguments containing the following attributes:
        - lignol : str
            Abbreviation for the lignin dimer derivative.
        - colvar : str
            Path to the PLUMED file containing the collective variables.
        - clust_labels : str
            Path to the file containing the cluster labels.
        - center_bound_clusters : str
            JSON file with groupings of clusters which are composed of configurations with the lignin dimer bound near the center of the cyclodextrin. Cluster group names are keys and cluster labels are provided as lists of integers under each cluster group key.
        - cluster_group : str
            Name of the cluster group to analyze.
        - outdir : str
            Directory to save the output plots to.

    Returns
    -------
    None
    """

    # Read normal distances from the PLUMED file
    distances = read_colvar(arguments.colvar, cv_list=["dnorm", "dhead", "dtail"])

    # Read cluster labels
    clust_labels = np.loadtxt(arguments.clust_labels, dtype=int)

    # Load center bound clusters JSON file
    with open(arguments.center_bound_clusters, "r", encoding="utf-8") as jf:
        json_data = json.load(jf)

    # Loop over each cluster label in the cluster group
    for cluster in json_data[arguments.cluster_group]:

        # Filter the normal distances for the cluster label
        # and concatenate them to cluster_distances DataFrame
        try:
            cluster_distances = distances[clust_labels == cluster]
        except NameError:
            cluster_distances = pd.concat([cluster_distances, distances[clust_labels == cluster]])

        # Create a figure
        plt.figure()

        # Plot KDEs for each column
        bandwidth = 0.15
        if arguments.lignol == "GG_BB":
            
            # Plot KDE for COM
            data = cluster_distances["dnorm"]
            sns.kdeplot(data, bw_method=bandwidth, label="COM")
            
            # Plot KDE for ends (head and tail concatenated)
            data = pd.concat([cluster_distances["dhead"], cluster_distances["dtail"]], ignore_index=True)
            sns.kdeplot(data, bw_method=bandwidth/2, label="ends")
            
        else:
            
            for column in cluster_distances.columns:
                data = cluster_distances[column]
                sns.kdeplot(data, bw_method=bandwidth, label=column)

        # Add legend and labels
        if arguments.lignol == "GG_BB":
            plt.legend()
        else:
            plt.legend(["COM", "head", "tail"])
        
        plt.xlabel("Normal distance (nm)")
        plt.ylabel("Density")
        
        xlims = plt.xlim()
        ylims = plt.ylim()
        if arguments.lignol == "GG_BB":
            text = f"COM Bandwidth = {bandwidth:.3f} nm\nends Bandwidth = {bandwidth/2:.3f} nm"
        else:
            text = f"Bandwidth = {bandwidth:.3f} nm"
        plt.text(
            xlims[0] + 0.01 * (xlims[1] - xlims[0]),
            ylims[1] - 0.01 * (ylims[1] - ylims[0]),
            text,
            ha="left",
            va="top",
        )

        # Save the figure and close it
        outfile = Path(
            arguments.outdir,
            f"cluster_distances_KDE_{arguments.cluster_group}.png",
        )
        plt.savefig(outfile, dpi=300, bbox_inches="tight")
        plt.close()

        # Remove the cluster_distances DataFrame
        del cluster_distances


if __name__ == "__main__":

    # Command-line inputs
    parser = argparse.ArgumentParser()
    parser.add_argument("--lignol", type=str, help="Lignin dimer abbreviation")
    parser.add_argument(
        "--colvar", type=str, help="Input PLUMED file with collective variables data for each frame"
    )
    parser.add_argument(
        "--clust_labels", type=str, help="Input file with cluster labels for each frame"
    )
    parser.add_argument(
        "--center_bound_clusters",
        type=str,
        help="JSON file with groupings of clusters which are composed of configurations with the lignin dimer bound near the center of the cyclodextrin",
    )
    parser.add_argument("--cluster_group", type=str, help="Cluster group name")
    parser.add_argument("--outdir", type=str, help="Output directory for KDE plot PNG file")
    args = parser.parse_args()

    main(args)
