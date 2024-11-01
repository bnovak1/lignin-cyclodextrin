"""
Calculate fractions of each cluster and cluster grouping.

This script reads cluster labels for each frame and a JSON file containing groupings of clusters
with configurations of the lignin dimer bound near the center of the cyclodextrin.
It calculates the fractions for each cluster and cluster grouping and saves the results
to a CSV file.

Command line arguments
----------------------
--labels : str
    File with integer cluster labels for each frame, one per line.
--center_bound_clusters : str
    JSON file with groupings of clusters which are composed of configurations with the lignin dimer bound near the center of the cyclodextrin. Cluster group names are keys and cluster labels are provided as lists of integers under each cluster group key.
--outfile : str
    Output CSV file with fractions of for each cluster and cluster grouping.

Example usage
--------
python cluster_fractions.py --labels labels.txt --center_bound_clusters clusters.json --outfile output.csv

Example JSON file content
-------------------------
{
    "head:sec": [
        4,
        5,
        6,
        7,
        8
    ],
    "center:sec": [
        2
    ],
    "tail:sec": [
        3
    ]
}
"""

import argparse
import json
import numpy as np
import pandas as pd

def main(labels_file, center_bound_clusters_file, outfile):
    """
    Main function to calculate fractions of each cluster and cluster grouping.

    Parameters
    ----------
    labels_file : str
        Path to the file with cluster labels for each frame.
    center_bound_clusters_file : str
        Path to the JSON file with groupings of clusters.
    outfile : str
        Path to the output CSV file.

    Returns
    -------
    None
    """
    # Load cluster labels
    labels = np.loadtxt(labels_file, dtype=int)

    # Load center bound clusters file
    with open(center_bound_clusters_file, "r", encoding="utf-8") as jf:
        json_data = json.load(jf)

    # Create list of center bound clusters
    center_bound_clusters = []
    for cluster_group in json_data:
        center_bound_clusters += json_data[cluster_group]

    # Total number of center bound configurations
    TOTAL_CENTER_BOUND = 0

    for cluster in center_bound_clusters:
        TOTAL_CENTER_BOUND += np.sum(labels == cluster)

    # Calculate fractions for each cluster and cluster grouping. Store in a DataFrame
    fractions = pd.DataFrame(columns=["Label", "Fraction"])

    for cluster_group in json_data:

        TOTAL_GROUP_CONFIGS = 0

        for cluster in json_data[cluster_group]:

            n_configs = np.sum(labels == cluster)
            TOTAL_GROUP_CONFIGS += n_configs
            fraction = round(n_configs / TOTAL_CENTER_BOUND, 4)
            fractions = pd.concat(
                [fractions, pd.DataFrame([{"Label": cluster, "Fraction": fraction}])], ignore_index=True
            )

        fractions = pd.concat(
            [
                fractions,
                pd.DataFrame(
                    [
                        {
                            "Label": cluster_group,
                            "Fraction": round(TOTAL_GROUP_CONFIGS / TOTAL_CENTER_BOUND, 4),
                        }
                    ]
                ),
            ],
            ignore_index=True,
        )

    # Save fractions to CSV file
    fractions.to_csv(outfile, index=False)

if __name__ == "__main__":
    
    # Command line arguments
    argparser = argparse.ArgumentParser()
    argparser.add_argument(
        "--labels", type=str, help="File with cluster labels for each frame, one per line"
    )
    argparser.add_argument(
        "--center_bound_clusters",
        type=str,
        help="JSON file with groupings of clusters which are composed of configurations with the lignin dimer bound near the center of the cyclodextrin",
    )
    argparser.add_argument(
        "--outfile",
        type=str,
        help="Output CSV file with fractions of for each cluster and cluster grouping",
    )
    args = argparser.parse_args()

    # Run main function
    main(args.labels, args.center_bound_clusters, args.outfile)
