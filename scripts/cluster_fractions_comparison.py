import argparse

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd


def main(arguments):
    """
    Compare and plot cluster fractions from case using 3 CVs and HDBSCAN with previous work using many CVs + PCA + DBSCAN.

    Parameters
    ----------
    arguments : argparse.Namespace
        Command line arguments containing the following attributes:
        - lignol : str
            Abbreviation for the lignin dimer derivative.
        - fractions : str
            Path to the CSV file containing the new cluster fractions.
        - fractions_previous : str
            Path to the CSV file containing the previous cluster fractions.
        - outfile : str
            Path to save the output plot.

    Returns
    -------
    None
    """

    # Load cluster fractions for the 3 CVs + HDBSCAN case
    fractions = pd.read_csv(arguments.fractions)
    fractions.set_index("Label", inplace=True)

    # Load cluster fractions from previous work using many CVs + PCA + DBSCAN
    # These are in relative to all configurations, not just the ones with the lignin dimer bound near the center of the cyclodextrin. They must be renormalized.
    fractions_previous = pd.read_csv(arguments.fractions_previous, sep=r"\s+")

    # Match the cluster group names. Renormalize the previous fractions for comparison.
    if arguments.lignol == "GG":

        fractions_previous.rename(index={"head": "head-sec", "tail": "tail-sec"}, inplace=True)
        fractions_previous.loc["center-sec"] = 0.0
        fractions_previous = fractions_previous.loc[["head-sec", "center-sec", "tail-sec"]]

        fractions = fractions.loc[["head-sec", "center-sec", "tail-sec"]]

    elif arguments.lignol == "TGG":

        fractions_previous.rename(
            index={"head": "head-sec", "tail": "tail-sec", "center": "center-sec"}, inplace=True
        )
        fractions_previous.loc["tail-pri"] = 0.0
        fractions_previous = fractions_previous.loc[["head-sec", "center-sec", "tail-sec"]]

        fractions = fractions.loc[["head-sec", "center-sec", "tail-sec"]]

    elif arguments.lignol == "GG_BB":

        fractions_previous.rename(index={"head": "center"}, inplace=True)
        fractions_previous.loc["end-sec"] = 0.0
        fractions_previous = fractions_previous.loc[["center", "end-sec"]]

        fractions = fractions.loc[["center", "end-sec"]]

    fractions_previous.rename(columns={"probability": "Fraction"}, inplace=True)
    fractions_previous /= fractions_previous.Fraction.sum()

    # Create a figure
    plt.figure()

    # Plot fractions and fractions_previous with offset bars
    plt.bar(np.arange(len(fractions)), fractions.Fraction, width=0.4, label="3 CVs + HDBSCAN")
    plt.bar(
        np.arange(len(fractions)) + 0.4,
        fractions_previous.Fraction,
        width=0.4,
        label="Many CVs + PCA + DBSCAN",
    )
    
    # Change x tick labels to cluster group names
    plt.xticks(np.arange(len(fractions)) + 0.2, fractions.index)

    # Add legend and labels
    plt.legend()
    plt.xlabel("Cluster group")
    plt.ylabel("Proportion")

    # Save the figure and close it
    plt.savefig(arguments.outfile, dpi=300, bbox_inches="tight")
    plt.close()


if __name__ == "__main__":

    parser = argparse.ArgumentParser()
    parser.add_argument("--lignol", type=str, help="Lignin dimer abbreviation")
    parser.add_argument(
        "--fractions", type=str, help="CSV file with cluster fractions for 3 CVs + HDBSCAN"
    )
    parser.add_argument(
        "--fractions_previous",
        type=str,
        help="File with cluster fractions from previous work using many CVs + PCA + DBSCAN",
    )
    parser.add_argument(
        "--outfile", type=str, help="Output bar plot PNG file with cluster fractions comparison"
    )
    args = parser.parse_args()

    main(args)
