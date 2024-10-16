
import argparse
from pathlib import Path
import sys

import joblib
import matplotlib.pyplot as plt
import numpy as np
import plotly.express as px
import plotly.io as pio

sys.path.append(Path("../scripts"))
from read_colvar_file import read_colvar

def main():
    """
    3D plot of the clusters for the best HDBSCAN clustering model.
    
    Parameters
    ----------
    lignol : str
        Lignol abbreviation.
    colvar : str
        Input file with collective variables data.
    model : str
        Name of the joblib file with the best clustering model.
    dir : str
        Input and output directory name.

    Returns
    -------
    None
    """
    # Command line inputs
    parser = argparse.ArgumentParser()
    parser.add_argument("--lignol", type=str, help="Lignol abbreviation")
    parser.add_argument("--colvar", type=str, help="Input file with collective variables data")
    parser.add_argument(
        "--model", type=str, help="Name of the joblib file with the best clustering model"
    )
    parser.add_argument("--dir", type=str, help="Input and output directory name")
    args = parser.parse_args()

    # Read data
    infile = Path(args.model)
    pipe = joblib.load(infile)
    clustering = pipe["clustering"]
    colvars = read_colvar(args.colvar, args.lignol)
    colvars["label"] = clustering.labels_.astype(int)
    idx_cluster = colvars["label"] >= 0
    colvars = colvars[idx_cluster]
    colvars["label"] = colvars["label"].astype("category")

    # Plot the clusters
    colors = ["red", "green", "blue", "purple", "orange", "cyan", "brown", "black", "pink", "yellow"]

    fig = px.scatter_3d(
        colvars,
        x="dnorm",
        y="dtang",
        z="orient",
        color_discrete_sequence=colors,
        color="label",
        opacity=1.0,
        size=np.ones(colvars.shape[0]),
        size_max=5,
    )

    axis_label_font = {"size": 16}
    tick_label_font = {"size": 13}

    if args.lignol == "GG_BB":
        fig.update_layout(
            scene={
                "xaxis": {
                    "title_text": "d<sub>norm</sub>",
                    "title_font": axis_label_font,
                    "tickfont": tick_label_font,
                },
                "yaxis": {
                    "title_text": "d<sub>tang</sub>",
                    "title_font": axis_label_font,
                    "tickfont": tick_label_font,
                },
                "zaxis": {
                    "title_text": "cos<sup>2</sup>θ",
                    "title_font": axis_label_font,
                    "tickfont": tick_label_font,
                },
                "camera": {
                    "up": {"x": 0, "y": 0, "z": 1},
                    "center": {"x": 0, "y": 0, "z": 0},
                    "eye": {"x": 1.5, "y": 1.5, "z": 1.5},
                },
            },
            showlegend=True,
        )
    else:
        fig.update_layout(
            scene={
                "xaxis": {
                    "title_text": "d<sub>norm</sub>",
                    "title_font": axis_label_font,
                    "tickfont": tick_label_font,
                },
                "yaxis": {
                    "title_text": "d<sub>tang</sub>",
                    "title_font": axis_label_font,
                    "tickfont": tick_label_font,
                },
                "zaxis": {
                    "title_text": "cos θ",
                    "title_font": axis_label_font,
                    "tickfont": tick_label_font,
                },
                "camera": {
                    "up": {"x": 0, "y": 0, "z": 1},
                    "center": {"x": 0, "y": 0, "z": 0},
                    "eye": {"x": 1.5, "y": 1.5, "z": 1.5},
                },
            },
            showlegend=True,
        )

    fig.update_traces(marker_line_color="rgba(0,0,0,0)")

    fig.update_layout(plot_bgcolor="rgba(0,0,0,0)")
    fig.write_html(Path(args.dir, "clusters.html"))

    fig.update_layout(showlegend=False)
    pio.write_image(fig, Path(args.dir, "clusters.png"), scale=4)

    # Save the cluster labels
    np.savetxt(Path(args.dir, "cluster_labels.dat"), clustering.labels_, fmt="%d")

    # Save the number of clusters
    nclusters = colvars["label"].nunique()
    np.savetxt(Path(args.dir, "nclusters.dat"), [nclusters], fmt="%d")

if __name__ == "__main__":
    main()