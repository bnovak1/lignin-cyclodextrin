"""
2D scatter plots and KDE plots of collective variables.
"""

import argparse
from pathlib import Path
import sys

import matplotlib
import matplotlib.pyplot as plt
import numpy as np
import plotly.express as px
import plotly.io as pio
from scipy.spatial import ConvexHull, Voronoi
import seaborn as sns

sys.path.append(Path("../scripts"))
from read_colvar_file import read_colvar

# Command line inputs
parser = argparse.ArgumentParser()
parser.add_argument("--lignol", type=str, help="Lignol abbreviation")
parser.add_argument("--colvar", type=str, help="Input file with collective variables data")
parser.add_argument("--outdir", type=str, help="Output directory name")
parser.add_argument("--dimension", type=int, help="Dimension of the plot(s)")
args = parser.parse_args()

# Read data
colvars = read_colvar(args.colvar, args.lignol)

if args.dimension == 2:

    # non-interactive backend
    matplotlib.use("Agg")

    # Normal distance vs tangential distance - scatter plot
    plt.plot(colvars["dnorm"], colvars["dtang"], "ko", markersize=1, alpha=0.5)
    plt.xlabel("$d_{norm}$")
    plt.ylabel("$d_{tang}$")
    plt.savefig(Path(args.outdir, "dnorm_dtang.png"), dpi=300, bbox_inches="tight")
    plt.close()

    # Color map
    color_map = plt.cm.get_cmap('nipy_spectral')
    reversed_color_map = color_map.reversed()

    # Normal distance vs tangential distance - KDE plot
    sns.kdeplot(
        data=colvars, x="dnorm", y="dtang",
        fill=True,
        thresh=0, levels=100,
        cmap=reversed_color_map, cbar=True,
    )
    plt.xlabel("$d_{norm}$")
    plt.ylabel("$d_{tang}$")
    plt.savefig(Path(args.outdir, "dnorm_dtang_kde.png"), dpi=300, bbox_inches="tight")
    plt.close()
    
    # Normal distance vs orientation - scatter plot
    plt.plot(colvars["dnorm"], colvars["orient"], "ko", markersize=1, alpha=0.5)
    plt.xlabel("$d_{norm}$")
    if args.lignol == "GG_BB":
        plt.ylabel("$\cos^2 \\theta$")
    else:
        plt.ylabel("$\cos \\theta$")
    plt.savefig(Path(args.outdir, "dnorm_orient.png"), dpi=300, bbox_inches="tight")
    plt.close()

    # Normal distance vs orientation - KDE plot
    sns.kdeplot(
        data=colvars, x="dnorm", y="orient",
        fill=True,
        thresh=0, levels=100,
        cmap=reversed_color_map, cbar=True,
    )
    plt.xlabel("$d_{norm}$")
    if args.lignol == "GG_BB":
        plt.ylabel("$\cos^2 \\theta$")
    else:
        plt.ylabel("$\cos \\theta$")
    plt.savefig(Path(args.outdir, "dnorm_orient_kde.png"), dpi=300, bbox_inches="tight")
    plt.close()
    
    # Tangential distance vs orientation - scatter plot
    plt.plot(colvars["dtang"], colvars["orient"], "ko", markersize=1, alpha=0.5)
    plt.xlabel("$d_{tang}$")
    if args.lignol == "GG_BB":
        plt.ylabel("$\cos^2 \\theta$")
    else:
        plt.ylabel("$\cos \\theta$")
    plt.savefig(Path(args.outdir, "dtang_orient.png"), dpi=300, bbox_inches="tight")
    plt.close()

    # Tangential distance vs orientation - KDE plot
    sns.kdeplot(
        data=colvars, x="dtang", y="orient",
        fill=True,
        thresh=0, levels=150,
        cmap=reversed_color_map, cbar=True,
    )
    plt.xlabel("$d_{tang}$")
    if args.lignol == "GG_BB":
        plt.ylabel("$\cos^2 \\theta$")
    else:
        plt.ylabel("$\cos \\theta$")
    plt.savefig(Path(args.outdir, "dtang_orient_kde.png"), dpi=300, bbox_inches="tight")
    plt.close()

elif args.dimension == 3:
    
    if args.lignol == "GG_BB":
        voro = Voronoi(colvars[["dnorm", "dtang"]].values)
    else:    
        voro = Voronoi(colvars[["dnorm", "dtang", "orient"]].values)
    
    # Compute volumes of Voronoi cells
    # https://stackoverflow.com/questions/19634993/volume-of-voronoi-cell-python
    vol = np.zeros(voro.npoints)
    for i, reg_num in enumerate(voro.point_region):
        indices = voro.regions[reg_num]
        if -1 in indices: # unbounded cells
            vol[i] = np.inf
        else:
            vol[i] = ConvexHull(voro.vertices[indices]).volume

    # Keep only the lower 50% of the volumes
    colvars["volume"] = vol
    idx = vol < np.percentile(vol, 50)
    df = colvars.loc[idx].copy()

    fig = px.scatter_3d(
        df,
        x="dnorm",
        y="dtang",
        z="orient",
        opacity=1.0,
        size=np.ones(df.shape[0]),
        size_max=3
    )
    
    axis_label_font = {"size": 16}
    tick_label_font = {"size": 13}
    
    if args.lignol == "GG_BB":
        orient_str = "cos<sup>2</sup> θ"
    else:
        orient_str = "cos θ"
        
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
                "title_text": orient_str, 
                "title_font": axis_label_font, 
                "tickfont": tick_label_font
            },
            "camera": {
                "up": {"x": 0, "y": 0, "z": 1},
                "center": {"x": 0, "y": 0, "z": 0},
                "eye": {"x": 1.5, "y": 1.5, "z": 1.5}
            }
        },
        showlegend = False
    )
    fig.update_traces(marker_line_color="rgba(0,0,0,0)")
    fig.update_layout(plot_bgcolor="rgba(0,0,0,0)")

    fig.write_html(Path(args.outdir, "high_density_scatter.html"))
    pio.write_image(fig, Path(args.outdir, "high_density_scatter.png"), scale=4)