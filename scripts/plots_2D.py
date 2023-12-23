import argparse
from pathlib import Path
import matplotlib.pyplot as plt
import pandas as pd

# Command line inputs
parser = argparse.ArgumentParser()
parser.add_argument("--lignol", type=str, help="Lignol abbreviation")
parser.add_argument("--colvar", type=str, help="Input file with collective variables data")
parser.add_argument("--outdir", type=str, help="Output directory name")
args = parser.parse_args()

# Read data
with open(args.colvar, "r", encoding="utf-8") as fid:
    if args.lignol == "GG_BB":
        column_names = fid.readline().split()[-3:]
    else:
        column_names = fid.readline().split()[-4:]

colvar_data = pd.read_csv(args.colvar, sep="\s+", header=None, skiprows=1, names=column_names)

# Plots
plt.plot(colvar_data["dnorm"], colvar_data["dtang"], "ko", markersize=1, alpha=0.5)
plt.xlabel("dnorm")
plt.ylabel("dtang")
plt.savefig(Path(args.outdir, "dnorm_dtang.png"), dpi=300, bbox_inches="tight")
plt.close()

if args.lignol == "GG_BB":
    
    Path(args.outdir, "dnorm_orient.png").touch()
    Path(args.outdir, "dtang_orient.png").touch()

else:
    
    plt.plot(colvar_data["dnorm"], colvar_data["orient"], "ko", markersize=1, alpha=0.5)
    plt.xlabel("dnorm")
    plt.ylabel("orient")
    plt.savefig(Path(args.outdir, "dnorm_orient.png"), dpi=300, bbox_inches="tight")
    plt.close()

    plt.plot(colvar_data["dtang"], colvar_data["orient"], "ko", markersize=1, alpha=0.5)
    plt.xlabel("dtang")
    plt.ylabel("orient")
    plt.savefig(Path(args.outdir, "dtang_orient.png"), dpi=300, bbox_inches="tight")
    plt.close()
