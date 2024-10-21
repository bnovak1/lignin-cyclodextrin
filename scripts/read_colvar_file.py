import pandas as pd


def read_colvar(colvar_file):
    """
    Read the colvar file and extract the collective variables; 'dnorm', 'dtang', and 'orient'

    Parameters
    ----------
    colvar_file : str
        Path to the colvar file.

    Returns
    -------
    pd.DataFrame
        A DataFrame containing the columns 'dnorm', 'dtang', and 'orient' from the colvar file.
    """

    # Read the colvar file
    colvars = pd.read_csv(
        colvar_file, sep="\s+", header=None, skiprows=1, names=["dnorm", "dtang", "orient"]
    )

    return colvars
