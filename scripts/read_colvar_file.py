import numpy as np
import pandas as pd


def read_colvar(colvar_file, cv_list=["dnorm", "dtang", "orient"]):
    """
    Read the colvar file and extract the specified collective variables in cv_list.

    Parameters
    ----------
    colvar_file : str
        Path to the colvar file.
    cv_list : list of str, default=["dnorm", "dtang", "orient"]
        List of collective variables to extract from the colvar file.

    Returns
    -------
    pd.DataFrame
        A DataFrame containing the columns in cv_list from the colvar file.
    """

    # Read the colvar file header to get the column names
    with open(colvar_file, "r", encoding="utf-8") as f:
        cols = f.readline().split()[2:]

    # Read the colvar file and save the data to a DataFrame
    colvars = pd.DataFrame(np.loadtxt(colvar_file, skiprows=1), columns=cols)

    # Return the requested collective variables
    return colvars[cv_list]
