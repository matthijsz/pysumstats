import pandas as pd
import numpy as np


def cov_matrix_from_phenotype_file(dataframe, phenotypes=None):
    '''Function to generate a covariance (cov_Z) matrix from a phenotype file.

    :param dataframe: pd.Dataframe containing phenotypic data
    :type dataframe: pd.Dataframe
    :param phenotypes: list of phenotypes to include
    :type phenotypes: list.
    :return: pd.Dataframe of covariance matrix.

    '''
    if phenotypes is None:
        phenotypes = list(dataframe.columns)
    cov_matrix = pd.DataFrame(0, index=phenotypes, columns=phenotypes)
    for p1 in phenotypes:
        cov_matrix.loc[p1, p1] = 1
        n1 = len(dataframe.loc[dataframe[p1].notna(), :])
        for p2 in [x for x in phenotypes if x != p1]:
            ns = len(dataframe.loc[(dataframe[p1].notna() & dataframe[p2].notna()), :])
            n2 = len(dataframe.loc[dataframe[p2].notna(), :])
            r = dataframe[[p1, p2]].corr().iloc[0, 1]
            cov_matrix.loc[p1, p2] = (r * ns) / np.sqrt((n1 * n2))
    return cov_matrix
