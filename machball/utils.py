#Copyright 2013 Argonne UChicago LLC
#This file is part of Machball
"""
This module contains a few utility functions and constants that are used
in Machball.

"""
import numpy as np


kB = 1.38e-23
"""Boltzmann constant in SI units"""

amu = 1.66e-27
"""Atomic mass unit in kg"""

densities = {
    "Al2O3": 3.95,
    "HfO2": 9.68,
    "TiO2": 3.78,
    "ZrO2": 5.68,
    "ZnO": 5.61,
}
"""Densities in g/cm3"""

molecularmasses = {
    "Al2O3": 101.96,
    "HfO2": 210.49,
    "TiO2": 79.87,
    "ZrO2": 123.22,
    "ZnO": 81.41,
}
"""Mass in mols per gram"""

def sitearea(M, density, gpc, nmol=1):
    """Average area of a surface site

    Calculate the average area of a surface site

    Parameters
    ----------

    M : float
        Molecular mass in atomic mass units
    density : float
        Density of the film, in g/cm3
    gpc : float
        Growth per cycle, in Angstroms
    nmol : int, optional (default 1)
        Number of precursor molecules per unit formula of the solid

    Returns
    -------
    float
        Average area of a surface site in sq. meters

    """

    masscm2 = density*gpc*1e-8
    molcm2 = masscm2/M*6.022e23
    return 1e-4/(nmol*molcm2)

def sitearea_fromqcm(M, mpc, nmol=1):
    """Average area of a surface site

    Calculate the average area of a surface site from qcm data

    Parameters
    ----------

    M : float
        Molecular mass in atomic mass units
    mpc : float
        Mass per cycle in  ng/cm2
    nmol : int, optional (default 1)
        Number of precursor molecules per unit formula of the solid

    Returns
    -------
    float
        Average area of a surface site in sq. meters

    """

    return M/(mpc*1e-5*6.022e23*nmol)


def vth(T, M):
    r"""Mean thermal velocity

    Calculate the mean thermal velocity, defined as:

    .. math::
        v_{th} = \sqrt{\frac{8k_BT}{M}}

    Parameters
    ----------

    T : float
        Temperature in K
    M : float
        Molecular mass in atomic mass units

    Returns
    -------
    float
        Mean thermal velocity

    """
    return np.sqrt(8*kB*T/(amu*M*np.pi))


def save_saturation(filename, times, cov, csv=False):
    """Save a saturation curve in a text file

    The saturation curve is saved as a two column text file.

    Parameters
    ----------

    filename : str
        Name of the file
    times : numpy.array
        Time column
    cov : numpy.array
        Coverage as a function  of time
    csv : bool, optional (default False)
        Save as comma separated (default is space separated)
    """

    ntimes = times.shape[0]
    data = np.zeros((ntimes,2))
    data[:,0] = times
    data[:,1] = cov
    if csv:
        np.savetxt(filename, data)
    else:
        np.savetxt(filename, data, delimiter=",")



def save_saturationprofile(filename, times, covs, csv=False):
    """Save a time series of coverage profiles as a text file

    The saturation curve is saved as a N+1 column text file, where N is
    the number of discrete zones in the feature

    Parameters
    ----------

    filename : str
        Name of the file
    times : numpy.array
        Time column
    cov : numpy.array
        2D array with coverage profiles as a function of time
    csv : bool, optional (default False)
        Save as comma separated (default is space separated)
    """

    ntimes, nzones = covs.shape

    data = np.zeros((ntimes,nzones+1))
    data[:,0] = times
    data[:,1:] = covs
    if csv:
        np.savetxt(filename, data)
    else:
        np.savetxt(filename, data, delimiter=",")
