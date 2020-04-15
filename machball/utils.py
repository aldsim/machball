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
