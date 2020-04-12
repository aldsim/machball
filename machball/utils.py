#Copyright 2019 Argonne UChicago LLC
#This file is part of MachBall
"""
Utility functions
"""
import numpy as np


# Boltzman constant
kB = 1.38e-23
# atomic mass unit in kg
amu = 1.66e-27

def vth(T, M):
    """Mean thermal velocity

    Parameters
    ----------

    T : float
        Temperature in K
    M : float
        Molecular mass in atomic mass units

    """
    return np.sqrt(8*kB*T/(amu*M*np.pi))
