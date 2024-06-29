#Copyright 2013 Argonne UChicago LLC
#This file is part of Machball

from . import viewfactors as vf

import numpy as np
import pickle


class Structure:
    """Implement a geometrical feature as an array of areas and view factors.

    A Structure is just an array of areas and view factors. It
    optionally allows to define regions as either individual indices or as
    slices.

    Parameters
    ----------
    areas : numpy.array
        1D array with the areas of each element in the structure
    qij : numpy.array
        2D array with the view factors for all the elements. qij[i,j] represents
        the view factor of section i from section j
    regions : dict of list of ints, optional
        a dictionary partitioning the structure into a series of regions
    elements : iterable object, optional
        a list of elements containing geometrical information on each structure

    """

    def __init__(self, areas, qij, regions=None, coordinates=None):
        self.N = len(areas)
        self.areas = areas
        self.qij = qij
        if regions is None:
            self.regions = {'default':[i for i in range(self.N)]}
        else:
            self.regions = regions.copy()

    def region(self, name):
        return self.regions[name]
    
    def save(self, filename, mode="pickle", areafile=None):
        return save_structure(filename, self, mode, areafile)


def save_structure(filename, st, mode="pickle", areafile=None):
    """Save a structure to a file

    There are two options: saving the structure as a pickle,
    which retains all metadata (entrypoints, coordinates, regions),
    or just the view factors and areas as numpy arrays. In this
    second case, if areafile is provided the view factors and areas are
    saved in different files. Otherwise, the areas are saved as a first
    column of a NxN+1 array.

    Parameters
    ----------
    filename : str
        file name
    st : Structure
        the structure to be saved
    mode : {'pickle', 'numpy'}
        If "pickle" saves it as pickle.
    areafile : str, optional
        Optional additional filename for storing the areas
        when saving as numpy array.

    """

    if mode=="pickle":
        with open(filename, "wb") as f:
            pickle.dump(st, f)
    elif mode=="numpy":
        if areafile is None:
            savedata = np.zeros((st.N, st.N+1))
            savedata[:,1:] = st.qij
            savedata[:,0] = st.areas
            np.savetxt(filename, savedata)
        else:
            np.savetxt(filename, st.qij)
            np.savetxt(areafile, st.areas)
    else:
        raise ValueError("mode %s not recognized" % mode)

    

def read_structure(filename, mode="pickle", areafile=None):
    """Read a Structure from file

    There are two options: reading the structure as a pickle or
    as numpy arrays. In this
    second case, if an optional areafile is provided the areas are
    retrieved from that file. Otherwise, the areas are saved as a first
    column of a NxN+1 array in filename.

    Parameters
    ----------

    filename : str
        file name
    mode : {'pickle', 'numpy'}
        If "pickle" saves it as pickle.
    areafile : str, optional
        Optional additional filename for storing the areas
        when saving as numpy array.

    Returns
    -------

    st : Structure
        structure read from file

    """

    if mode == "pickle":
        with open(filename, "rb") as f:
            st = pickle.load(f)
        return st
    elif mode == "numpy":
        data = np.loadtxt(filename)
        if areafile is None:
            return Structure(data[:,0], data[:,1:])
        else:
            areas = np.loadtxt(areafile)
            return Structure(areas, data)
    else:
        raise ValueError("mode %s not recognized" % mode)
