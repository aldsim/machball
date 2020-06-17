#Copyright 2013 Argonne UChicago LLC
#This file is part of Machball

import numpy as np
import machball.ballistic.viewfactors as vf

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
        2D array with the view factors for all the elements
    regions : dict
        a dictionary partitioning the structure into a series of regions

    """

    def __init__(self, areas, qij, regions=None):
        self.N = len(areas)
        self.areas = areas
        self.qij = qij
        if regions == None:
            self.regions = {}
        else:
            self.regions = regions.copy()

    def region(self, name):
        return self.regions[name]

    def get_p0(self, entrypoints):
        p0 = np.zeros(self.N)
        for i in entrypoints:
            p0[i] = self.areas[i]
        p0 /= sum(p0)
        return p0

def read_structure(filename, areafile=None):
    """Read a Structure from file

    If areafile is not defined, it assumes that
    the first column contains the normalized areas, with
    the remaining columns the view factors.

    Otherwise, it reads the first column of areafile and assumes
    that the data contained in filename are just the view factors.
    """

    data = np.loadtxt(filename)
    if areafile is None:
        return Structure(data[:,0], data[:,1:])
    else:
        areas = np.loadtxt(areafile)[:,0]
        return Structure(areas, data)



class Via(Structure):
    """Implement a circular via

    A Via is composed of `Nz` + 2 elements, with the first and last elements
    corresponding to the opening and the base, respectively.

    The view factors are calculated assuming a cosine law distribution, which
    corresponds to the behavior expected both from diffuse reemision and a
    non-directional flux of incident species.

    Parameters
    ----------
    AR : float
        Aspect ratio, defined as the depth to diameter ratio
    Nz : int
        Number of vertical sections in the discretized wall

    """

    def __init__(self, AR, Nz):

        N = Nz + 2
        regions = {"top":0, "bottom":(N-1), "wall":slice(1,N-1)}
        areas, qij = create_via(AR, Nz)
        Structure.__init__(self, areas, qij, regions)


class Trench(Structure):
    """Implement a rectangular trench

    A trench is composed of `Nz` + 2 elements, with the first and last
    elements corresponding to the opening and the base, respectively.

    The view factors are calculated assuming a cosine law distribution, which
    corresponds to the behavior expected both from diffuse reemision and a
    non-directional flux of incident species.

    Parameters
    ----------
    AR : float
        Aspect ratio, defined as the depth to width ratio
    Nz : int
        Number of vertical sections in the discretized wall

    """

    def __init__(self, AR, Nz):

        N = Nz + 2
        regions = {"top":0, "bottom":(N-1), "wall":slice(1,N-1)}
        areas, qij = create_trench(AR, Nz)
        Structure.__init__(self, areas, qij, regions)



def create_via(AR, Nz):
    """Return the areas and view factor of a circular via, where
    the vertical wall is divided into identical sections.

    Parameters
    ----------
    AR : float
        Aspect ratio, defined as the depth to diameter ratio
    Nz : int
        Number of vertical sections in the discretized wall

    Returns
    -------

    (numpy.array, numpy.array)
        Tuple with the areas (1D array), and view factors (2D array)

    """

    N = int(Nz + 2)
    qij = np.zeros((N,N))
    S0 = 0.25*np.pi
    areas = S0*np.ones(N)
    dz = AR/Nz
    areas[1:(N-1)] = np.pi*dz

    for i in range(Nz):
        qij[0,i+1] = vf.cylinderwall_to_disk(0.5, 0.5, i*dz, (i+1)*dz)
        qij[i+1,0] = areas[i+1]/areas[0]*qij[0,i+1]
        qij[N-2-i,-1] = qij[i+1,0]
        qij[-1,N-2-i] = qij[0,i+1]
    qij[-1,0] = 1-sum(qij[:,0])
    qij[0,-1] = qij[-1,0]

    for i in range(Nz):
        qij[i+1,i+1] = vf.cylinderwall_to_itself(0.5, dz)

    qdj = np.zeros(Nz-1)
    for i in range(Nz-1):
        qdj[i] = vf.cylindersection_to_section(0.5, dz, i+1)

    for i in range(Nz-1):
        for j in range(i+1,Nz):
            qij[i+1,j+1] = qdj[j-i-1]
            qij[j+1,i+1] = qij[i+1,j+1]

    return areas, qij



def create_trench(AR, Nz):
    """Return the areas and view factor of a rectangular trench, where
    the vertical wall is divided into identical sections.

    Parameters
    ----------
    AR : float
        Aspect ratio, defined as the width to diameter ratio
    Nz : int
        Number of vertical sections in the discretized wall

    Returns
    -------

    (numpy.array, numpy.array)
        Tuple with the areas (1D array), and view factors (2D array)

    """

    N = int(Nz + 2)
    qij = np.zeros((N,N))
    S0 = 1
    areas = np.ones(N)
    dz = AR/Nz
    areas[1:(N-1)] = 2*dz

    for i in range(Nz):
        qij[i+1,0] = vf.strip_to_wall(dz, i+1)
        qij[0,i+1] = qij[i+1,0]/(2*dz)
        qij[N-2-i,-1] = qij[i+1,0]
        qij[-1,N-2-i] = qij[0,i+1]
    qij[-1,0] = 1-sum(qij[:,0])
    qij[0,-1] = qij[-1,0]

    for i in range(Nz):
        qij[i+1,i+1] = vf.strip_to_strip(dz, 0)

    qdj = np.zeros(Nz-1)
    for i in range(Nz-1):
        qdj[i] = vf.strip_to_strip(dz, i+1)

    for i in range(Nz-1):
        for j in range(i+1,Nz):
            qij[i+1,j+1] = qdj[j-i-1]
            qij[j+1,i+1] = qij[i+1,j+1]

    return areas, qij

def create_taperedvia(AR, dr, Nseg):
    """Return the areas and view factor of a rectangular trench, where
    the vertical wall is divided into identical sections.

    Parameters
    ----------
    AR : float
        Aspect ratio, defined as the width to diameter ratio
    Nz : int
        Number of vertical sections in the discretized wall

    Returns
    -------

    (numpy.array, numpy.array)
        Tuple with the areas (1D array), and view factors (2D array)

    """

    dh = AR/Nseg
    tana = 0.5*dr/AR

    hi = np.arange(Nseg+1)*dh
    ri = 0.5*(1-tana*hi)
    Ai  = ri*ri
    Si = (ri[0:-1]+ri[1:])*dh

    areas = np.zeros(Nseg+2)
    areas[0] = Ai[0]
    areas[-1] = Ai[-1]
    areas[1:-1] = Si
    areas = areas/areas[0]
    Ai = Ai/Ai[0]

    Fij = np.ones((Nseg+1,Nseg+1))
    for i in range(Nseg):
        for j in range(i+1, Nseg+1):
            Fij[j,i] = vf.disk_to_shrinkcone(ri[i], dh*(j-i), tana)
            Fij[i,j] = Ai[i]*Fij[j,i]/Ai[j]

    Qij = np.zeros((Nseg+2, Nseg+2))

    for i in range(1,Nseg):
        for j in range(i+1,Nseg+1):
            q0 = Ai[i+1]*(Fij[j-1,i]-Fij[j,i])-Ai[i]*(Fij[j-1,i-1]-Fij[j,i-1])
            Qij[j,i] = q0/areas[i]
            Qij[i,j] = q0/areas[j]

    Qij[1,0] = 1-Fij[1,0]
    for i in range(2,Nseg+1):
        Qij[i,0] = Fij[i-1,0]-Fij[i,0]
    for i in range(1,Nseg):
        Qij[i,-1] = Fij[i,-1]-Fij[i-1,-1]
    Qij[-2,-1] = 1-Fij[-2,-1]
    for  i in range(1,Nseg+1):
        Qij[0,i] = areas[0]/areas[i]*Qij[i,0]
        Qij[-1,i] = areas[-1]/areas[i]*Qij[i,-1]
    Qij[-1,0] = 1-np.sum(Qij[:-1,0])
    Qij[0,-1] = 1-np.sum(Qij[1:,-1])

    for i in range(1,Nseg+1):
        Qij[i,i] = 1-np.sum(Qij[:,i])

    return areas, Qij


if __name__ == "__main__":
    a, q = create_via(10,20)
    print(np.sum(q, axis=0), np.sum(q, axis=1),  np.sum(q))
    a, q = create_trench(10,20)
    print(np.sum(q, axis=0), np.sum(q, axis=1), np.sum(q))
    a, q = create_taperedvia(20,0.5,40)
    print(np.sum(q, axis=0), np.sum(q, axis=1), np.sum(q))
