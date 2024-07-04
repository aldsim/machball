#Copyright 2013 Argonne UChicago LLC
#This file is part of Machball

from . import viewfactors as vf
from ..structure import Structure

import numpy as np
import pickle


class Feature(Structure):
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
    entrypoints : list of ints or strings, optional
        a list containing the indices or name of regions corresponding to entry points of
        the feature. Defaults to 0.
    regions : dict of list of ints, optional
        a dictionary partitioning the structure into a series of regions
    coordinates : iterable object, optional
        a list of coordinates that can be used for exporting data

    """

    def __init__(self, areas, qij, entrypoints=[0], regions=None, coords=None):
        super().__init__(areas, qij, regions, coords)
        self.set_entrypoints(entrypoints)

    def set_entrypoints(self, entrypoints):
        elist = []
        for el in entrypoints:
            if el in self.regions:
                elist.extend(self.regions[el])
            else:
                elist.append(el)

        self._entrypoints = set(elist)
        self._entryarea = sum([self.areas[i] for i in self.entrypoints])
        self._growthsections = [i for i in range(self.N) if i not in self.entrypoints]

    @property
    def entrypoints(self):
        return self._entrypoints
    
    @property
    def entryarea(self):
        return self._entryarea

    @property
    def growthsections(self):
        return self._growthsections

    def get_p0(self, entrypoints=[], flux=None):
        """Calculate the starting probability from a series
        of entrypoints.

        Parameters
        ----------
        entrypoints : list with entrypoints. Use default if empty

        """

        if entrypoints == []:
            entrypoints = self.entrypoints
        if flux is None:
            flux = {e:1.0 for e in entrypoints}

        p0 = np.zeros(self.N)
        for el,f in flux.items():
            p0[el] = f*self.areas[el]
        p0 /= sum(p0)
        return self.qij @ p0


    def solve_markov(self, sprobs, flux=None):
        Q = np.zeros((self.N,self.N))

        if not isinstance(sprobs, np.array):
            probs = np.ones(self.N)
            if isinstance(sprobs, dict):
                probs = np.ones(self.N)
                for k, v in sprobs.items():
                    probs[k] = v
            else:
                for i in self.growthsections:
                    probs[i] = sprobs
            sprobs = probs

        for i in range(self.N):
            Q[:,i] = (1-sprobs[i])*self.qij[:,i]

        ImQ = np.identity(self.N) - Q
        p0 = self.get_p0(self.entrypoints, flux)
        ptrans = np.linalg.solve(ImQ, p0)
        prob_outcome = sprobs * ptrans
        return prob_outcome, ptrans, self.entryarea*prob_outcome/self.areas


class Via(Feature):
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
        coords = np.zeros(N)
        coords[regions["wall"]] = AR/Nz*(0.5+np.arange(Nz))
        coords[regions["top"]] = 0
        coords[regions["bottom"]] = AR
        
        areas, qij = create_via(AR, Nz)
        Feature.__init__(self, areas, qij, regions=regions, coords=coords)

class Trench(Feature):
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
        Feature.__init__(self, areas, qij, regions=regions)



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



def create_trench(AR, Nseg):
    """Return the areas and view factor of a rectangular trench, where
    the vertical wall is divided into identical sections.

    Parameters
    ----------
    AR : float
        Aspect ratio, defined as the width to diameter ratio
    Nseg : int
        Number of vertical sections in the discretized wall

    Returns
    -------

    (numpy.array, numpy.array)
        Tuple with the areas (1D array), and view factors (2D array)

    """

    d0 = 1.0
    dh = d0*AR/Nseg

    Ai  = d0*np.ones(Nseg+1)
    Si = 2*dh*np.ones(Nseg)


    areas = np.zeros(Nseg+2)
    areas[0] = Ai[0]
    areas[-1] = Ai[-1]
    areas[1:-1] = Si
    areas = areas/areas[0]
    Ai = Ai/Ai[0]

    Fij = np.ones((Nseg+1,Nseg+1))
    for i in range(Nseg):
        for j in range(i+1, Nseg+1):
            Fij[j,i] = vf.strip_to_strip2(d0, dh*(j-i))
            Fij[i,j] = Fij[j,i]

    Qij = np.zeros((Nseg+2, Nseg+2))

    for i in range(1,Nseg):
        for j in range(i+1,Nseg+1):
            q0 = Ai[i+1]*(Fij[j-1,i]-Fij[j,i])-Ai[i]*(Fij[j-1,i-1]-Fij[j,i-1])
            Qij[j,i] = q0/areas[i]
            Qij[i,j] = q0/areas[j]

    Qij[1,0] = 1-Fij[1,0]
    for i in range(2,Nseg+1):
        Qij[i,0] = Fij[i-1,0]-Fij[i,0]
    Qij[-1,0] = 1-np.sum(Qij[:-1,0])

    for i in range(1,Nseg+1):
        Qij[i,-1] = Fij[i,-1]-Fij[i-1,-1]
    Qij[-2,-1] = 1-Fij[-2,-1]
    for  i in range(1,Nseg+1):
        Qij[0,i] = areas[0]/areas[i]*Qij[i,0]
        Qij[-1,i] = areas[-1]/areas[i]*Qij[i,-1]
    Qij[0,-1] = 1-np.sum(Qij[1:,-1])

    for i in range(1,Nseg+1):
        Qij[i,i] = 1-np.sum(Qij[:,i])

    return areas, Qij

def create_taperedvia(AR, dr, Nseg):
    """Return the areas and view factor of a rectangular trench, where
    the vertical wall is divided into identical sections.

    Parameters
    ----------
    AR : float
        Aspect ratio, defined as the width to top diameter ratio
    dr : Float
        Ratio between the bottom and the top diameter.
    Nseg : int
        Number of vertical segments in the discretized wall

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


def recper_to_cross(dh, W, n):
    if n == 0:
        Fwide = vf.rect_to_side(W, dh, 1.0)
        Fside = vf.rect_to_side(1.0, dh, W)
    else:
        ya = (0, W)
        yb = (0, W)
        xa = (n*dh, (n+1)*dh)
        zb = (0, 1)
        Fwide = vf.perp_rec_to_rec(xa, ya, yb, zb)
        ya = (0, 1)
        yb = (0, 1)
        xa = (n*dh, (n+1)*dh)
        zb = (0, W)
        Fside = vf.perp_rec_to_rec(xa, ya, yb, zb)
    return (2*W*Fwide + 2*Fside)/(2*W+2)


def create_rect(AR, W, Nseg):
    """Return the areas and view factor of a rectangular trench, where
    the vertical wall is divided into identical sections.

    Parameters
    ----------
    AR : float
        Aspect ratio, defined as the depth to height ratio
    W : float
        Normalized width, defined as the width to height ratio
    Nseg : int
        Number of vertical sections in the discretized wall

    Returns
    -------

    (numpy.array, numpy.array)
        Tuple with the areas (1D array), and view factors (2D array)

    """

    d0 = 1.0
    dh = d0*AR/Nseg

    Ai  = d0*W*np.ones(Nseg+1)
    Si = (2*dh*W+2*dh)*np.ones(Nseg)

    areas = np.zeros(Nseg+2)
    areas[0] = Ai[0]
    areas[-1] = Ai[-1]
    areas[1:-1] = Si

    F0 = []

    for n in range(Nseg+1):
        F0.append(recper_to_cross(dh, W, n))

    Qij = np.zeros((Nseg+2, Nseg+2))
    for i in range(1, Nseg+1):
        for j in range(i+1, Nseg+1):
            Qij[j,i] = F0[j-i-1]-F[j-i]
            Qij[i,j] = Qij[j,i]
        Qij[0,i] = F0[i]
        Qij[i,0] = areas[i]*F0[i]/areas[0]
        Qij[Nseg+1,i] = F0[Nseg-i]
        Qij[i, Nseg+1] = areas[i]*F0[Nseg-i]/areas[0]
        Qij[i,i] = 1 - np.sum(Qij[:,i])
    Qij[-1, 0] = np.sum(Qij[1:-1,0])
    Qij[0,-1] - Qij[-1,0]

    return areas/areas[0], Qij


if __name__ == "__main__":
    a, q = create_via(10,20)
    print(np.sum(q, axis=0), np.sum(q, axis=1),  np.sum(q))
    a, q = create_trench(10,20)
    print(np.sum(q, axis=0), np.sum(q, axis=1), np.sum(q))
    a, q = create_taperedvia(20,0.5,40)
    print(np.sum(q, axis=0), np.sum(q, axis=1), np.sum(q))
