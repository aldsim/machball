#Copyright 2013 Argonne UChicago LLC
#This file is part of Machball
r"""
This module contains the implementations of different ALD processes.
For now, it is limited to an ideal self-limited kinetics with an
optional surface recombination pathway.

"""

from .utils import vth, kB
from .ballistic import solve_ideal, dose_ideal

import numpy as np


class ALDIdeal:
    r"""Model for an ideal self-limited process

    ALDIdeal implements a self-limited process as an irreversible
    first order langmuir kinetics. The evolution of the surface
    coverage with time is given by:

    .. math::

        \frac{d\Theta}{dt} = s_0 \frac{1}{4}v_{th}\frac{p_0}{k_BT}\beta_0
        (1-\Theta)

    Parameters
    ----------

    beta0 : float
        Bare sticking probability
    MM : float
        Molecular mass in amu
    T : float
        Process temperature in K
    p0 : float
        Precursor pressure in Pa
    s0 : float
        Area of a surface site in sq. meters
    betarec : float, optional
        Recombination probability

    """

    def __init__(self, beta0=1e-2, MM=150, T=473,
            p0=1e2, s0=10e-19, betarec=0):

        self._beta0 = beta0
        self._MM = MM
        self._T = T
        self._p0 = p0
        self._s0 = s0
        self.betarec = betarec

        self._update()

    def _update(self):

        self.vth = vth(self.T, self.MM)
        self.n0 = self.p0/(kB*self.T)
        self.t0 = 1.0/(0.25*self.vth*self.n0*self.s0)

    @property
    def beta0(self):
        return self._beta0

    @beta0.setter
    def beta0(self, beta0):
        self._beta0 = beta0

    @property
    def MM(self):
        return self._MM

    @MM.setter
    def MM(self, MM):
        self._MM = MM
        self._update()

    @property
    def T(self):
        return self._T

    @T.setter
    def T(self, T):
        self._T = T
        self._update()

    @property
    def s0(self):
        return self._s0

    @s0.setter
    def s0(self, s0):
        self._s0 = s0
        self._update()

    @property
    def p0(self):
        return self._p0

    @p0.setter
    def p0(self, p0):
        self._p0 = p0
        self._update()

    def dose(self, t, st=None):
        if st is None:
            return self.coverage_flat(t)
        else:
            return dose_ideal(st, t/self.t0, self.beta0,
                self.betarec)

    def saturation(self, st=None, endcov=0.95, verbose=True):
        if st is None:
            return self.saturation_flat()
        else:
            return self.saturation_ballistic(st, endcov, verbose)

    def coverage_flat(self, t):
        """Return the surface coverage on a flat surface

        Calculate the surface coverage on a flat surface for
        a dose time of duration `t`

        Parameters
        ----------

        t : float
            Dose time in seconds

        Returns
        -------

        float
            Fractional surface coverage

        """
        return 1-np.exp(-self.beta0*t/self.t0)

    def saturation_flat(self):
        """Return the saturation curve on a flat surface

        Returns
        -------

        Tuple
            Tuple of numpy arrays, containing the dose time in seconds
            and the fractional surface coverage at the surface

        """

        t = np.arange(0.1,6,0.1)*self.t0/self.beta0
        c = self.coverage_flat(t)
        return t, c


    def saturation_ballistic(self, st, endcov=0.95, verbose=True):
        """Calculate the evolution of the coverage profile inside a structure

        Solve the ballistic transport of precursor inside a structure `st`,
        such as a trench or via.

        Parameters
        ----------
        st : Structure
            The structure model to be modeled
        covend : float, optional
            The final surface coverage at the bottom of the structure

        Returns
        -------
        Tuple
            A pair of numpy arrays, representing the dose times (in seconds),
            and a 2D array containing the fractional surface coverage for
            each element of the structure for every time step.

        """

        times, cov = solve_ideal(st, self.beta0, self.betarec, endcov=endcov,
            verbose=verbose)
        times = self.t0*times
        return times, cov
