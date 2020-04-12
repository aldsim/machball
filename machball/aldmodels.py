#Copyright 2013 Argonne UChicago LLC
#This file is part of MachBall
"""
Implements two models of self-limited processes
"""

from .utils import vth, kB
from .ballistic import solve_ideal

class ALDIdeal:
    """Model for an ideal self-limited process

    ALDIdeal implements a self-limited process as an irreversible
    first order langmuir kinetics. The evolution of the surface
    coverage with time is given by:

    ..math::
        \frac{d\Theta}{dt} = s_0 \frac{1}{4}v_{th}\frac{p_0}{k_BT}\beta_0
        (1-\Theta)

    """

    def __init__(self, beta0, MM, T, p0, s0, betarec=0):
        """
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

        self.beta0 = beta0
        self.MM = MM
        self.T = T
        self.p0 = p0
        self.s0 = s0
        self.betarec = betarec

        self._update()


    def _update(self):

        self.vth = vth(self.T, self.MM)
        self.n0 = self.p0/(kB*self.T)
        self.t0 = 1.0/(0.25*self.vth*self.n0*self.s0)

    def coverage_flat(self, t):
        return 1-np.exp(-self.beta0*t/self.t0)

    def saturation_flat(self):
        t = np.arange(100)*0.05*self.t0
        c = self.coverage_flat(t)
        return t, c

    def saturation_ballistic(self, st, endcov=0.95):
        """Calculates the evolution of the coverage profile inside a structure

        Solves the ballistic transport of precursor inside a nanostructured
        substrate or a feature, such as a trench or via.

        Parameters
        ----------
        st : Structure
            The structure model to be modeled
        covend : float, optional
            The final surface coverage at the bottom of the structure
        """

        self._update()
        times, cov = solve_ideal(self.beta0, self.betarec, endcov)
        times = times*t0
        return times, cov

# TODO: ALDSoftSat
