#Copyright 2013 Argonne UChicago LLC
#This file is part of Machball
"""
Solve the ballistic transport inside nanostructures using the
Markov chain approach

"""

import numpy as np

def solve_constant(st, beta0):
    """
    Solve the steady state reactive transport for a constant reaction probability.

    Parameters
    ----------
        st : Structure
            nanostructure in which the ballistic transport takes place
        beta0 : float
            sticking probability.

    Returns
    -------
        (nd.array, float)
            Tuple of 1D array representing the outcome probabilities and the
            effective sticking probability

    """

    sprobs = beta0*np.ones(st.N)
    for i in st.entrypoints:
        sprobs[i] = 1.0
    probs = ballistic_markov(sprobs, st, st.entrypoints, transients=False)
    p = 1-sum([probs[i] for i in st.entrypoints])
    growth = np.array([probs[i] for i in range(st.N) if i not in st.entrypoints])
    return growth, p


def solve_ideal(st, beta0, betarec=0, endcov=0.95, ep=0.05, wt=0.01,
        verbose=True, return_betaeff=False):
    r"""
    Solve the reactive transport inside a nanostructure for an ideal
    self-limited process

    The kinetic model is a first order irreversible Langmuir, where
    molecules react with an available surface sites with a reaction
    probability beta0. We can optionally define a coverage
    independent recombination probability betarec.

    The simulation uses normalized time using the impingement rate of gas phase
    molecule into a single site as normalization units:

    .. math::
        \frac{1}{t_0}= \frac{1}{4} s_0 v_{th} \frac{p}{k_BT} s_0

    Parameters
    ----------

    st : Structure
        structure to be modeled
    beta0 : float
        bare sticking probability for the self-limited process
    betarec : float, optional
        recombination probability
    endcov : float, optional
        termination condition, given by the final minimum surface coverage
    ep : float, optional
        maximum coverage variation allowed by the implicit method used to
        solve the evolution of surface coverage.
    wt : float, optional
        different in average surface coverage between write intervals
    verbose : bool, optional
        if `true`, outputs intermediate steps to stdout
    return_betaeff : bool, optional
        if `true`, returns also the effective sticking probability

    Returns
    -------

    Tuple(numpy.array)
        returns a tuple containing the normalized times, the surface coverage
        at every point of the feature and time, and, if `return_betaeff` is true,
        the effective sticking probability.

    """

    av = np.ones(st.N-1)
    dtl = []
    cl = []
    betaeffl = []
    t = 0
    endav = 1-endcov
    totalarea = np.sum(st.areas[1:])

    sav = np.sum(av*st.areas[1:])/totalarea

    while np.amax(av) > endav:
        betatot = beta0*av + betarec
        sprobs = np.ones(st.N)
        sprobs[1:] = betatot
        flux = ballistic_markov(sprobs, st, [0])
        sticking = 1-flux[0]
        persite = st.areas[0]*beta0*flux[1:]/st.areas[1:]
        maxp = np.amax(persite)
        dt = ep/maxp

        av = av/(1+persite*dt)
        t += dt
        nsav = np.sum(av*st.areas[1:])/totalarea
        if sav-nsav > wt:
            dtl.append(t)
            cl.append(1-av)
            betaeff = 1-flux[0]
            if verbose:
                print(t, nsav, betaeff)
            if return_betaeff:
                betaeffl.append(betaeff)
            sav = nsav
    dtl.append(t)
    cl.append(1-av)
    if return_betaeff:
        betaeff = 1-flux[0]
        betaeffl.append(betaeff)
        return np.array(dtl), np.array(cl), np.array(betaeffl)
    else:
        return np.array(dtl), np.array(cl)


def solve_ideal2(st, beta0, betarec=0, endmode="average", startav=None, endcov=0.95, ep=0.05, wt=0.01,
        verbose=True, return_betaeff=False):
    r"""
    Solve the reactive transport inside a nanostructure for an ideal
    self-limited process

    The kinetic model is a first order irreversible Langmuir, where
    molecules react with an available surface sites with a reaction
    probability beta0. We can optionally define a coverage
    independent recombination probability betarec.

    The simulation uses normalized time using the impingement rate of gas phase
    molecule into a single site as normalization units:

    .. math::
        \frac{1}{t_0}= \frac{1}{4} s_0 v_{th} \frac{p}{k_BT} s_0

    Parameters
    ----------

    st : Structure
        structure to be modeled
    beta0 : float
        bare sticking probability for the self-limited process
    betarec : float, optional
        recombination probability
    endcov : float, optional
        termination condition, given by the final minimum surface coverage
    ep : float, optional
        maximum coverage variation allowed by the implicit method used to
        solve the evolution of surface coverage.
    wt : float, optional
        different in average surface coverage between write intervals
    verbose : bool, optional
        if `true`, outputs intermediate steps to stdout
    return_betaeff : bool, optional
        if `true`, returns also the effective sticking probability

    Returns
    -------

    Tuple(numpy.array)
        returns a tuple containing the normalized times, the surface coverage
        at every point of the feature and time, and, if `return_betaeff` is true,
        the effective sticking probability.

    """

    surface_ids = [i for i in range(st.N) if i not in st.entrypoints]
    nsurf = len(surface_ids)

    if startav is None:
        av = np.ones(nsurf)
    else:
        av = startav
        assert(len(av)==nsurf)

    dtl = []
    cl = []
    betaeffl = []
    t = 0
    endav = 1-endcov
    s_areas = np.array([st.areas[i] for i in surface_ids])
    totalarea = np.sum(s_areas)
    entryarea = sum(st.areas[i] for i in st.entrypoints)

    sav = np.sum(av*s_areas)/totalarea
    nsav = np.sum(av*s_areas)/totalarea

    if endmode == "average":
        cond = lambda av, sav, t:  sav > endav
    elif endmode == "time":
        cond  = lambda av, sav, t : t < endcov
    else:
        cond = lambda av, sav, t: np.amax(av) > endav

    while cond(av, sav, t):
        betatot = beta0*av + betarec
        sprobs = np.ones(st.N)
        for i, k in enumerate(surface_ids):
            sprobs[k] = betatot[i]
        flux = ballistic_markov(sprobs, st, st.entrypoints)
        fluxout = sum(flux[i] for i in st.entrypoints)
        probs = np.array([flux[i] for i in surface_ids])
        betaeff = 1-fluxout
        persite = entryarea*beta0*probs/s_areas
        maxp = np.amax(persite)
        dt = ep/maxp

        av = av/(1+persite*dt)
        t += dt
        nsav = np.sum(av*s_areas)/totalarea
        if sav-nsav > wt:
            dtl.append(t)
            cl.append(1-av)
            if verbose:
                print(t, nsav, betaeff)
            if return_betaeff:
                betaeffl.append(betaeff)
            sav = nsav
    dtl.append(t)
    cl.append(1-av)
    if return_betaeff:
        betaeffl.append(betaeff)
        return np.array(dtl), np.array(cl), np.array(betaeffl)
    else:
        return np.array(dtl), np.array(cl)


def dose_ideal(st, tdose, beta0, betarec=0, startav=None, ep=0.05,
        verbose=True):

    r"""
    Return the surface coverage inside a structure after a dose of an ideal
    self-limited process

    The kinetic model is a first order irreversible Langmuir, where
    molecules react with an available surface sites with a reaction
    probability beta0. We can optionally define a coverage
    independent recombination probability betarec.

    The simulation uses normalized time using the impingement rate of gas phase
    molecule into a single site as normalization units:

    .. math::
        \frac{1}{t_0}= \frac{1}{4} s_0 v_{th} \frac{p}{k_BT} s_0

    Parameters
    ----------

    st : Structure
        structure to be modeled
    tdose : float
        normalized dose time
    beta0 : float
        bare sticking probability for the self-limited process
    betarec : float, optional
        recombination probability
    startav : numpy.array, optional
        fraction of available sites at the beginning of the dose
    ep : float, optional
        maximum coverage variation allowed by the implicit method used to
        solve the evolution of surface coverage.
    verbose : bool, optional
        if `true`, outputs intermediate steps to stdout

    Returns
    -------

    Tuple(numpy.array)
        returns a tuple containing the normalized times, the surface coverage
        at every point of the feature and time, and, if `return_betaeff` is true,
        the effective sticking probability.

    """


    surface_ids = [i for i in range(st.N) if i not in st.entrypoints]
    nsurf = len(surface_ids)

    if startav is None:
        av = np.ones(nsurf)
    else:
        av = startav
        assert(len(av)==nsurf)

    dtl = []
    cl = []
    betaeffl = []
    t = 0
    s_areas = np.array([st.areas[i] for i in surface_ids])
    totalarea = np.sum(s_areas)
    entryarea = sum(st.areas[i] for i in st.entrypoints)

    sav = np.sum(av*s_areas)/totalarea
    nsav = np.sum(av*s_areas)/totalarea

    while t < tdose:
        betatot = beta0*av + betarec
        sprobs = np.ones(st.N)
        for i, k in enumerate(surface_ids):
            sprobs[k] = betatot[i]
        flux = ballistic_markov(sprobs, st, st.entrypoints)
        fluxout = sum(flux[i] for i in st.entrypoints)
        probs = np.array([flux[i] for i in surface_ids])
        betaeff = 1-fluxout
        persite = entryarea*beta0*probs/s_areas
        maxp = np.amax(persite)
        dt = ep/maxp
        if t+dt > tdose:
            dt = tdose-t
        av = av/(1+persite*dt)
        t += dt
        nsav = np.sum(av*s_areas)/totalarea
    return 1-av



def evolve_slowsat(st, beta1, beta2, f2, betarec=0, endmode="average",
    endcov=0.95, ep=0.05, wt=0.01, verbose=True, return_betaeff=False):
    r"""
    Solve the reactive transport inside a nanostructure for a
    soft-saturating self-limited process with two reaction pathways.

    The kinetic model is a first order irreversible Langmuir with two
    types of surface sites, with the second pathway comprising
    a fraction `f2` of the surface sites. Each pathway is
    characterized by its own reaction probability `beta1` and
    `beta2`, and the evolution of its respective surface coverages
    is given by a first order irreversible Langmuir kinetics.
    We can optionally define a coverage
    independent recombination probability `betarec`.

    The simulation uses normalized time using the impingement rate of gas phase
    molecule into a single site as normalization units:

    .. math::
        \frac{1}{t_0}= \frac{1}{4} s_0 v_{th} \frac{p}{k_BT} s_0

    Parameters
    ----------

    st : Structure
        structure to be modeled
    beta1 : float
        reaction probability for the first self-limited process
    beta2 : float
        reaction probability for the second self-limited process
    f2 : float
        fraction of sites for the second pathway. It should be a number
        between 0 and 1.
    betarec : float, optional
        recombination probability
    endcov : float, optional
        termination condition, given by the final minimum surface coverage
    ep : float, optional
        maximum coverage variation allowed by the implicit method used to
        solve the evolution of surface coverage.
    wt : float, optional
        different in average surface coverage between write intervals
    verbose : bool, optional
        if `true`, outputs intermediate steps to stdout
    return_betaeff : bool, optional
        if `true`, returns also the effective sticking probability

    Returns
    -------

    Tuple(numpy.array)
        returns a tuple containing the normalized times, the surface coverage
        at every point of the feature and time for each of the process and,
        if `return_betaeff` is true, the effective sticking probability.

    """

    av1 = np.ones(st.N-1)
    av2 = np.ones(st.N-1)
    dtl = []
    cl1 = []
    cl2 = []
    betaeffl = []
    t = 0
    endav = 1-endcov
    totalarea = np.sum(st.areas[1:])
    sav1 = np.sum(av1*st.areas[1:])/totalarea
    sav2 = np.sum(av2*st.areas[1:])/totalarea
    sav = (1-f)*sav1 + f*sav2

    while np.amax(av1+av2) > endav:
        betaald1 = beta1*av1*(1-f2)
        betaald2 = beta2*av2*f2
        betatot = betaald1 + betaald2 + betarec
        sprobs = np.ones(st.N)
        sprobs[1:] = betatot
        flux = ballistic_markov(sprobs, st, [0])
        persite1 = (1-f)*st.areas[0]*beta1*flux[1:]/st.areas[1:]
        persite2 = (1-f)*st.areas[0]*beta2*flux[1:]/st.areas[1:]
        maxp1 = np.amax(persite1)
        maxp2 = np.amax(persite2)
        dt = ep/max(maxp1, maxp2)
        av1 = av1/(1+persite1*dt)
        av2 = av2/(1+persite2*dt)
        t += dt
        nsav1 = np.sum(av1*st.areas[1:])/totalarea
        nsav2 = np.sum(av2*st.areas[1:])/totalarea
        nsav = (1-f)*nsav1 + f*nsav2
        if sav-nsav > wt:
            dtl.append(t)
            cl1.append(1-av1)
            cl2.append(1-av2)
            betaeff = 1-flux[0]
            if verbose:
                print(t, nsav, betaeff)
            if return_betaeff:
                betaeffl.append(betaeff)
            sav = nsav

    dtl.append(t)
    cl1.append(1-av1)
    cl2.append(1-av2)
    if return_betaeff:
        betaeff = 1-flux[0]
        betaeffl.append(betaeff)
        return np.array(dtl), np.array(cl1), np.array(cl2), np.array(betaeffl)
    else:
        return np.array(dtl), np.array(cl1), np.array(cl2)


def ballistic_markov(sprobs, st, entrypoints, transients=True):
    """
    Solve the ballistic transport inside a structure.

    The function solves the process as a Markov chain, where surfaces
    and entry points contains absorbing states (in the Markov chain
    sense) that terminates the stochastic process if a particle reacts
    or reaches one of the entry points (i.e. leaves the feature).

    The initial condition
    is calculated from the view factors of the entry points, with the
    probability of entering through a particular entry point being
    proportional to its surface area. This is equivalent to considering
    that the pressure is the same at all entry points of the structure.

    The probability of a particle reaching the interior of the feature
    is calculated from the view factor of each item in entrypoints.
    The probability is proportional to the area of each region in entrypoint.
    This is equivalent to considering that the pressure is the same at all
    entry points of the structure.

    In order to properly model the reactive transport process, the sticking
    probabilities of each of the entry points should be set to 1, to
    indicate no reentry condition.

    The outcome of the Markov chain process is an array of probabilities,
    indicating the probability that a particle reacts with that section
    of the structure or, in the case of the entry points, leaves the
    structure through that specific section.
    Consequently, one minus the sum of all the outcome probabilities for all
    the entrypoints defines the effective sticking probability, the probability
    that an incoming particle reacts somewhere inside the Structure.


    Parameters
    ----------
    sprobs : numpy.array
        Sticking probabilities at each point of the feature. Entry points should
        have sticking probabilities equal to one.
    st : Structure
        A structure variable, with well-defined areas and view factors for each
        of its elements.
    entrypoints : [int]
        A list of indices indicating the entry points of the structure
    transient : bool, optional
        If true, returns the total flux reaching it section normalized to
        the total area of the entry points. Otherwise, returns the outcome
        probabilities.

    Returns
    -------
    numpy.array
        if `transient` is true returns the total flux reaching each section
        normalized to the total area of the entry points. Otherwise, returns
        the outcome probabilities.

    """

    Q = np.zeros((st.N,st.N))
    for i in range(st.N):
        Q[:,i] = (1-sprobs[i])*st.qij[:,i]
    ImQ = np.identity(st.N) - Q
    p0 = st.get_p0(entrypoints)
    ptrans = np.linalg.solve(ImQ, p0)
    if transients:
        outcome = ptrans
    else:
        outcome = sprobs * ptrans
    return outcome

if __name__ == "__main__":

    from structure import Via
    via = Via(10,20)
    solve_constant(via, 0.01)
