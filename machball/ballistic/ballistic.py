"""
Solve the ballistic transport inside nanostructures using the
Markov chain approach

"""

import numpy as np

def solve_constant(st, beta0):
    """
    Solve the steady state reactive transport for a constant reaction probability beta0

    st is a Structure. The algorithm assumes that the opening is located at index 0
    beta0 is either a constant or an array of size N-1, where N is the number of
    discrete regions in the structure.
    """

    sprobs = np.ones(st.N)
    sprobs[1:] = beta0
    return ballistic_markov(sprobs, st, [0])


def solve_ideal(st, beta0, betarec=0, endcov=0.95, ep=0.05):
    """
    Solve the time-dependent ALD growth inside a nanostructure st with two
    surface components.

    The kinetic model is a first order irreversible Langmuir, where
    molecules react with an available surface sites with a reaction
    probability beta0. We can optionally define a coverage
    independent recombination probability betarec.

    The simulation uses normalized time using the impingement rate of gas phase
    molecule into a single site as normalization units: 1/4 vth n s_0, with n
    the number of molecules per unit volume, related to the precursor
    partial pressure p through: p = n kB T.
    """
    av = np.ones(st.N-1)
    dtl = [0]
    cl = [1-av]
    t = 0
    endav = 1-endcov
    while np.amax(av) > endav:
        betaald = beta0*av
        betatot = betaald + betarec
        pald = betaald/betatot
        prec = betarec/betatot
        sprobs = np.ones(st.N)
        sprobs[1:] = betatot
        flux = ballistic_markov(sprobs, st, [0])
        aldflux = pald*flux[1:]
        recflux = prec*flux[1:]
        persite = st.areas[0]*aldflux/st.areas[1:]
        maxp = np.amax(persite)
        dt = ep/maxp
        t += dt
        print(t, 1-np.mean(av), np.amax(av))
        av = av/(1+persite*dt)
        dtl.append(t)
        cl.append(1-av)
    return dtl, cl, aldflux, recflux


def evolve_slowsat(st, beta1, beta2, f2, betarec=0, endcov=0.95, ep=0.05):
    """
    Solve the time-dependent ALD growth inside a nanostructure st with two surface
    components.

    It consider two parallel reaction channels, with the second one comprising a
    fraction f of the total number of sites.

    In both cases, the kinetic model is a first order irreversible Langmuir kinetics,
    where molecules react with an available surface sites with reaction probabilities
    beta1 and beta2 for each zone. We can optionally define a coverage
    independent recombination probability betarec.

    The simulation uses normalized time using the impingement rate of gas phase
    molecule into a single site as normalization units: 1/4 vth n s_0, with n
    the number of molecules per unit volume, related to the precursor
    partial pressure p through: p = n kB T.

    The use of a second reaction pathway provides a simple way of modeling systems
    that are soft-saturating.
    """

    av1 = np.ones(st.N-1)
    av2 = np.ones(st.N-1)
    dtl = [0]
    cl1 = [1-av1]
    cl2 = [1-av2]
    beffl = []
    t = 0
    endav = 1-endcov
    while np.amax(av1) > endav or np.amax(av2) > endav:
        betaald1 = beta1*av1*(1-f2)
        betaald2 = beta2*av2*f2
        betatot = betaald1 + betaald2 + betarec
        pald1 = betaald1/betatot
        pald2 = betaald2/betatot
        prec = betarec/betatot
        sprobs = np.ones(st.N)
        sprobs[1:] = betatot
        flux = ballistic_markov(sprobs, st, [0])
        beffl.append(1-flux[0])
        aldflux1 = pald1*flux[1:]
        aldflux2 = pald2*flux[1:]
        recflux = prec*flux[1:]
        persite1 = (1-f2)*st.areas[0]*aldflux1/st.areas[1:]
        persite2 = f2*st.areas[0]*aldflux2/st.areas[1:]
        maxp1 = np.amax(persite1)
        maxp2 = np.amax(persite2)
        dt = ep/max(maxp1, maxp2)
        t += dt
        print(t, 1-np.mean(av1), np.amax(av1), 1-np.mean(av2), np.amax(av2), beffl[-1])
        av1 = av1/(1+persite1*dt)
        av2 = av2/(1+persite2*dt)
        dtl.append(t)
        cl1.append(1-av1)
        cl2.append(1-av2)
    return dtl, cl1, cl2, beffl


def ballistic_markov(sprobs, st, entrypoints):
    """
    Solve the ballistic transport using a Markov chain method. It takes
    a structure and a reaction probability array. The initial condition
    is calculated assuming that particles enter the feature through a subset
    of regions defined by entrypoints.

    The probability of a particle reaching the interior of the feature
    is calculated from the view factor of each item in entrypoints.
    The probability is proportional to the area of each region in entrypoint.
    This is equivalent to considering that the pressure is the same at all
    entry points of the structure.

    In order to properly model a reactive transport process, sprobs must
    be set to one for all entry elements. However, this function does not
    check for that and instead it is responsibility of the caller to check
    that sprobs is properly defined.

    The outcome of the process
    is the probability that a molecule reacts in each region of
    the Structure or, in the case of an entry point, that it leaves the feature.

    Consequently, one minus the sum of all the outcome probabilities for all
    the entrypoints defines the effective sticking probability, the probability
    that an incoming particle reacts somewhere inside the Structure.
    """

    Q = np.zeros((st.N,st.N))
    for i in range(st.N):
        Q[:,i] = (1-sprobs[i])*st.qij[:,i]
    ImQ = np.identity(st.N) - Q
    p0 = st.qij @ st.get_p0(entrypoints)
    ptrans = np.linalg.solve(ImQ, p0)
    outcome = sprobs * ptrans
    return outcome

if __name__ == "__main__":

    from structure import Via
    via = Via(10,20)
    solve_constant(via, 0.01)
