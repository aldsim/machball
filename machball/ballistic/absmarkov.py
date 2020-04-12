# Copyright Angel Yanguas-Gil, Jeffrey W. Elam, 2013

"""
absmarkov is a general implementation of an absorbing markov chain

An absorbing Markov chain process is defined by the existence of Nt
transient states and Na absorbing states. Transient states are defined
as states from which the system can escape. Absorbing states are end-game
states from which the system cannot escape.

Transitions between these states are controlled by transition probabilities.
A Nt x Nt matrix is required to describe the transition probabilities
between transient states. Likewise, a NaxNs matrix is required to describe
the transition from transient to absorbing states.

Since these matrices are sparse, absmarkov uses a dictionary to codify
the transition probabilities using a (i,j) tuple to codify
a transition from i to j: p(i->j) = tprob[(i,j)].

"""

import numpy as np

class AbsMarkovChain:

    def __init__(self, Nt, Na, tprob, aprob):
        self.Nt = Nt
        self.Na = Na
        self.tprob = tprob
        self.aprob = aprob
        self.ImQ = np.identity(self.Nt) - self.tprob

    def __call__(self, ptrans, pabs=None):
        print(self.ImQ, ptrans)
        ptrans = np.linalg.solve(self.ImQ, ptrans)
        pout = self.aprob @ ptrans
        if pabs is None:
            return pout
        else:
            return pout + pabs


def dicttomatrix(trans, n1, n2):
    """Take a sparse matrix represented as a dictionary and return the
    corresponding matrix.

    trans is a dictionary, in which non-zero ij elements of the matrix are
    stored using the (i,j) key, so that trans[(i,j)] = pij is the
    probability of the system undergoing the i->j transition.

    n1, n2 are the dimensions of the matrix, n1xn2
    """

    m = np.zeros((n1 ,n2))
    for k, v in trans.items():
        i, f = k
        m[i,f] = v
    return m

def dicttovector(spv, n):
    """Take a sparse vector represented as a dictionary and return the
    corresponding vector.

    spv is a dictionary so that v[k] = spv[k] for all the keys

    n is the vector size.
    """

    m = np.zeros(n)
    for k, v in spv.items():
        m[k] = v
    return m

def make_start(probstart, nt):
    """Take a dictionary of starting probabilities and return the
    corresponding 1xnt matrix

    """
    m = dicttovector(probstart, nt)
    return m


def make_transient(probdict, n):
    """Create the transition probability matrix between transient
    states. Take a dictionary of transition probabilities where
    keys are tuples (i,f) indicating the transition probability from
    the initial to the final state. n is the number of transient
    states in the system.

    Return a matrix object
    """

    return dicttomatrix(probdict, n, n)


def make_absorbing(probdict, nt, nabs):
    """Create the transition probability matrix between transient
    and absorbing states. Take a dictionary of transition
    probabilities where the keys are (it, fa) pairs connecting
    the initial transient state with the final absorbing state.
    nt and nabs are the number of transition and absorbing states
    in the system

    Return the corresponding nt x nabs matrix.

    """

    return dicttomatrix(probdict, nt, nabs)


def create_absmarkov(trprob, abprob, Nt, Na):
    """Return an absorbing Markov chain object from the corresponding
    transition probabilities.

    trprob is a dictionary codifying the transition probability from
        i to j codified as trprob[(i,j)]
    abprob is a dictionary codifying the transition probability from
        i to a codified as abprob[(i,a)]
    Nt is the number of transient states
    Na is the number of absorbing states

    Both transient and absorbing states are labeled from 0 to Nt and Na,
    respectively.

    """

    trmat = make_transient(trprob, Nt)
    abmat = make_absorbing(abprob, Nt, Na)
    return AbsMarkovChain(Nt, Na, trmat, abmat)


if __name__ == '__main__':

    print("Example: gambler's ruin with p=q=0.5")
    d = {(0,1):0.5, (1,0):0.5}
    absdict = {(0,0):0.5, (1,1):0.5}
    sdict = {0:1}
    s = make_start(sdict,2)
    markov = create_absmarkov(d, absdict, 2, 2)
    prob = markov(s)
    print("Outcome:", prob)
