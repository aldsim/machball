========
Tutorial
========

Machball computes the evolution in surface coverage inside
nanostructures for self-limited processes such as atomic
layer deposition (ALD) or atomic layer etching (ALE).

This tutorial assumes that the user has a working `Python`
environment and that Machball is already installed in the
system.


Ballistic transport
===================

In order to compute the ballistic transport inside a nanostructure,
Machball needs information about the view factors connecting different
points of a structure. Machball codifies this information in its
`Structure` class.

.. autoclass:: machball.ballistic.Structure

`Structure` implements a structure discretized in a finite number
of sections. A `Structure` object requires three key pieces of information:
the view factors `qij` codified as a 2D `numpy` array, the area
of the different elements, and a list of the sections representing
the opening of the feature.

Specific structures
-------------------

Machball implements a few examples of the most commonly used nanostructures,
such as circular vias, infinite rectangular trenches, and tapered circular
vias.

.. autoclass:: machball.ballistic.Via

.. autoclass:: machball.ballistic.Trench

These classes automatically calculate the view factors for the user.
For instance, if we want to work with a 100 aspect ratio circular
via discretized in 200 equally spaced segments, we just use::

    from machball.ballistic import Via
    st = Via(100, 200)


Saving and loading structures
-----------------------------

It is possible to save and load the structures from file. Machball currently
supports two formats: a `pickle` format using Python's `pickle` module and
a txt `numpy` array format:

.. autofunction:: machball.ballistic.save_structure

.. autofunction:: machball.ballistic.read_structure

We can use these functions as follows::

    save_structure("via.pickle", st, mode="pickle")
    st2 = read_structure("via.pickle", mode="pickle")

    save_structure("via.dat", st, mode="numpy")
    st3 = read_structure("via.dat", mode="numpy")

    save_structure("via.dat", st, mode="numpy", areafile="via_areas.dat")
    st4 = read_structure("via.dat", mode="numpy", areafile="via_areas.dat")

The key difference between the pickle and numpy modes is that when saving
and loading numpy arrays we lose all the additional metadata. This means
that we have to manually set the entrypoints unless it uses the default
value (first section in the array is the entrypoint).


Defining an ALD process
=======================

Machball implements an ideal self-limited process through
its class ``ALDIdeal``, which models the self-limited
adsorption of a precursor molecule as a first order
irreversible Langmuir kinetics.

The first step is to import the class::

    from machball import ALDIdeal

A self limited process is then defined as::

    ald = ALDIdeal(beta0=1e-2, MM=100, T=473,
        p0=100, s0=10e-20, betarec=0)

The parameters that we need to pass to define a self-limited
process includes the reaction probability ``beta0``,
the molecular mass in atomic mass units ``MM``, the process temperature
``T`` in Kelvin, the precursor pressure ``p0`` in Pa, the area of
an absorption site ``s0`` in square meters, and finally an optional
parameter codifying a recombination probability ``betarec`` (zero
by default).

Ballistic transport
-------------------

Once it has been installed, Machball is trivial to use. Here is
a simple snippet to model the transport inside a circular via::

    from machball import ALDIdeal
    from machball.ballistic import Via
    ald = ALDIdeal(1e-2, 100, 473, 10, 10e-20, betarec=0)
    st = Via(50, 100) # Aspect ratio, and number of segments
    dose_times, coverages = ald.saturation_ballistic(st)

The variables `dose_times` and `coverages` are 1D array with dose
times and a 2D array with the fractional coverage profile inside
the via for each of the times.
