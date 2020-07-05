Getting started
===============

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
