# What is Machball?

Machball models the reactive transport
inside high aspect ratio and nanostructured materials for
self-limited processes such as atomic layer deposition (ALD)
and atomic layer etching (ALE).

# Documentation

Machball's documentation can be found under docs and in readthedocs.

# Install instructions

Using pip:
```
pip install machball
```

From github:
```
git clone https://github.com/aldsim/machball.git
cd machball/
pip install -e .
```

# Quickstart

You can model ALD inside a high aspect ratio feature with just
a few lines of code:
```
from machball import ALDIdeal
from machball.ballistic import Via
ald = ALDIdeal(1e-2, 100, 473, 10, 10e-20, betarec=0)
st = Via(50, 100) # Aspect ratio, and number of segments
dose_times, coverages = ald.saturation_ballistic(st)
```

# Authors

Machball was developed at Argonne National Laboratory. Currently
the team comprises:

* Angel Yanguas-Gil, <ayg@anl.gov>, Lead and founder
* Jeffrey W Elam

# Citing

If you are referencing MachBall in a publication, please cite
the following paper:

* Ballistic transport model:

    A. Yanguas-Gil and J. W. Elam, **A Markov chain approach to
    simulate Atomic Layer Deposition chemistry and transport inside
    nanostructured substrates**, Theoretical Chemistry Accounts
    133, Article number: 1465 (2014). http://dx.doi.org/10.1007/s00214-014-1465-x


# Acknowledgements

Machball development was partially funded through Argonne's
Laboratory Directed Research and Development program.

# Copyright and license

Copyright (2013), UChicago Argonne, LLC

Machball is distributed under the terms of the BSD License. A
copy of the license can be found [here](https://github.com/aldsim/machball/blob/master/LICENSE)

Argonne Software Number: SF-13-072
