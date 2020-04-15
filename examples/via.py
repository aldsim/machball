#Copyright 2013 Argonne UChicago LLC
#This file is part of Machball

from machball import ALDIdeal
from machball.ballistic import Via

import numpy as np

ald = ALDIdeal(1e-2, 100, 473, 10, 10e-20, betarec=0)
x, y = ald.saturation_flat()
#print(x[-1])
arl = np.array([1, 2, 3, 5, 10, 20, 30, 50, 100, 200, 300, 500])
times = []
for AR in arl:
    print(AR)
    N = 2*AR
    st = Via(AR, N)
    x2, y2 = ald.saturation_ballistic(st, endcov=y[-1], verbose=False)
    times.append(x2[-1]/x[-1])
times = np.array(times)

for a, t in zip(arl, times):
    print(a, t)
