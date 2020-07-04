#Copyright 2013 Argonne UChicago LLC
#This file is part of Machball

from machball.ballistic import Via, save_structure, read_structure

st = Via(20, 40)

save_structure("via.pickle", st, mode="pickle")
st2 = read_structure("via.pickle", mode="pickle")

save_structure("via.dat", st, mode="numpy")
st3 = read_structure("via.dat", mode="numpy")

save_structure("via.dat", st, mode="numpy", areafile="via_areas.dat")
st4 = read_structure("via.dat", mode="numpy", areafile="via_areas.dat")

print(st.qij[0,1])
print(st2.qij[0,1])
print(st3.qij[0,1])
print(st4.qij[0,1])
