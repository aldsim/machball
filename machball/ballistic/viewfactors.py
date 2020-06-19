#Copyright 2013 Argonne UChicago LLC
#This file is part of Machball

"""
A collection of view factors for different geometries

"""
import numpy as np

def base_to_cylinderwall(R1, H1):
    H = H1/(2*R1)
    return 2*H*(np.sqrt(1+H*H)-H)

def base_to_base(R1, H1):
    return 1-base_to_cylinderwall(R1,H1)

def cylinderwall_to_itself(R1, H1):
    r = R1/H1
    rh = (np.sqrt(4*r*r+1)-1)/r
    return 1-0.5*rh

def cylinderwall_to_disk(Rc, Rd, h1, h2):
    R = Rc/Rd
    H1 = h1/Rd
    H2 = h2/Rd
    r2 = R*R
    X1 = 1+H1*H1 + r2
    X2 = 1+H2*H2 + r2
    return (X1-X2-np.sqrt(X1*X1-4*r2)+np.sqrt(X2*X2-4*r2))/(4*R*(H2-H1))

def cylindersection_to_section(R, dz, n):
    snear = cylinderwall_to_disk(R, R, dz*(n-1), dz*n)
    sfar =  cylinderwall_to_disk(R, R, dz*n, dz*(n+1))
    return snear-sfar

def strip_to_strip2(d, h):
    h = h/d
    return 1.0/(h+np.sqrt(1+h*h))

def strip_to_wall(dz, n):
    a = 1 + dz*dz*(n+1)**2
    b = 1 + dz*dz*n**2
    return dz - dz*dz*(1+2*n)/(np.sqrt(a)+np.sqrt(b))

def strip_to_strip(dz, n):
    if n > 0:
        return wall_to_strip(dz,n)-wall_to_strip(dz,n+1)
    else:
        return 1-2*wall_to_strip(dz, 1)

def wall_to_strip(dz, n):
    return 2*dz*strip_to_wall(dz, n)

def disk_to_disk(r, h):
    r = r/h
    r2 = 4*r*r
    return r2/(1+np.sqrt(r2+1))**2

def disk1_to_disk2(r1, r2, h):
    r22 = r2*r2
    r11 = r2*r1
    d = 4*r22/r11
    x = 1 + (h*h + r22)/(r11)
    return d/(x + np.sqrt(x-d))

def disk_to_shrinkcone(r1, h, tana):
    xi = h*tana/r1
    d = 4*(1-xi)**2
    x = 1  + h*h/(r1*r1) + (1-xi)**2
    return 0.5*d/(x+np.sqrt(x*x - d))
