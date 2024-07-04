#Copyright 2013 Argonne UChicago LLC
#This file is part of Machball

"""
A collection of view factors for different geometries

"""
import numpy as np

def cylinder_base_to_wall(R1, H1):
    """Compute the view factor of a cylinder base to its wall
    """
    H = H1/(2*R1)
    return 2*H*(np.sqrt(1+H*H)-H)

def cylinder_base_to_base(R1, H1):
    """Compute the view factor of a cylinder base to the opposite base
    """
    return 1-cylinder_base_to_wall(R1,H1)

def cylinder_wall_to_itself(R1, H1):
    """Compute the view factor of a cylinder wall to itself
    """
    r = R1/H1
    rh = (np.sqrt(4*r*r+1)-1)/r
    return 1-0.5*rh

def cylinder_wall_to_disk(Rc, Rd, h1, h2):
    R = Rc/Rd
    H1 = h1/Rd
    H2 = h2/Rd
    r2 = R*R
    X1 = 1+H1*H1 + r2
    X2 = 1+H2*H2 + r2
    return (X1-X2-np.sqrt(X1*X1-4*r2)+np.sqrt(X2*X2-4*r2))/(4*R*(H2-H1))

def cylinder_section_to_section(R, dz, n):
    snear = cylinder_wall_to_disk(R, R, dz*(n-1), dz*n)
    sfar =  cylinder_wall_to_disk(R, R, dz*n, dz*(n+1))
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

def square_to_square(W, H):
    w = W/H
    w2 = w*w
    x = np.sqrt(1+w2)
    y = x*np.arctan(w/x)-np.arctan(w)
    return (np.log(x**4/(1+2*w2)) + 4*w*y)/(np.pi*w2)


def rect_to_side(W, L, H):
    h = H/L
    w = W/L
    h2 = h*h
    w2 = w*w
    a = (1+h2)*(1+w2)/(1+h2+w2)
    b = w2*(1+h2+w2)/((1+w2)*(h2+w2))
    c = h2*(1+h2+w2)/((1+h2)*(h2+w2))
    d = np.sqrt(h2+w2)
    F = h*np.arctan(1/h) + w*np.arctan(1/w) - d*np.arctan(1/d)
    F += 0.25*np.log(a*b**w2*c**h2)
    F /= np.pi*w
    return F

def par_rec_to_rec(xa, ya, xb, yb, z):
    F = 0
    for i in range(2):
        for j in range(2):
            for k in range(2):
                for l in range(2):
                    x = xa[i]
                    y = ya[j]
                    u = xb[k]
                    v = yb[l]
                    sgn = 1 if (i+j+k+l) %2 == 0 else -1
                    F += sgn*Bpar(x, y, u, v, z)
    return F/(2*np.pi*(xa[1]-xa[0])*(ya[1]-ya[0]))

def Bpar(x, y, u, v, z):
    dx = x-u
    dy = y-v
    dx2 = dx*dx
    dy2 = dy*dy
    z2 = z*z
    p = np.sqrt(dx2 + z2)
    q = np.sqrt(dy2 + z2)
    return dy*p*np.arctan(dy/p) + dx*q*np.arctan(dx/q) - z2/2*np.log(dx2 + dy2 + z2)

def perp_rec_to_rec(xa, ya, yb, zb):
    F = 0
    for i in range(2):
        for j in range(2):
            for k in range(2):
                for l in range(2):
                    sgn = 1 if (i+j+k+l) %2 == 0 else -1
                    F += sgn*Bperp(xa[i], ya[j], yb[k], zb[l])
    return F/(2*np.pi*(xa[1]-xa[0])*(ya[1]-ya[0]))
  

def Bperp(x, y, y2, z):
    C = np.sqrt(x*x + z*z)
    D = (y-y2)/C
    C2 = C*C
    D2 = D*D
    return (y-y2)*C*np.arctan(D)-C2/4*(1-D2)*np.log(C2*(1+D2))
