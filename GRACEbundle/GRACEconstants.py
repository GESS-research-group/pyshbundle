# -*- coding: utf-8 -*-
"""
Created on Sat May  9 18:49:45 2020

@author: bv18488
"""
clight = 2.99792458e8	      # speed of light [m/s]
G      = 6.67259e-11         # gravitational constant [m^3 /(kg s^2)]
au     = 149.597870691e9     # astronomical unit [m]

# GRS80 defining constants:
ae     = 6378137             # semi-major axis of ellipsoid [m]
GM     = 3.986005e14         # geocentric grav. constant [m^3 / s^2]
J2     = 1.08263e-3          # earth's dyn. form factor (= -C20 unnormalized)
Omega  = 7.292115e-5         # mean ang. velocity [rad/s]

# GRS80 derived constants:
flat = 1/298.257222101     # flattening
J4   = -0.237091222e-5     # -C40 unnormalized
J6   =  0.608347e-8        # -C60 unnormalized
J8   = -0.1427e-10         # -C80 unnormalized

