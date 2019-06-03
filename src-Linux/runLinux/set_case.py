# FIXME 'srcPy' should be 'solstice_py' or something like that
from srcPy.gen_YAML import gen_YAML
import numpy as N
import os

# the sun
# =========

# S1. DNI
DNI=1000 # W/m2

# S2. sunshape
sunshape='pillbox' # or 'buie'
sunsize=0.2664 # the half angle of the pillbox sunshape distribution, in degree 
               # or CSR value of Buie sunshape

# S3. sun position
# e.g. summer solstice, solar noon
azimuth=270. # from East to North, deg
elevation =78.  # 0 is horizontal, deg

# S4. number of rays for the ray-tracing simulation
num_rays=2000000

# the field
# ==========

# F1.Layout
layout=N.loadtxt('./demo_layout.csv', delimiter=',', skiprows=2)
hst_pos=layout[:,:3]
hst_foc=layout[:,3] # F2.5
hst_aims=layout[:,4:] # F4.

# F2. Heliostat
hst_w=10. # m
hst_h=10. # m
rho_refl=0.95 # mirror reflectivity
slope_error=2.e-3 # radians

# the receiver
# ============

# R1. shape
receiver='flat'

# R2. Size
rec_w=8. # width, m
rec_h=6. # height, m

# R3. tilt angle
tilt=0. # deg

# R4. position
loc_x=0. # m
loc_y=0. # m
loc_z=62.# m

# R5. Abosrptivity
rec_abs=0.9
# receiver mesh, for binning the flux distribution
rec_mesh=100

# ===============================================================
# the lines below will automatically create the ray-tracing scene
# based on the settings above

# the folder where files will be saved is controlled by the 'run' script.
casefolder=os.environ['CASEFOLDER']

if not os.path.exists(casefolder):
    os.makedirs(casefolder) 
rec_param=N.r_[rec_w, rec_h, rec_mesh, loc_x, loc_y, loc_z, tilt]
gen_YAML(DNI, sunshape, sunsize, hst_pos, hst_foc, hst_aims,hst_w, hst_h, rho_refl, slope_error, receiver, rec_param, rec_abs, casefolder, spectral=False, medium=0., OneHeliostat=False )
N.savetxt(casefolder+'/azimuth.input', [azimuth])
N.savetxt(casefolder+'/elevation.input', [elevation])
N.savetxt(casefolder+'/rays.input', [num_rays],fmt="%s")
N.savetxt(casefolder+'/mirror.input', [rho_refl])

