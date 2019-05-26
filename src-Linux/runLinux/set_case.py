from srcPy.gen_YAML import gen_YAML

import numpy as N
import os

solstice_dir='/home/ye/Solstice-0.8.1-GNU-Linux64'
casefolder='./example'
if not os.path.exists(casefolder):
    os.makedirs(casefolder)

#
# the sun
#
# e.g. spring equinox, solar noon
azimuth=270. # from East to North
elevation =53.  # 0 is horizontal

sunshape='pillbox' # or 'buie'
sunsize=0.2664 #deg or CSR value
DNI=1000 # W/m2
num_rays=2000000
#
# the field
#
layout=N.loadtxt('./PS10Layout.csv', delimiter=',', skiprows=2)
hst_pos=layout[:,:3]
hst_foc=layout[:,3]
hst_aims=layout[:,4:]
hst_w=12.925
hst_h=9.575
rho_refl=0.95 # mirror reflectivity
slope_error=2.e-3 # rad
tower_h=115.
tower_r=3.
#
# the receiver
#
receiver='flat'
rec_abs=1.
rec_w=10.
rec_h=10.
rec_mesh=100
loc_x=0.
loc_y=0.
loc_z=115.
tilt=0. # deg
#
#
rec_param=N.r_[rec_w, rec_h, rec_mesh, loc_x, loc_y, loc_z, tilt]
gen_YAML(DNI, sunshape, sunsize, hst_pos, hst_foc, hst_aims,hst_w, hst_h, rho_refl, slope_error, tower_h, tower_r, receiver, rec_param, rec_abs, casefolder, spectral=False, medium=0 )
#
N.savetxt(casefolder+'/azimuth.input', [azimuth])
N.savetxt(casefolder+'/elevation.input', [elevation])
N.savetxt(casefolder+'/rays.input', [num_rays],fmt="%s")
N.savetxt(solstice_dir+'/src-Linux/runLinux/casedir.input', [casefolder], fmt='%s')





