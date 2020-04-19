#
# This is an example script to simulate a central receiver system (CRS)
# using Solstice via solsticepy
#
import solsticepy
from solsticepy.master import Master
import numpy as np
import os

#==================================================
# INPUT PARAMETERS

# whether run a new case (True) or load pre-existing input.yaml inputs (False):
new_case=True 

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
azimuth=270.   # from East to North, deg
elevation =78. # 0 is horizontal, deg
latitude=34.   # latitude of the crs plant
# S4. number of rays for the ray-tracing simulation
num_rays=2000000
#
# the field
# ==========
# F1.Layout
layoutfile='./demo_layout.csv'
hst_field=True # simulate the full heliostat field (if False: give single heliostat information)

# F2. Heliostat
hst_w=10. # m
hst_h=10. # m
rho_refl=0.95 # mirror reflectivity
slope_error=2.e-3 # radians
# F3. Tower
tower_h=0.01 # tower height
tower_r=0.01 # tower radius
#
# the receiver
# ============
# R1. shape
receiver='flat' # 'flat' or 'stl'
# R2. Size
rec_w=8. # width, m
rec_h=6. # height, m
# R3. tilt angle
tilt=0.  # deg
# R4. position
loc_x=0. # m
loc_y=0. # m
loc_z=62.# m
# R5. Abosrptivity
rec_abs=0.9

if receiver=='flat':
    # receiver mesh, for binning the flux distribution
    rec_mesh=100
elif receiver=='stl':
    stlfile='./demo_plane.stl'

# set the folder name for saving the output files
# False: an automatically generated name  
# or
# string: name of the user defined folder
# (Note: you must set this folder name if you want to use new_case==False)
userdefinedfolder=False


# NO NEED TO CHANGE THE CONTENT BELOW
# ===============================================================
# the ray-tracing scene will be generated 
# based on the settings above

if hst_field:
	# extract the heliostat positions from the loaded CSV file
    layout=np.loadtxt(layoutfile, delimiter=',', skiprows=2)
    hst_pos=layout[:,:3]
    hst_foc=layout[:,3] # F2.5
    hst_aims=layout[:,4:] # F4.
    one_heliostat=False
else:
    one_heliostat=True

if new_case==False:
	assert os.path.isdir(userdefinedfolder)

if userdefinedfolder:
    casefolder=userdefinedfolder
else:
	# define a unique case folder for the user
    snum = 0
    suffix = ""
    while 1:
        import datetime,time
        dt = datetime.datetime.now()
        ds = dt.strftime("%a-%H-%M")
        casefolder = os.path.join(os.getcwd(),'case-%s-%s%s'%(os.path.basename(__file__),ds,suffix))
        if os.path.exists(casefolder):
            snum+=1
            suffix = "-%d"%(snum,)
            if snum > 200:
                raise RuntimeError("Some problem with creating casefolder")
        else:
            # good, we have a new case dir
            break

if receiver =='flat':
    rec_param=np.r_[rec_w, rec_h, rec_mesh, loc_x, loc_y, loc_z, tilt]
elif receiver =='stl':
    rec_param=np.array([rec_w, rec_h, stlfile, loc_x, loc_y, loc_z, tilt])

master=Master(casedir=casefolder)
outfile_yaml = master.in_case('input.yaml')
outfile_recv = master.in_case('input-rcv.yaml')

if new_case:
	# generate the YAML file from the input parameters specified above
    solsticepy.gen_yaml(DNI, sunshape, sunsize, hst_pos, hst_foc, hst_aims,hst_w, hst_h
		, rho_refl, slope_error, receiver, rec_param, rec_abs
		, outfile_yaml=outfile_yaml, outfile_recv=outfile_recv
		, hemisphere='North', tower_h=tower_h, tower_r=tower_r,  spectral=False
		, medium=0, one_heliostat=one_heliostat)

# run Solstice using the generate inputs, and run all required post-processing
master.run(azimuth, elevation, num_rays, rho_refl,DNI)

# annual solution (see instructions)
#master.run_annual(nd=5, nh=5, latitude=latitude, num_rays=num_rays, num_hst=len(hst_pos),rho_mirror=rho_refl, dni=DNI)

