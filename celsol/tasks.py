from __future__ import absolute_import, unicode_literals
from .celery import app
import solsticepy
import os
import numpy as np

@app.task
def add(x, y):
    return x + y


@app.task
def mul(x, y):
    return x * y


@app.task
def xsum(numbers):
    return sum(numbers)


@app.task
def solstice_eta_opt(azimuth,elevation,topdir=None):
	"""
	azimuth=270.   # from East to North, deg
	elevation =78. # 0 is horizontal, deg
	"""

	sun = solsticepy.Sun(dni=1000, sunshape='pillbox', half_angle_deg=0.2664)

	# S3. sun position
	# e.g. summer solstice, solar noon
	latitude=34.   # latitude of the crs plant
	# S4. number of rays for the ray-tracing simulation
	num_rays=2000000
	#
	# the field
	# ==========
	# F1.Layout
	layoutfile='/home/john/solstice-scripts/example/demo_layout.csv'

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
	# receiver mesh, for binning the flux distribution
	rec_mesh=100

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

	# extract the heliostat positions from the loaded CSV file
	layout=np.loadtxt(layoutfile, delimiter=',', skiprows=2)
	hst_pos=layout[:,:3]
	hst_foc=layout[:,3] # F2.5
	hst_aims=layout[:,4:] # F4.
	one_heliostat=False

	# define a unique case folder for the user
	snum = 0
	suffix = ""
	while 1:
		import datetime,time
		dt = datetime.datetime.now()
		ds = dt.strftime("%a-%H-%M")
		if topdir is None:
			topdir = os.getcwd()
		casefolder = os.path.join(topdir,'case-%s-%s%s-%d'%(os.path.basename(__file__),ds,suffix,os.getpid()))
		if os.path.exists(casefolder):
		    snum+=1
		    suffix = "-%d"%(snum,)
		    if snum > 200:
		        raise RuntimeError("Some problem with creating casefolder")
		else:
		    # good, we have a new case dir
		    break

	rec_param=np.r_[rec_w, rec_h, rec_mesh, loc_x, loc_y, loc_z, tilt]

	master=solsticepy.Master(casedir=casefolder)
	outfile_yaml = master.in_case('input.yaml')
	outfile_recv = master.in_case('input-rcv.yaml')

	solsticepy.gen_yaml(sun, hst_pos, hst_foc, hst_aims,hst_w, hst_h
		, rho_refl, slope_error, receiver, rec_param, rec_abs
		, outfile_yaml=outfile_yaml, outfile_recv=outfile_recv
		, hemisphere='North', tower_h=tower_h, tower_r=tower_r,  spectral=False
		, medium=0, one_heliostat=one_heliostat)

	# run Solstice using the generate inputs, and run all required post-processing
	eta = master.run(azimuth, elevation, num_rays, rho_refl,sun.dni,gen_vtk=False)

	return (eta.n, eta.s)


