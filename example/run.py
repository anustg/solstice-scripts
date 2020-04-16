
# TODO
# Suggested principle: the user should copy and modify this script. The script
# Can start its life in the c:\Program Files\solstice-0.9.0\ folder, and then
# be copied by the user to their home directory. All output files will be
# placed in a subdirectory in where the script is located, with a name that
# includes the script name, plus the date/time. Students need only double-click
# the script file to run it in Windows. We may need to add a way of
# alerting the user before the script window closes.

import os, sys, subprocess, glob, datetime
import numpy as np

#-------------------------------------------------------------------------------
# INPUT PARAMETERS / SETUP

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
#
# the field
# ==========

# F1.Layout
# file is assumed to be located in the same directory as this script!
layout=np.loadtxt('demo_layout.csv', delimiter=',', skiprows=2)
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
receiver='stl' # 'flat', 'cylinder' or 'stl'

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

if receiver=='flat':
    # receiver mesh, for binning the flux distribution
    rec_mesh=100

elif receiver=='stl':
    stlfile='./demo_plane.stl'

#-------------------------------------------------------------------------------
# GENERATE THE INPUT FILES, RUN THE ANALYSIS

#TODO: we must decide whether it is better to wrap up any of this part into
# the SolsticePy scripts which are imported. Depends whether the students
# will want to customise any of the following, or not.

import colorama
colorama.init()
import solsticepy

snum = 0
suffix = ""

# check that the script itself is user-writable
if not os.access(__file__,os.W_OK):
	raise RuntimeError("Your script '%s' is not writeable. Create a copy of \
the file somewhere in your home directory in order to run simulations. \
The script will create data-output directories wherever your script is located."%(os.path.basename(__file__),))

# create a subfolder based on the name of the script (plus a timestamp)
while 1:
	import datetime,time
	dt = datetime.datetime.now()
	ds = dt.strftime("%a-%H:%M")
	case_dir = os.path.join(os.getcwd(),'case-%s-%s%s'%(os.path.basename(__file__),ds,suffix))
	if os.path.exists(case_dir):
		snum+=1
		suffix = "-%d"%(snum,)
		if snum > 200:
			raise RuntimeError("Some problem with creating case_dir")
	else:
		# good, we have a new case dir
		break

def SPROG(n):
	return solsticepy.find_prog(n)
SOLSTICE = SPROG('solstice')

def in_case(fn):
	return os.path.join(case_dir,fn)

def yellow(text):
	return colorama.Fore.YELLOW + colorama.Style.BRIGHT + text + colorama.Style.RESET_ALL
def green(text):
	return colorama.Fore.GREEN + colorama.Style.BRIGHT + text + colorama.Style.RESET_ALL


def run_prog(prog,args,output_file=None):
	args1 = [str(a) for a in args]
	sys.stderr.write("Running '%s' with args: %s\n" % (prog," ".join(args1)))
	if output_file is not None:
		# any error will cause an exception (and we capture the output to a file)
		res = subprocess.check_output([prog]+args1)
		with open(output_file,'w') as f:
			f.write(res.decode('ascii'))
	else:
		# any error will cause an exception...
		subprocess.check_call([prog]+args1)


# FIXME I haven't tested the yaml=True case, but it would cause problems
# because the YAML_IN file is assumed to be located in the case_dir.

yaml=False # whether run an existing yaml file directly?
           # NOTE, the yaml files must be in the casefolder

YAML_IN = in_case('input.yaml')
RECV_IN = in_case('input-rcv.yaml')

# create the ray-tracing scene based on the settings above

hst_pos=layout[:,:3]
hst_foc=layout[:,3] # F2.5
hst_aims=layout[:,4:] # F4.

if not os.path.exists(case_dir):
	os.makedirs(case_dir)
	assert os.path.isdir(case_dir)

sys.stderr.write("Case directory is '%s'\n" % (yellow(case_dir),))

if receiver =='flat':
    rec_param=np.r_[rec_w, rec_h, rec_mesh, loc_x, loc_y, loc_z, tilt]
elif receiver =='stl':
    rec_param=np.array([rec_w, rec_h, stlfile, loc_x, loc_y, loc_z, tilt])


if not yaml:
	sys.stderr.write("Generating YAML input files...\n");
	solsticepy.gen_yaml(DNI, sunshape, sunsize, hst_pos, hst_foc, hst_aims,hst_w, hst_h
		, rho_refl, slope_error, receiver, rec_param, rec_abs
		, tower_h=tower_h, tower_r=tower_r, spectral=False, medium=0.
		, one_heliostat=False
		, outfile_yaml = YAML_IN
		, outfile_recv = RECV_IN
	)

# run solstice
#

# main raytrace
run_prog(SOLSTICE,['-D%s,%s'%(azimuth,elevation),'-v','-n',num_rays,'-R',RECV_IN,'-fo',in_case('simul'),YAML_IN])

# post processing

run_prog(SOLSTICE,['-D%s,%s'%(azimuth,elevation),'-g','format=obj:split=geometry','-fo',in_case('geom'),YAML_IN])
run_prog(SOLSTICE,['-D%s,%s'%(azimuth,elevation),'-q','-n','100','-R',RECV_IN,'-p','default',YAML_IN], output_file=in_case('solpaths'))

# post process (inside our case_dir):
popd = os.getcwd()
try:
	os.chdir(case_dir)

	# Read "simul" results and produce a text file with the raw results
	run_prog(SPROG('solppraw'),[in_case('simul')])

	# Read "simul" results and produce receiver files (.vtk) of incoming and/or absorbed solar flux per-primitive
	run_prog(SPROG('solmaps'),[in_case('simul')])

	# Read "geom" and "simul" file results and produce primaries and receivers files (.vtk), and .obj geometry files
	run_prog(SPROG('solpp'),[in_case('geom'),in_case('simul')])

	# Read "solpaths" file and produce readable file (.vtk) by paraview to visualize the ray paths
	run_prog(SPROG('solpaths'),[in_case('solpaths')])

	# create summary data files 'result-raw.csv' and 'result-formatted.csv'
	eta = solsticepy.process_raw_results(in_case('simul'), case_dir,rho_refl)
	sys.stderr.write('\n' + yellow("Total efficiency: %s\n"%(repr(eta),)))

finally:
	os.chdir(popd)

sys.stderr.write(green("Completed successfully.\n"))

