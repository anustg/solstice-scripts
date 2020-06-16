import numpy as np
import platform
import os, sys, subprocess, glob, datetime
import colorama
colorama.init()

from .process_raw import *
from .find_solstice import *
from .cal_sun import *

def yellow(text):
    return colorama.Fore.YELLOW + colorama.Style.BRIGHT + text + colorama.Style.RESET_ALL

def green(text):
    return colorama.Fore.GREEN + colorama.Style.BRIGHT + text + colorama.Style.RESET_ALL

def SPROG(name):
    return find_prog(name)

def run_prog(name,args,output_file=None,verbose=True):
	prog = SPROG(name)
	args1 = [str(a) for a in args]
	if verbose: 
		sys.stderr.write("Running '%s' with args: %s\n" % (name," ".join(args1)))
	if output_file is not None:
		# any error will cause an exception (and we capture the output to a file)
		res = subprocess.check_output([prog]+args1)
		with open(output_file,'w') as f:
			f.write(res.decode('ascii'))
	else:
		# any error will cause an exception...
		subprocess.check_call([prog]+args1)

class Master:

	def __init__(self, casedir='.'):
		"""Set up the Solstice simulation, i.e. establishing the case folder, calling the Solstice program and post-processing the results

		``Argument``
		  * casedir (str): the case directory
		"""
		self.casedir=casedir
		if not os.path.exists(self.casedir):
		    os.makedirs(self.casedir)
		    assert os.path.isdir(casedir)
		sys.stderr.write("Case directory is '%s'\n" % (yellow(self.casedir),))


	def in_case(self, folder, fn):
		"""Joining a file name with the case directory

		``Argument``

		  * fn (str): file name

		``Return``

		  * a joining directory of the file in the case directory
		"""
		
		if not os.path.exists(folder):
		    os.makedirs(folder)
		    assert os.path.isdir(folder)

		return os.path.join(folder,fn)

	def run(self, azimuth, elevation, num_rays, rho_mirror, dni, folder, gen_vtk=True, printresult=True):

		"""Run an optical simulation (one sun position) using Solstice 

		* `azimuth` (float): the azimuth angle of the ray-tracing simulation in Solstice, counted from East towards to North
		* `elevation` (float): the elevation angle of the ray-tracing simulation in Solstice
		* `num_rays` (int): number of rays to be cast in the ray-tracing simulation
		* `rho_mirror` (float): reflectivity of mirrors, required for results post-processing 
		* `dni` (float): the direct normal irradiance (W/m2), required to obtain performance of individual heliostat
		* `gen_vtk` (boolean): if True, generate .vtk files for rendering in Paraview
			
		Returns: no return value (results files are created and written)
		"""

		YAML_IN = self.in_case(self.casedir, 'input.yaml')
		RECV_IN = self.in_case(self.casedir, 'input-rcv.yaml')

		# main raytrace
		run_prog("solstice",['-D%s,%s'%(azimuth,elevation),'-v','-n',num_rays,'-R',RECV_IN,'-fo',self.in_case(folder, 'simul'),YAML_IN])


		if gen_vtk:
			dirn = os.getcwd()
			try:
				os.chdir(folder)

				# Read "simul" results and produce a text file with the raw results
				run_prog('solppraw',[self.in_case(folder, 'simul')])

				# post processing
				run_prog("solstice",['-D%s,%s'%(azimuth,elevation),'-g','format=obj:split=geometry','-fo',self.in_case(folder, 'geom'),YAML_IN])

				# run a short raytrace to produce some ray paths
				run_prog("solstice",['-D%s,%s'%(azimuth,elevation),'-q','-n','100','-R',RECV_IN,'-p','default',YAML_IN], output_file=self.in_case(folder, 'solpaths'))

				# Read "simul" results and produce receiver files (.vtk) of incoming and/or absorbed solar flux per-primitive
				run_prog('solmaps',[self.in_case(folder, 'simul')])

				# Read "geom" and "simul" file results and produce primaries and receivers files (.vtk), and .obj geometry files
				run_prog('solpp',[self.in_case(folder, 'geom'),self.in_case(folder, 'simul')])

				# Read "solpaths" file and produce readable file (.vtk) by paraview to visualize the ray paths
				run_prog('solpaths',[self.in_case(folder, 'solpaths')])
			finally:
				os.chdir(dirn)

		eta, performance_hst=process_raw_results(self.in_case(folder, 'simul'), folder,rho_mirror,dni)

		if printresult:
			sys.stderr.write('\n' + yellow("Total efficiency: {:f}\n".format(eta)))
			sys.stderr.write(green("Completed successfully.\n"))

		return eta, performance_hst

	def run_annual(self, nd, nh, latitude, num_rays, num_hst,rho_mirror,dni,gen_vtk=False):

		"""Run a list of optical simulations to obtain annual performance (lookup table) using Solstice 

		``Arguments``

		  * nd (int): number of rows in the lookup table (discretisation of the declination angle)
		  * nh (int): number of columns in the lookup table (discretisation of the solar hour angle)
		  * latitude (float): the latitude angle of the plan location (deg)
		  * num_rays (int): number of rays to be cast in the ray-tracing simulation
		  * num_hst (int): number of heliostats 
		  * rho_mirror (float): reflectivity of mirrors, required for results post-processing 
		  * dni (float): the direct normal irradiance (W/m2), required to obtain performance of individual heliostat
		  * gen_vtk (bool): True - perform postprocessing for visualisation of  each individual ray-tracing scene (each sun position), False - no postprocessing for visualisation 


		``Return``

		  * No return value (results files are created and written)
		"""

		YAML_IN = self.in_case(self.casedir, 'input.yaml')
		RECV_IN = self.in_case(self.casedir, 'input-rcv.yaml')

		sun=SunPosition()
		AZI, ZENITH,table,case_list=sun.annual_angles(latitude, casefolder=self.casedir, nd=nd, nh=nh)
		case_list=case_list[1:]
		SOLSTICE_AZI, SOLSTICE_ELE=sun.convert_convention('solstice', AZI, ZENITH)

		# performance of individual heliostat is recorded
		# TODO note, DNI is not varied in the simulation, 
		# i.e. performance is not dni-weighted
		ANNUAL=np.zeros((num_hst, 9))    
		run=np.r_[0]

		for i in range(len(case_list)):     
			c=int(case_list[i,0].astype(float))
			if c not in run:
				azimuth=SOLSTICE_AZI[c-1]
				elevation=SOLSTICE_ELE[c-1]
				#if np.sin(elevation*np.pi/180.)>=1.e-5:
				#	dni=1618.*np.exp(-0.606/(np.sin(elevation*np.pi/180.)**0.491))
				#else:
				#	dni=0.

			   
				sys.stderr.write("\n"+green('Sun position: %s \n'%c))
				print('azimuth: %.2f'% azimuth, ', elevation: %.2f'%elevation)

				onesunfolder=os.path.join(self.casedir,'sunpos_%s'%(c))

				# run solstice
				if elevation<1.: # 1 degree
				    efficiency_total=ufloat(0,0)
				    performance_hst=np.zeros((num_hst, 9))  
				else:
					efficiency_total, performance_hst=self.run(azimuth, elevation, num_rays, rho_mirror, dni, folder=onesunfolder, gen_vtk=False, printresult=False)

					sys.stderr.write(yellow("Total efficiency: {:f}\n".format(efficiency_total)))

				ANNUAL+=performance_hst
			else:
				ANNUAL+=performance_hst


			for a in range(len(table[3:])):
				for b in range(len(table[0,3:])):
				    val=re.findall(r'\d+',    table[a+3,b+3])
				    if val==[]:
				        table[a+3,b+3]=0
				    else:
				        if c==float(val[0]):
				            table[a+3,b+3]=efficiency_total.nominal_value

			run=np.append(run,c)   

		annual_title=np.array(['Q_solar','Q_cosine', 'Q_shade', 'Q_hst_abs', 'Q_block', 'Q_atm', 'Q_spil', 'Q_refl', 'Q_rcv_abs']) 
		ANNUAL=np.vstack((annual_title, ANNUAL))
		np.savetxt(self.casedir+'/lookup_table.csv', table, fmt='%s', delimiter=',')
		np.savetxt(self.casedir+'/result-heliostats-annual-performance.csv', ANNUAL, fmt='%s', delimiter=',')
		sys.stderr.write("\n"+green("Lookup table saved.\n"))
		sys.stderr.write(green("Completed successfully.\n"+"\n"))
		return table, ANNUAL


