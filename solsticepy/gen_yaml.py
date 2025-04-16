from __future__ import print_function
import numpy as np
import matplotlib.pyplot as plt

#for python 2:
#from builtins import super

from .data_spectral import SolarSpectrum, MirrorRhoSpectrum
from .cal_layout import multi_aperture_pos
import sys

def yamltransform(pos,rot):
    return "transform: { translation: [%e,%e,%e], rotation: [%e,%e,%e] }" % (*pos,*rot)

class Sun:
	"""Sun parameters for solstice-input

	Example:
	>>> sun = Sun(dni=1000, sunshape='buie', csr=0.2)
	>>> print(sun.yaml(spectrum = "*solar_spectrum"))

	"""
	def __init__(self,dni=1000,sunshape=None,csr=0.01,half_angle_deg=0.2664,std_dev=0.2):
		"""Define sun parameters for Solstice input file.

		`dni`: Direct normal irradance (W/m2)
		`sunshape`: Sunshape: can be None, ``'pillbox'``,``'gaussian'`` or ``'buie'``
		`csr`: circumsolar ratio (ONLY in case of ``'buie'``)
		`half_angle_deg`: sun angular size (in DEGREES, half-angle) (ONLY in case of ``'pillbox'``)
		`std_dev`: standard deviation of the angular dsn ratio (ONLY in case of ``'gaussian'``)
		"""
		self.dni = dni
		self.sunshape = sunshape
		if sunshape is not None:
			assert sunshape in ['buie','pillbox','gaussian']
			self.sunshape = sunshape
			if sunshape == "pillbox":
				self.half_angle_deg = half_angle_deg
			elif sunshape == "buie":
				self.csr = csr
			elif sunshape == "gaussian":
				self.std_dev = std_dev
	def yaml(self,spectrum = None):
		"""YAML representation of the sun for solstice-input.

		`spectrum`: YAML value of the solar spectrum. Would normally be set to the string ``"*solar_spectrum"``.
		"""
		# FIXME surely we find a smarter way to do this serialization AND deseralization with YAML?
		s = "{dni: %15.8e" % (self.dni,)
		if spectrum is not None:
			s += ", spectrum = %s" % (spectrum,)
		if self.sunshape is not None:
			if self.sunshape=='pillbox':
				s += ", pillbox: {half_angle: %6.4f}" % (self.half_angle_deg,)   
			elif self.sunshape=='buie':
				s += ", buie: {csr: %6.4f}" % (self.csr,) 
			elif self.sunshape=='gaussian':
				s += ", gaussian: {std_dev: %6.4f}" % (self.std_dev,)
		s += "}"
		return s


def gen_yaml(sun, hst_pos, hst_foc, hst_aims,hst_w, hst_h
		, rho_refl, slope_error, receiver, rec_param, rec_abs
		, outfile_yaml, outfile_recv
		, hemisphere='North', tower_h=0.01, tower_r=0.01,  spectral=False
		, medium=0, one_heliostat=False
):
	"""Generate the heliostat field and receiver YAML input files for Solstice ray-tracing simulation.

	1. the sun
	  * `sun` (`Sun` object): parameters relating to the solar source

	2. the field
	  * `hst_pos` (nx3 numpy array): heliostat positions (x, y, z) (first of the 'field' parameters)
	  * `hst_foc` (nx1 numpy array): heliostat focal length
	  * `hst_aims` (nx3 numpy array): heliostat aiming point (ax, ay, az)
	  * `hst_w` (float): heliostat mirror width  (in x direction)
	  * `hst_h` (float): heliostat mirror height (in y direction)
	  * `hst_z` (float): heliostat center height (in z direction)
	  * `rho_refl` (float): reflector reflectivity
	  * `slope_error` (float): reflector surface slope error rms, radians
	  * `tower_h` (float): tower height (m)
	  * `tower_r` (float): tower radius (a cylindrical shape tower) (m)
	3. the receiver
	  * `receiver` (str): ``'flat'``, ``'cylinder'``, or ``'stl' or 'multi-aperture'`` (first of the 'receiver' parameters)
	  * `rec_abs` (float): receiver absorptivity
	  * `rec_param` (numpy array or str): each element contains the geometrical parameter of the corresponding receiver.
	4. others
	  * `spectral` (bool): True - simulate the spectral dependent performance (first of the 'other' parameters)
	  * `medium` (float): if the atmosphere is surrounded by non-participant medium, medium=0; otherwise it is the extinction coefficient in m-1
	  * `one_heliosat` (boolean): if `True`, implements ray tracing from just one heliostat.
	  	
	Returns: nothing (requested files are created and written)

	Note that the parameters are in groups that relate to the `sun`, the `field` and the `receiver` then `others`. 

	Note also the type for `rec_param` should be as follows.
	  * if ``receiver == 'flat'``: np.array([width, height, grid_w, grid_h,, x, y, z, tilt angle (deg))]
	  * if ``receiver == 'cylinder'``: np.array([radius, height, grid_circ, grid_h, x, y, z, tilt angle (deg)])
	  * if ``receiver == 'stl'``: the directory of the stl file
	  * if ``receiver == 'multi-aperture'``:  np.array([width, height, grid_w, grid_h,, x, y, z, tilt angle (deg),num_aperture, gamma (deg) ])
	"""
	# FIXME Parameters should be named according to what they are, eg
	# the parameter should be called 'csr', not 'sunsize', to avoid confusion.
	# We can still improve our calling convention, to make this library easier
	# to use and more maintainable.

	sys.stderr.write("Generating YAML file...\n")

	iyaml='' # the input yaml file

	# 
	### Section (1)
	# set the spectral data: 
	# solar radiative intensity, refractive indexes, extinction coefficients, reflectivities
	#------------------------------
	if spectral:
		I_sun=SolarSpectrum()
		# CREATE the spectrum for the sun
		iyaml+='- spectrum: &solar_spectrum  \n'
		for i in range(0,len(I_sun)-1):
		    iyaml+='  - {wavelength: %e, data: %e }\n' % (I_sun[i][0],I_sun[i][1])  
		i = len(I_sun)-1
		iyaml+='  - {wavelength: %e, data: %e }\n' % (I_sun[i][0],I_sun[i][1]) 
		iyaml+='\n'

		# CREATE the spectrum for the reflectivity (mirror)
		mirror_rho= MirrorRhoSpectrum()
		mirror_ref=mirror_rho
		for i in range(0,len(mirror_rho)):
		    mirror_ref[i][0] = mirror_rho[len(mirror_rho)-1-i][0]/1000.
		    mirror_ref[i][1] = mirror_rho[len(mirror_rho)-1-i][1]/100.
		mirror_ref.append([4,0.9])
		iyaml+='- spectrum: &%s  \n' % 'ref_mirror'
		for i in range(0,len(mirror_ref)-1):
		    iyaml+='  - {wavelength: %15.8e, data: %15.8e }\n' % (float(mirror_ref[i][0]),float(mirror_ref[i][1])) 
		i = len(mirror_ref)-1
		iyaml+='  - {wavelength: %15.8e, data: %15.8e }\n' % (float(mirror_ref[i][0]),float(mirror_ref[i][1])) 
		iyaml+='\n'

	# 
	### Section (2)
	# set the medium types: 
	# air, glass, vacuum, etc. gathering spectral data
	#------------------------------
	#

	#
	# Creation of the sun and atmosphere
	#
	if spectral:
		spectrum = "*solar_spectrum"
	else:
		spectrum = None
	
	iyaml += "- sun: %s\n" % (sun.yaml(spectrum),)

	if medium>1e-99:
		iyaml+='- atmosphere: {extinction: %e}\n'%medium 
		iyaml+='\n'

		   
	# 
	### Section (3)
	# set the materials
	# (gathering media)
	# occultant material, mirror specular material, receiver material, virtual target
	#------------------------------
	#
	# CREATE an occultant material
	iyaml+='- material: &material_black\n'
	iyaml+='   front:\n'
	iyaml+='     matte: {reflectivity: 0.}\n' # front    
	iyaml+='   back:\n'
	iyaml+='     matte: {reflectivity: 0. }\n' # and back reflectivity
	iyaml+='\n'
	#
	# CREATE a specular material
	iyaml+='- material: &material_mirror\n'
	iyaml+='   front:\n'
	if spectral:
		iyaml+='     mirror: {reflectivity: *ref_mirror, slope_error: %15.8e }\n' % (slope_error ) 
	else:
		iyaml+='     mirror: {reflectivity: %6.4f, slope_error: %15.8e }\n' % (rho_refl, slope_error) 
	iyaml+='   back:\n'
	iyaml+='     matte: {reflectivity: 0. }\n'
	iyaml+='\n'
	#
	# CREATE a material for the target
	r_f = 1.-rec_abs # front
	r_b = 1.-rec_abs # and back reflectivity
	iyaml+='- material: &material_target\n'
	iyaml+='   front:\n'
	iyaml+='     matte: {reflectivity: %6.4f }\n' % r_f    
	iyaml+='   back:\n'
	iyaml+='     matte: {reflectivity: %6.4f }\n' % r_b
	iyaml+='\n'
	#
	# CREATE a virtual material for the calculation of spillage
	iyaml+='- material: &material_virtual\n'
	iyaml+='   virtual:\n'
	iyaml+='\n'


	# 
	### Section (4)
	# set the geometries
	# (gathering shapes and materials)
	# the tower, the receiver, the heliostat
	#------------------------------
	#
	# Tower Geometry
	# (cylindrical shape)
	#
	slices = 10 # slices for the envelop circle
	iyaml+='- geometry: &tower_g\n' 
	iyaml+='  - material: *material_black\n' 
	#iyaml+='    transform: { translation: %s, rotation: %s }\n' % ([0, 0, h_tow*0.5], [0, 90, 0]) 
	iyaml+='    cylinder: {height: %7.3f, radius: %7.3f, slices: %d }\n' % (tower_h, tower_r, slices) 
	iyaml+='\n'
	#
	# Receiver Geometry
	#
	if receiver=='flat':
		geom, rec_entt, rcv = flat_receiver(rec_param, hemisphere)
		iyaml+=geom

	elif receiver=='cylinder':
		geom, rec_entt, rcv = cylindrical_receiver(rec_param, hemisphere)
		iyaml+=geom

	elif receiver=='stl':
		rec_entt, rcv=STL_receiver(rec_param, hemisphere)

	elif receiver=='multi-aperture':
		geom, rec_entt, rcv =multi_aperture_receiver(rec_param, hemisphere)
		iyaml+=geom
	#
	# Heliostats Geometry
	#
	if one_heliostat:
		hst_x=np.r_[hst_pos[0]]
		hst_y=np.r_[hst_pos[1]]
		hst_z=np.r_[hst_pos[2]]
		aim_x=np.r_[hst_aims[0]] 
		aim_y=np.r_[hst_aims[1]]
		aim_z=np.r_[hst_aims[2]]
		num_hst=1
		hst_foc=np.r_[hst_foc]
	else:
		hst_x=hst_pos[:,0]
		hst_y=hst_pos[:,1]
		hst_z=hst_pos[:,2]
		aim_x=hst_aims[:,0]
		aim_y=hst_aims[:,1]
		aim_z=hst_aims[:,2]
		num_hst=len(hst_x)
	slices = 4 # slices for the envelop circle
	
	# CREATE a reflective facet (mirror)
	for i in range(0,num_hst):
		name_hst_g = 'hst_g_'+str(i)
		iyaml+='- geometry: &%s\n' % name_hst_g 
		iyaml+='  - material: *material_mirror\n' 
		#iyaml+='    transform: { translation: %s, rotation: %s }\n' % ([hst_x[i], hst_y[i], hst_z[i]], [0, 0, 0]) )
		iyaml+='    parabol: \n'
		iyaml+='      focal: %s\n' % hst_foc[i]
		iyaml+='      clip: \n'  
		iyaml+='      - operation: AND \n'
		iyaml+='        vertices: [ [%e, %e], [%e, %e], [%e, %e], [%e, %e] ]\n' % (-hst_w*0.5, -hst_h*0.5, -hst_w*0.5, hst_h*0.5, hst_w*0.5, hst_h*0.5, hst_w*0.5,-hst_h*0.5)
		iyaml+='      slices: %d\n' % slices  

	# CREATE the pylon "pylon_g" geometry cylindrical shape
	h_pyl = 0.001 # pylon height
	r_pyl = 0.2 # pylon radius
	slices = 4 # slices for the envelop circle
	iyaml+='- geometry: &pylon_g\n'
	iyaml+='  - material: *material_black\n' 
	iyaml+='    '+yamltransform(pos=[0,0,-h_pyl*3],rot=[0,90,0]) + '\n'
	iyaml+='    cylinder: {height: %7.3f, radius: %7.3f, slices: %d }\n' % (h_pyl,r_pyl,slices) 
	#   

	# 
	### Section (5)
	# set the templates
	# (programming objects gathering geometries or pivot and geometries)
	#------------------------------
	# CREATE the heliostat templates
	for i in range(0,num_hst):    
		name_hst_t = 'hst_t_'+str(i)
		iyaml+='- template: &%s\n' % name_hst_t 
		name_hst_n = 'hst_'+ str(i)
		iyaml+='    name: %s\n' % name_hst_n 
		iyaml+='    primary: 0\n'   
		iyaml+='    geometry: *pylon_g\n'
		iyaml+='    children: \n' 
		iyaml+='    - name: pivot\n'
		iyaml += '      zx_pivot: {target: {position: [%.6f, %.6f, %.6f]}}\n' % (aim_x[i], aim_y[i], aim_z[i])
		iyaml+='      children: \n'
		iyaml+='      - name: reflect_surface\n'
		iyaml+='        primary: 1\n'
		iyaml+='        transform: {rotation: [-90,0,0]} \n'   
		name_hst_g = 'hst_g_'+str(i)
		iyaml+='        geometry: *%s\n' % name_hst_g 

	# 
	### Section (6)
	# set the entities
	# (gather templates to be created and active in the scene)
	#------------------------------
	#
	# receiver entities
	iyaml+=rec_entt
	#
	# tower entities
	iyaml+='\n- entity:\n'
	iyaml+='    name: tower_e\n'
	iyaml+='    primary: 0\n' 
	iyaml+='    ' + yamltransform(pos=[0,-tower_r, tower_h*0.5],rot=[0,0,0]) + '\n'
	iyaml+='    geometry: *%s\n' % 'tower_g'    
	#
	# heliostat entities from the template
	for i in range(0,num_hst):
		name_e ='H_'+str(i)
		name_hst_t = 'hst_t_'+str(i)
		iyaml+='\n- entity:\n'
		iyaml+='    name: %s\n' % name_e
		iyaml+='    ' + yamltransform(pos=[hst_x[i], hst_y[i], hst_z[i]],rot=[0,0,0]) + '\n'
		iyaml+='    children: [ *%s ]\n' % name_hst_t    

	with open(outfile_yaml,'w') as f:
		f.write(iyaml)

	with open(outfile_recv,'w') as f:
		f.write(rcv) 


def flat_receiver(rec_param, hemisphere='North'):
	"""
	hemisphere : 'North' or 'South' hemisphere of the earth where the field located
		        if North: the field is in the positive y direction
		        if South: the field is in the negtive y direction
		        this will influence:
		         (1) the setting of the receiver tilt angle, 
		             if the front surface always facing to the field is desirable
		         (2) the position of the virtual target
	"""
	rec_w=float(rec_param[0])
	rec_h=float(rec_param[1])
	slices=int(rec_param[2]) # it assumes equal number of slices in x and y directions
	x=float(rec_param[4])
	y=float(rec_param[5])
	z=float(rec_param[6])
	tilt=float(rec_param[7])
	# receiver tilt angle:
	# 0 is vertical
	# the standby posiion of a plane in solstice is normal points to the +z axis
	# rotation anagle, positive is anti-clockwise

	geom=''

	geom+='- geometry: &target_g\n'
	geom+='  - material: *material_target\n'
	geom+='    plane: \n'
	geom+='      clip: \n' 
	geom+='      - operation: AND \n'
	geom+='        vertices: [ [%e, %e], [%e, %e], [%e, %e], [%e, %e] ] \n' % (-rec_w*0.5, -rec_h*0.5, -rec_w*0.5, rec_h*0.5, rec_w*0.5, rec_h*0.5, rec_w*0.5,-rec_h*0.5)
	geom+='      slices: %d\n' % slices 
	geom+='\n'

	# CREATE a receiver entity from "target_g" geometry (primary = 0)
	entt=''
	entt+='\n- entity:\n'
	entt+='    name: target_e\n'
	entt+='    primary: 0\n'
	if hemisphere=='North':
		entt+='    ' + yamltransform(pos=[x, y , z],rot=[-90.-tilt, 0, 0]) + '\n'
	else:
		entt+='    ' + yamltransform(pos=[x, y , z],rot=[90.+tilt, 0, 0]) + '\n'
	entt+='    geometry: *%s\n' % 'target_g'

	# CREATE a virtual target entity from "target_g" geometry (primary = 0)
	 
	slices = 4
	entt+='\n- entity:\n'
	entt+='    name: virtual_target_e\n'
	entt+='    primary: 0\n'
	if hemisphere=='North':
		entt+='    ' + yamltransform(pos=[x, y-5, z],rot=[-90.-tilt, 0, 0]) + '\n'
	else:
		entt+='    ' + yamltransform(pos=[x, y+5 , z],rot=[90.+tilt, 0, 0]) + '\n'
	entt+='    geometry: \n' 
	entt+='      - material: *%s\n' % 'material_virtual' 
	entt+='        plane: \n'
	entt+='          clip: \n'    
	entt+='          - operation: AND \n'
	entt+='            vertices: [ [%e, %e], [%e, %e], [%e, %e], [%e, %e] ] \n' % (-rec_w*10., -rec_h*10., -rec_w*10., rec_h*10., rec_w*10., rec_h*10., rec_w*10.,-rec_h*10.)
	entt+='          slices: %d\n' % slices  

	rcv=''
	rcv+='- name: target_e \n' 
	rcv+='  side: %s \n' % 'FRONT_AND_BACK'
	rcv+='  per_primitive: %s \n' % 'INCOMING_AND_ABSORBED'
	rcv+='- name: virtual_target_e \n'
	rcv+='  side: %s \n' % 'FRONT'
	rcv+='  per_primitive: %s \n' % 'INCOMING'

	return geom, entt, rcv



def cylindrical_receiver(rec_param, hemisphere='North'):
	'''
	hemishpere : 'North' or 'South' hemisphere of the earth where the field located
		        if North: the field is in the positive y direction
		        if South: the field is in the negtive y direction
		        this will influence:
		         (1) the setting of the receiver tilt angle, 
		             if the front surface always facing to the field is desirable
		         (2) the position of the virtual target
	'''
	rec_r=float(rec_param[0]/2.)
	rec_h=float(rec_param[1])
	slices=int(rec_param[2]) # number of elements in the circumferetial direction
	stacks=int(rec_param[3]) # number of elements in the vertical direction
	x=float(rec_param[4])
	y=float(rec_param[5])
	z=float(rec_param[6])

	geom=''
	geom+='- geometry: &%s\n' % 'target_g'
	geom+='  - material: *%s\n' % 'material_target'
	geom+='    cylinder: \n'
	geom+='      height: %s\n'%rec_h 
	geom+='      radius: %s\n'%rec_r 
	geom+='      slices: %d\n' % slices 
	geom+='      stacks: %d\n' % stacks 
	geom+='\n'

	# CREATE a receiver entity from "target_g" geometry (primary = 0)
	entt=''
	entt+='\n- entity:\n'
	entt+='    name: target_e\n'
	entt+='    primary: 0\n'
	entt+='    ' + yamltransform(pos=[x, y , z],rot=[0, 0, 0]) + '\n'
	entt+='    geometry: *%s\n' % 'target_g'

	# CREATE a virtual target entity from "target_g" geometry (primary = 0)
	Vsize=100.

	slices = 4
	entt+='\n- entity:\n'
	entt+='    name: virtual_target_e\n'
	entt+='    primary: 0\n'
	entt+='    ' + yamltransform(pos=[x, y, z+rec_h/2.+1],rot=[-180., 0, 0]) + '\n'
	entt+='    geometry: \n' 
	entt+='      - material: *%s\n' % 'material_virtual' 
	entt+='        plane: \n'
	entt+='          clip: \n'    
	entt+='          - operation: AND \n'
	entt+='            vertices: [ [%e, %e], [%e, %e], [%e, %e], [%e, %e] ]\n' % (-rec_h*Vsize, -rec_h*Vsize, -rec_h*Vsize, rec_h*Vsize, rec_h*Vsize, rec_h*Vsize, rec_h*Vsize,-rec_h*Vsize)
	entt+='          slices: %d\n' % slices  

	rcv=''
	rcv+='- name: target_e \n' 
	rcv+='  side: %s \n' % 'FRONT_AND_BACK'
	rcv+='  per_primitive: %s \n' % 'INCOMING_AND_ABSORBED'
	rcv+='- name: virtual_target_e \n'
	rcv+='  side: %s \n' % 'FRONT'
	rcv+='  per_primitive: %s \n' % 'INCOMING'

	return geom, entt, rcv

    

def STL_receiver(rec_param, hemisphere='North'):
	'''
	hemishpere : 'North' or 'South' hemisphere of the earth where the field located
		        if North: the field is in the positive y direction
		        if South: the field is in the negtive y direction
		        this will influence:
		         (1) the setting of the receiver tilt angle, 
		             if the front surface always facing to the field is desirable
		         (2) the position of the virtual target
	'''

	rec_w=float(rec_param[0]) # for creating the virtual target
	rec_h=float(rec_param[1])
	stlfile=rec_param[2] # directory of the stl file
	x=float(rec_param[3])
	y=float(rec_param[4])
	z=float(rec_param[5])
	tilt=float(rec_param[6]) # need to figure out the initial mesh orientation

	# CREATE a receiver entity from a STL file 
	entt=''
	entt+='\n- entity:\n'
	entt+='    name: STL_receiver_e\n'
	entt+='    primary: 0\n'
	if hemisphere=='North':
		entt+='    ' + yamltransform(pos=[x, y, z],rot=[-90.-tilt, 0, 0]) + '\n'
	else:
		# if it is the mesh model of the bladed receiver at CSIRO
		entt+='    ' + yamltransform(pos=[x, y, z],rot=[180.+tilt, 0, 0]) + '\n'
	entt+='    geometry:\n'
	entt+='    - material: *material_target\n'
	entt+='    ' + yamltransform(pos=[0,0,0],rot=[0, 0, 0]) + '\n'
	entt+="      stl : {path: %s }  \n"%(stlfile)


	# CREATE a virtual target entity from "target_g" geometry (primary = 0)

	slices = 4
	entt+='\n- entity:\n'
	entt+='    name: virtual_target_e\n'
	entt+='    primary: 0\n'
	if hemisphere=='North':
		entt+='    ' + yamltransform(pos=[x, y-5., z],rot=[-90.-tilt, 0, 0]) + '\n'
	else:
		entt+='    ' + yamltransform(pos=[x, y+5., z],rot=[90.+tilt, 0, 0]) + '\n'
	entt+='    geometry: \n' 
	entt+='      - material: *%s\n' % 'material_virtual' 
	entt+='        plane: \n'
	entt+='          clip: \n'    
	entt+='          - operation: AND \n'
	entt+='            vertices: [ [%e, %e], [%e, %e], [%e, %e], [%e, %e] ] \n' % (-rec_w*10., -rec_h*10., -rec_w*10., rec_h*10., rec_w*10., rec_h*10., rec_w*10.,-rec_h*10.)
	entt+='          slices: %d\n' % slices  

	rcv=''
	rcv+='- name: STL_receiver_e \n' 
	rcv+='  side: %s \n' % 'FRONT_AND_BACK'
	rcv+='  per_primitive: %s \n' % 'INCOMING_AND_ABSORBED'
	rcv+='- name: virtual_target_e \n'
	rcv+='  side: %s \n' % 'FRONT'
	rcv+='  per_primitive: %s \n' % 'INCOMING'

	return entt, rcv

def multi_aperture_receiver(rec_param, hemisphere='North', plot=False):
	"""
	hemisphere : 'North' or 'South' hemisphere of the earth where the field located
		        if North: the field is in the positive y direction
		        if South: the field is in the negtive y direction
		        this will influence:
		         (1) the setting of the receiver tilt angle, 
		             if the front surface always facing to the field is desirable
		         (2) the position of the virtual target
	"""
	# rec_w and rec_h is the size of one aperture
	# rec_grid_w and rec_gird_h is the number of elements of one aperture
	# rec_z is a list of the elevation height of the center of apertures
	# rec_tilt is the tilt angle of each aperture (the default is facing to the horizon)
	# num_aperture is the number of apertures
	# gamma is the angular range of the multi-aperture configration 


	rec_w=rec_param[0]
	rec_h=rec_param[1]
	rec_grid_w=int(rec_param[2])
	rec_grid_h=int(rec_param[3])

	rec_z=rec_param[4]
	rec_tilt=rec_param[5]
	# receiver tilt angle:
	# 0 is vertical
	# the standby posiion of a plane in solstice is normal points to the +z axis
	# rotation anagle, positive is anti-clockwise
	num_aperture=int(rec_param[6]) 
	gamma=rec_param[7]  # angular range of the multi-aperture configration (deg)


	geom=''
	entt=''
	vir_z=0.
	for i in range(num_aperture):

		geom+='- geometry: &%s\n' % 'target_g_%.0f\n'%(i)
		geom+='  - material: *%s\n' % 'material_target'
		geom+='    plane: \n'
		geom+='      clip: \n' 
		geom+='      - operation: AND \n'
		geom+='        vertices: [ [%e, %e], [%e, %e], [%e, %e], [%e, %e] ]\n' % (-rec_w[i]*0.5, -rec_h[i]*0.5, -rec_w[i]*0.5, rec_h[i]*0.5, rec_w[i]*0.5, rec_h[i]*0.5, rec_w[i]*0.5,-rec_h[i]*0.5)
		geom+='      slices: %d\n' % rec_grid_w 
		geom+='\n'

		ang_pos, xc, yc=multi_aperture_pos(rec_w, gamma, num_aperture, i)

		zc=rec_z[i]		
		vir_z+=zc

		# CREATE a receiver entity from "target_g" geometry (primary = 0)

		entt+='\n- entity:\n'
		entt+='    name: target_e_%.0f\n'%(i)
		entt+='    primary: 0\n'
		if hemisphere=='North':
			entt+='    ' + yamltransform(pos=[xc, yc, zc],rot=[-90.-rec_tilt, 90.-ang_pos,0]) + '\n'
		else:
			entt+='    ' + yamltransform(pos=[-xc, -yc, zc],rot=[90.+rec_tilt, 90.-ang_pos,0]) + '\n'
		entt+='    geometry: *%s\n' % 'target_g_%.0f\n'%(i)

	vir_z/=float(num_aperture)

	# CREATE a virtual target entity from "target_g" geometry (primary = 0)
	slices = 16
	radius=vir_z*0.5
	entt+='\n- entity:\n'
	entt+='    name: virtual_target_e\n'
	entt+='    primary: 0\n'
	entt+='    transform: { translation: %s}\n' % ([0., 0., vir_z])
	entt+='    geometry: \n' 
	entt+='      - material: *%s\n' % 'material_virtual' 
	entt+='        sphere: \n'
	entt+='          radius: %s\n' % radius   
	entt+='          slices: %d\n' % slices  

	rcv=''
	for i in range(num_aperture):
		rcv+='- name: target_e_%.0f \n'%(i)
		rcv+='  side: %s \n' % 'FRONT_AND_BACK'
		rcv+='  per_primitive: %s \n' % 'INCOMING_AND_ABSORBED'
	rcv+='- name: virtual_target_e \n'
	rcv+='  side: %s \n' % 'FRONT'
	rcv+='  per_primitive: %s \n' % 'INCOMING'

	return geom, entt, rcv


