from __future__ import print_function
import numpy as np
import matplotlib.pyplot as plt

#for python 2:
#from builtins import super

from .data_spectral import SolarSpectrum, MirrorRhoSpectrum
from .cal_layout import multi_aperture_pos
import sys

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
		`half_angle_deg`: sun angular size (in DEGREES, half-angle) (ONLY in case of ``'pillbox'``)
		`csr`: circumsolar ratio (ONLY in case of ``'buie'``)
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


def gen_yaml(sun, hst_pos, hst_foc, hst_aims, hst_w, hst_h
		, rho_refl, slope_error, cant, bands, receiver, rec_param, rec_abs
		, outfile_yaml, outfile_recv
		, hemisphere='North', tower_h=0.01, tower_r=0.01,  spectral=False
		, medium=0, one_heliostat=False
		, fct_w=0, fct_h=0, fct_gap=0, n_row=0, n_col=0, shape='curved'):
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
	  * `cant` (bool): True, multi-facets canted heliostats, False, ideally shaped single facet heliostat
	  * `bands` (nparray or None): a 2D numpy array ((band range, focal length)) that specifies 
						the distance range and focal length for each canting band, it is <= than this
						, or None for slant range canting
	  * `fct_w` (float): facet width (if cant==True) 
	  * `fct_h` (float): facet height (if cant==True) 
	  * `fct_gap` (float): facet gaps (if cant==True)  
	  * `n_row` (int): number of rows for the facet arrangement (if cant==True)  
	  * `n_col` (int): number of cols for the facet arrangement (if cant==True)  
	  * `shape` (str): "parabol" or "flat" or "parabolic-cylinder"shaped heliostats/facets  
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
	  * if ``receiver == 'cylinder'``: np.array([diameter, height, grid_circ, grid_h, x, y, z, tilt angle (deg)])
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
		    iyaml+='  - {wavelength: %s, data: %s }\n' % (I_sun[i][0],I_sun[i][1])  
		i = len(I_sun)-1
		iyaml+='  - {wavelength: %s, data: %s }\n' % (I_sun[i][0],I_sun[i][1]) 
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
		iyaml+='- atmosphere: {extinction: %s}\n'%medium 
		iyaml+='\n'

		   
	# 
	### Section (3)
	# set the materials
	# (gathering media)
	# occultant material, mirror specular material, receiver material, virtual target
	#------------------------------
	#
	# CREATE an occultant material
	r_f = 0. # front
	r_b = 0. # and back reflectivity
	iyaml+='- material: &%s\n' % 'material_black'
	iyaml+='   front:\n'
	iyaml+='     matte: {reflectivity: %6.4f }\n' % r_f    
	iyaml+='   back:\n'
	iyaml+='     matte: {reflectivity: %6.4f }\n' % r_b
	iyaml+='\n'
	#
	# CREATE a specular material
	r_f= rho_refl # front
	r_b = 0.      # and back reflectivity
	iyaml+='- material: &%s\n' % 'material_mirror'
	iyaml+='   front:\n'
	if spectral:
		iyaml+='     mirror: {reflectivity: *%s, slope_error: %15.8e }\n' % ('ref_mirror', slope_error ) 
	else:
		iyaml+='     mirror: {reflectivity: %6.4f, slope_error: %15.8e }\n' % (r_f, slope_error) 

	iyaml+='   back:\n'
	iyaml+='     matte: {reflectivity: %6.4f }\n' % r_b 
	iyaml+='\n'
	#
	# CREATE a material for the target
	r_f = 1.-rec_abs # front
	r_b = 1.-rec_abs # and back reflectivity
	iyaml+='- material: &%s\n' % 'material_target'
	iyaml+='   front:\n'
	iyaml+='     matte: {reflectivity: %6.4f }\n' % r_f    
	iyaml+='   back:\n'
	iyaml+='     matte: {reflectivity: %6.4f }\n' % r_b
	iyaml+='\n'
	#
	# CREATE a virtual material for the calculation of spillage
	iyaml+='- material: &%s\n' % 'material_virtual'
	iyaml+='   virtual:\n'
	iyaml+='\n'


	# 
	### Section (4) & (5)
	# set the geometries
	# (gathering shapes and materials)
	# the tower, the receiver, the heliostat
	# 
	# set the templates
	# (programming objects gathering geometries or pivot and geometries)
	#------------------------------

	#
	# Tower Geometry
	# (cylindrical shape)
	#
	slices = 10 # slices for the envelop circle
	iyaml+='- geometry: &%s\n' % 'tower_g' 
	iyaml+='  - material: *%s\n' % 'material_black' 
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
		hst_x=np.r_[hst_pos[0, 0]]
		hst_y=np.r_[hst_pos[0, 1]]
		hst_z=np.r_[hst_pos[0, 2]]
		aim_x=np.r_[hst_aims[0, 0]] 
		aim_y=np.r_[hst_aims[0 ,1]]
		aim_z=np.r_[hst_aims[0, 2]]
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
	pts_hst = [ [-hst_w*0.5, -hst_h*0.5], [-hst_w*0.5, hst_h*0.5], [hst_w*0.5, hst_h*0.5], [hst_w*0.5,-hst_h*0.5] ]

	# CREATE the pylon "pylon_g" geometry cylindrical shape
	h_pyl = hst_z[0] # pylon height
	r_pyl = 0.2 # pylon radius
	slices = 4 # slices for the envelop circle
	iyaml+='- geometry: &%s\n' % 'pylon_g' 
	iyaml+='  - material: *%s\n' % 'material_black' 
	iyaml+='    transform: { translation: %s, rotation: %s }\n' % ([0, 0, -h_pyl/2.], [0, 0, 0]) 
	iyaml+='    cylinder: {height: %7.3f, radius: %7.3f, slices: %d }\n' % (h_pyl*0.8,r_pyl,slices) 
	#   

	if bands.any()==None:
		if one_heliostat:
			bands=np.array([hst_foc, hst_foc])
			bands=bands.reshape(1,2)
		else:
			dist=(hst_w+hst_h)/2.
			min_foc=np.min(hst_foc)
			max_foc=np.max(hst_foc)
			bands=np.arange(min_foc+dist, max_foc+dist, dist)
			bands=np.append(bands, bands)
			bands=bands.reshape(2, int(len(bands)/2))
			bands=bands.T 


	if cant==True: # multi facets

		iyaml+="""
- geometry: &facet_g_s
  - material: *material_mirror
    plane:
      clip: 
      - operation: AND 
        vertices: [[-0.00001, -0.00001], [-0.00001, 0.00001], [0.000015, 0.00001], [0.00001, -0.00001]]
      slices: 1\n\n"""

		if shape=='flat':
			iyaml+="- geometry: &facet_g\n"
			iyaml+="  - material: *material_mirror\n"
			iyaml+='    plane: \n'
			iyaml+="      clip:\n"
			iyaml+="      - operation: AND\n" 
			iyaml+="        vertices: [[%s, %s], [%s, %s], [%s, %s], [%s, %s]]\n"%(-fct_w/2., -fct_h/2.,-fct_w/2., fct_h/2., fct_w/2., fct_h/2., fct_w/2., -fct_h/2.)
			iyaml+="      slices: 1\n\n"	



		elif shape=='parabolic-cylinder': # curved facets 
			for i in range(len(bands)):
				foc=bands[i,1]				
				iyaml+="- geometry: &facet_g_band_%s\n"%i
				iyaml+="  - material: *material_mirror\n"
				iyaml+="    parabolic-cylinder:\n"
				iyaml+="      focal: %s\n"%foc
				iyaml+="      clip:\n"
				iyaml+="      - operation: AND\n" 
				iyaml+="        vertices: [[%s, %s], [%s, %s], [%s, %s], [%s, %s]]\n"%(-fct_w/2., -fct_h/2.,-fct_w/2., fct_h/2., fct_w/2., fct_h/2., fct_w/2., -fct_h/2.)
				iyaml+="      slices: 4\n\n"
				
		else: # parabol curved facets 
			for i in range(len(bands)):
				foc=bands[i,1]				
				iyaml+="- geometry: &facet_g_band_%s\n"%i
				iyaml+="  - material: *material_mirror\n"
				iyaml+="    parabol:\n"
				iyaml+="      focal: %s\n"%foc
				iyaml+="      clip:\n"
				iyaml+="      - operation: AND\n" 
				iyaml+="        vertices: [[%s, %s], [%s, %s], [%s, %s], [%s, %s]]\n"%(-fct_w/2., -fct_h/2.,-fct_w/2., fct_h/2., fct_w/2., fct_h/2., fct_w/2., -fct_h/2.)
				iyaml+="      slices: 4\n\n"
				


		for b in range(len(bands)):
			iyaml+="- template: &facets_t_band_%s\n"%b
			iyaml+="    name: facets\n"
			iyaml+="    primary: 0\n"
			iyaml+="    transform: {translation: [0,0,0], rotation: [0,0,0]}\n"
			iyaml+="    geometry: *facet_g_s\n"
			iyaml+="    children:\n"

			foc=bands[b,1]		
			#print('foc', foc)
			data=heliostat_canted_facets(hst_w, hst_h, fct_w, fct_h, fct_gap, n_row, n_col, foc, shape)
			fct_x=data[:,0].reshape(n_row, n_col)
			fct_z=data[:,1].reshape(n_row, n_col)
			rotx=data[:,2].reshape(n_row, n_col)
			roty=data[:,3].reshape(n_row, n_col)

			for j in range(n_col):
				for i in range(n_row):			
					iyaml+='        - name: facet_%s_c%s_r%s\n'%(b, j, i)
					iyaml+='          primary: 1\n'
					iyaml+='          transform: {translation: [%.4f,0,%.4f], rotation: [%s,%s,0]}\n'%(fct_x[i,j], fct_z[i,j], rotx[i,j]-90, roty[i,j])
					if shape=='flat':
						iyaml+='          geometry: *facet_g\n'
					else:
						iyaml+='          geometry: *facet_g_band_%s\n'%b


		# heliostat entities from the template
		for i in range(num_hst):

			iyaml+='- template: &hst_t_%s\n'%i
			iyaml+='    name: hst_%s\n'%i
			iyaml+='    primary: 0\n'
			iyaml+='    geometry: *pylon_g\n'
			iyaml+='    children:\n'     
			iyaml+='      - name: pivot\n'
			iyaml+='        transform: { translation: [0,0, 0], rotation: [0,0,0] }\n'
			iyaml+='        zx_pivot: \n'
			iyaml+='          spacing: 0\n'
			iyaml+='          ref_point: [0,0,0]\n'
			iyaml+='          target: {position: [%s, %s, %s]}\n' % (aim_x[i],aim_y[i],aim_z[i]) 
			foc=hst_foc[i]
			idx=np.argmin(abs(foc-bands[:,0]))
			if foc-bands[idx, 0]>0:
				idx+=1
			#idx= np.where(bands[foc<=bands][0]==bands)[0][0]
			iyaml+='        children: [ *facets_t_band_%s ]\n\n'%idx

			iyaml+='- entity:\n'
			iyaml+='    name: H_%s\n'%i
			iyaml+='    transform: { translation: [%s, %s, %s], rotation: [0, 0, 0] }\n'%(hst_x[i], hst_y[i], hst_z[i])
			iyaml+='    children: [ *hst_t_%s ]\n'%i

	else: # single facet
		if shape=='flat':
			iyaml+="- geometry: &hst_g\n"
			iyaml+="  - material: *material_mirror\n"
			iyaml+='    plane: \n'
			iyaml+="      clip:\n"
			iyaml+="      - operation: AND\n" 
			iyaml+="        vertices: [[%s, %s], [%s, %s], [%s, %s], [%s, %s]]\n"%(-fct_w/2., -fct_h/2.,-fct_w/2., fct_h/2., fct_w/2., fct_h/2., fct_w/2., -fct_h/2.)
			iyaml+="      slices: 1\n\n"


			for i in range(num_hst):
			# CREATE the heliostat templates   
				name_hst_t = 'hst_t'
				iyaml+='- template: &%s\n' % name_hst_t 
				name_hst_n = 'hst_'
				iyaml+='    name: %s\n' % name_hst_n 
				iyaml+='    primary: 0\n'   
				iyaml+='    geometry: *pylon_g\n'
				iyaml+='    children: \n' 
				iyaml+='    - name: pivot\n'
				iyaml+='      zx_pivot: {target: {position: %s}} \n' % ([aim_x[i],aim_y[i],aim_z[i]]) 
				iyaml+='      children: \n'
				iyaml+='      - name: reflect_surface\n'
				iyaml+='        primary: 1\n'
				iyaml+='        transform: {rotation: [-90,0,0]} \n' 

				#idx= np.where(bands[foc<=bands][0]==bands)[0][0]
				name_hst_g = 'hst_g'
				iyaml+='        geometry: *%s\n\n' % name_hst_g 


				name_e ='H_'+str(i)
				name_hst_t = 'hst_t'
				iyaml+='\n- entity:\n'
				iyaml+='    name: %s\n' % name_e
				iyaml+='    transform: { translation: %s, rotation: %s }\n' % ([hst_x[i], hst_y[i], hst_z[i]], [0, 0, 0]) 
				iyaml+='    children: [ *%s ]\n' % name_hst_t 	

		else:

			for i in range(len(bands)):
				name_hst_g = 'hst_g_band_'+str(i)
				iyaml+='- geometry: &%s\n' % name_hst_g 
				iyaml+='  - material: *%s\n' % 'material_mirror' 
				iyaml+='    parabol: \n'
				iyaml+='      focal: %s\n' % bands[i,1]
				iyaml+='      clip: \n'  
				iyaml+='      - operation: AND \n'
				iyaml+='        vertices: %s\n' % pts_hst
				iyaml+='      slices: %d\n\n' % slices 
			summary=np.array(['x','y','z', 'band foc'])
			for i in range(num_hst):
				# CREATE the heliostat templates   
				name_hst_t = 'hst_t_'+str(i)
				iyaml+='- template: &%s\n' % name_hst_t 
				name_hst_n = 'hst_'+ str(i)
				iyaml+='    name: %s\n' % name_hst_n 
				iyaml+='    primary: 0\n'   
				iyaml+='    geometry: *pylon_g\n'
				iyaml+='    children: \n' 
				iyaml+='    - name: pivot\n'
				iyaml+='      zx_pivot: {target: {position: %s}} \n' % ([aim_x[i],aim_y[i],aim_z[i]]) 
				iyaml+='      children: \n'
				iyaml+='      - name: reflect_surface\n'
				iyaml+='        primary: 1\n'
				iyaml+='        transform: {rotation: [-90,0,0]} \n' 
				foc=hst_foc[i]
				if len(bands)==1:
					idx=0
				else:
					idx=np.argmin(abs(foc-bands[:,0]))
					if foc-bands[idx, 0]>0:
						idx+=1
				#idx= np.where(bands[foc<=bands][0]==bands)[0][0]
				name_hst_g = 'hst_g_band_'+str(idx)
				iyaml+='        geometry: *%s\n\n' % name_hst_g 
				summary=np.append(summary, (hst_x[i], hst_y[i], hst_z[i], bands[idx, 1]))
			summary=summary.reshape(int(len(summary)/4), 4)
			np.savetxt('heliostat-info-summary.csv', summary, fmt='%s', delimiter=',')


			for i in range(num_hst):
				name_e ='H_'+str(i)
				name_hst_t = 'hst_t_'+str(i)
				iyaml+='\n- entity:\n'
				iyaml+='    name: %s\n' % name_e
				iyaml+='    transform: { translation: %s, rotation: %s }\n' % ([hst_x[i], hst_y[i], hst_z[i]], [0, 0, 0]) 
				iyaml+='    children: [ *%s ]\n' % name_hst_t 

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
	iyaml+='    transform: { translation: %s, rotation: %s }\n' % ([0, -tower_r, tower_h*0.5], [0, 0, 0]) 
	iyaml+='    geometry: *%s\n' % 'tower_g'    
	#

	with open(outfile_yaml,'w') as f:
		f.write(iyaml)

	with open(outfile_recv,'w') as f:
		f.write(rcv) 

def heliostat_canted_facets(hst_w, hst_h, fct_w, fct_h, gap, n_row, n_col, foc, shape='curved'):
	"""
	Ideally on-axis canted (4fy=x^2+z^2)

	Arguments:
  	* `hst_w` (float): width of the heliostat
  	* `hst_h` (float): height of the heliostat
  	* `fct_w` (float): width of a facet
  	* `fct_h` (float): height of a facet
  	* `gap` (float): gap between facets (assumes equal gaps in x/y directions)
  	* `n_row` (int): number of rows for the facets arrangement 
  	* `n_col` (int): number of columns for the facets arrangement
  	* `foc` (float): focal length 
  	* `shape` (str): 'flat' or 'curved' facets 
	"""
	data=np.array([])
	for i in range(n_row):
		for j in range(n_col):
			z=-hst_h/2.+fct_h/2.+i*(gap+fct_h)
			x=-hst_w/2.+fct_w/2.+j*(gap+fct_w)

			# the heliostat center is (0, 0, 0)
			# the focal point is (0, foc, 0)
			# facet location
			O=np.r_[x, 0, z] 
			OH=np.r_[0,1,0] # original hst norm 
			OX=np.r_[1,0,0]	

			n0=np.r_[-x/2./foc, 1, -z/2./foc] # facet pointing direction           
			n0=n0/np.linalg.norm(n0)			

			n1=np.cross(OH, OX)
			n2=np.cross(n0, OX)	

			rotx=np.arccos(np.dot(n1, n2)/(np.linalg.norm(n1)*np.linalg.norm(n2)))*180./np.pi
			if z>0:
				rotx=-rotx

			v_proj=np.cross(n2, (np.cross(OH, n2)/np.linalg.norm(n2)))/np.linalg.norm(n2)
			norm=np.linalg.norm(v_proj)*np.linalg.norm(n0)

			if norm<1e-12:
				cosy=np.dot(n0, np.r_[0,0,1])
				roty=-np.arccos(cosy)*180./np.pi
			else:
				cosy=np.dot(v_proj, n0)/norm
				if 1-abs(cosy)<1e-12:
					roty=0
				else:
					roty=np.arccos(cosy)*180./np.pi

			if x>0:
				roty=-roty
			elif x<0 and z==0:
				roty=-roty
			

			data=np.append(data, (x, z, rotx, roty))

			#print(x, 0, z, n0[0], n0[1], n0[2])
	
	data=data.reshape(int(len(data)/4), 4)		
	
	return data


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
	rec_w=rec_param[0]
	rec_h=rec_param[1]
	slices=rec_param[2] # it assumes equal number of slices in x and y directions
	x=rec_param[4]
	y=rec_param[5]
	z=rec_param[6]
	tilt=rec_param[7]
	# receiver tilt angle:
	# 0 is vertical
	# the standby posiion of a plane in solstice is normal points to the +z axis
	# rotation anagle, positive is anti-clockwise

	geom=''
	pts=[ [-rec_w*0.5, -rec_h*0.5], [-rec_w*0.5, rec_h*0.5], [rec_w*0.5, rec_h*0.5], [rec_w*0.5,-rec_h*0.5] ]

	geom+='- geometry: &%s\n' % 'target_g'
	geom+='  - material: *%s\n' % 'material_target'
	geom+='    plane: \n'
	geom+='      clip: \n' 
	geom+='      - operation: AND \n'
	geom+='        vertices: %s\n' % pts
	geom+='      slices: %d\n' % slices 
	geom+='\n'

	# CREATE a receiver entity from "target_g" geometry (primary = 0)
	entt=''
	entt+='\n- entity:\n'
	entt+='    name: target_e\n'
	entt+='    primary: 0\n'
	if hemisphere=='North':
		entt+='    transform: { translation: %s, rotation: %s }\n' % ([x, y, z], [-90.-tilt, 0, 0]) 
	else:
		entt+='    transform: { translation: %s, rotation: %s }\n' % ([x, y, z], [90.+tilt, 0, 0]) 
	entt+='    geometry: *%s\n' % 'target_g'

	# CREATE a virtual target entity from "target_g" geometry (primary = 0)
	pts = [ [-rec_w*10., -rec_h*10.], [-rec_w*10., rec_h*10.], [rec_w*10., rec_h*10.], [rec_w*10.,-rec_h*10.] ]
	slices = 4
	entt+='\n- entity:\n'
	entt+='    name: virtual_target_e\n'
	entt+='    primary: 0\n'
	if hemisphere=='North':
		entt+='    transform: { translation: %s, rotation: %s }\n' % ([x, y+5., z], [-90.-tilt, 0, 0])
	else:
		entt+='    transform: { translation: %s, rotation: %s }\n' % ([x, y-5., z], [90.+tilt, 0, 0])
	entt+='    geometry: \n' 
	entt+='      - material: *%s\n' % 'material_virtual' 
	entt+='        plane: \n'
	entt+='          clip: \n'    
	entt+='          - operation: AND \n'
	entt+='            vertices: %s\n' % pts
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
	rec_r=rec_param[0]/2.
	rec_h=rec_param[1]
	slices=rec_param[2] # number of elements in the circumferetial direction
	stacks=rec_param[3] # number of elements in the vertical direction
	x=rec_param[4]
	y=rec_param[5]
	z=rec_param[6]

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

	entt+='    transform: { translation: %s, rotation: %s }\n' % ([x, y, z], [0., 0., 0.]) 

	entt+='    geometry: *%s\n' % 'target_g'

	# CREATE a virtual target entity from "target_g" geometry (primary = 0)
	Vsize=100.
	pts = [ [-rec_h*Vsize, -rec_h*Vsize], [-rec_h*Vsize, rec_h*Vsize], [rec_h*Vsize, rec_h*Vsize], [rec_h*Vsize,-rec_h*Vsize] ]
	slices = 4
	entt+='\n- entity:\n'
	entt+='    name: virtual_target_e\n'
	entt+='    primary: 0\n'

	entt+='    transform: { translation: %s, rotation: %s }\n' % ([x, y, z-rec_h/2.-1], [-180., 0, 0])

	entt+='    geometry: \n' 
	entt+='      - material: *%s\n' % 'material_virtual' 
	entt+='        plane: \n'
	entt+='          clip: \n'    
	entt+='          - operation: AND \n'
	entt+='            vertices: %s\n' % pts
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

	rec_w=rec_param[0].astype(float) # for creating the virtual target
	rec_h=rec_param[1].astype(float)
	stlfile=rec_param[2] # directory of the stl file
	x=rec_param[3].astype(float)
	y=rec_param[4].astype(float)
	z=rec_param[5].astype(float)
	tilt=rec_param[6].astype(float) # need to figure out the initial mesh orientation

	# CREATE a receiver entity from a STL file 
	entt=''
	entt+='\n- entity:\n'
	entt+='    name: STL_receiver_e\n'
	entt+='    primary: 0\n'
	if hemisphere=='North':

		entt+='    transform: { translation: %s, rotation: %s }\n' % ([x, y, z], [-90.-tilt, 0, 0]) 
	else:
		# if it is the mesh model of the bladed receiver at CSIRO
		entt+='    transform: { translation: %s, rotation: %s }\n' % ([x, y, z], [180.+tilt, 0, 0]) 
	entt+='    geometry:\n'
	entt+='    - material: *material_target\n'
	entt+='      transform: {translation: [0, 0, 0], rotation: [0, 0, 0]}\n'
	entt+="      stl : {path: %s }  \n"%(stlfile)


	# CREATE a virtual target entity from "target_g" geometry (primary = 0)
	pts = [ [-rec_w*10., -rec_h*10.], [-rec_w*10., rec_h*10.], [rec_w*10., rec_h*10.], [rec_w*10.,-rec_h*10.] ]
	slices = 4
	entt+='\n- entity:\n'
	entt+='    name: virtual_target_e\n'
	entt+='    primary: 0\n'
	if hemisphere=='North':
		entt+='    transform: { translation: %s, rotation: %s }\n' % ([x, y-5., z], [-90.-tilt, 0, 0])
	else:
		entt+='    transform: { translation: %s, rotation: %s }\n' % ([x, y+5., z], [90.+tilt, 0, 0])
	entt+='    geometry: \n' 
	entt+='      - material: *%s\n' % 'material_virtual' 
	entt+='        plane: \n'
	entt+='          clip: \n'    
	entt+='          - operation: AND \n'
	entt+='            vertices: %s\n' % pts
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
	rec_grid_w=rec_param[2]
	rec_grid_h=rec_param[3]

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

		pts=[ [-rec_w[i]*0.5, -rec_h[i]*0.5], [-rec_w[i]*0.5, rec_h[i]*0.5], [rec_w[i]*0.5, rec_h[i]*0.5], [rec_w[i]*0.5,-rec_h[i]*0.5] ]

		geom+='- geometry: &%s\n' % 'target_g_%.0f\n'%(i)
		geom+='  - material: *%s\n' % 'material_target'
		geom+='    plane: \n'
		geom+='      clip: \n' 
		geom+='      - operation: AND \n'
		geom+='        vertices: %s\n' % pts
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
			entt+='    transform: { translation: %s, rotation: %s }\n' % ([xc, yc, zc], [-90.-rec_tilt, 90.-ang_pos,0]) 
		else:
			entt+='    transform: { translation: %s, rotation: %s }\n' % ([-xc, -yc, zc], [90.+rec_tilt, 90.-ang_pos,0]) 
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



#------------------------------

