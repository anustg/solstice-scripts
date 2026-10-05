from __future__ import print_function
import numpy as np
import matplotlib.pyplot as plt

#for python 2:
#from builtins import super

from .data_spectral import SolarSpectrum, MirrorRhoSpectrum
from .cal_layout import multi_aperture_pos
import sys
from .geometry import get_rotation

def yamltransform(pos,rot):
    return "transform: { translation: [%e,%e,%e], rotation: [%e,%e,%e] }" % (*pos,*rot)

def _heliostat_entity(index, position, aim, geometry, target_aligned=False):
    """Write a single-facet heliostat directly in its placed entity."""
    text = '- entity:\n'
    text += '    name: H_%d\n' % index
    text += '    ' + yamltransform(position, [0, 0, 0]) + '\n'
    text += '    children:\n'
    text += '      - name: hst_%d\n' % index
    if not target_aligned:
        text += '        primary: 0\n'
        text += _inline_geometry(_pylon_geometry(), 8)
    text += '        children:\n'
    text += '          - name: pivot\n'
    text += '            zx_pivot:\n'
    if target_aligned:
        text += '              target_aligned: true\n'
    text += '              target:\n'
    text += '                position: [%e, %e, %e]\n' % tuple(aim)
    text += '            children:\n'
    text += '              - name: reflect_surface\n'
    text += '                primary: 1\n'
    if not target_aligned:
        text += '                transform: {rotation: [-90,0,0]}\n'
    text += _inline_geometry(geometry, 16) + '\n'
    return text


def _quadricxy_coefficients(hst_foc, idx, num_hst):
	"""Return quadricxy coefficients for one heliostat.

	Accepted hst_foc formats for shape='quadricxy':
	- scalar or 1D length n: focal length, with fx=fy=focal
	- 2D n x 2: focal lengths [fx, fy]
	- 2D n x 6: coefficients [ax2, ay2, axy, ax, ay, ac]
	"""
	foc = np.asarray(hst_foc, dtype=float)
	if foc.ndim == 0:
		fx = fy = float(foc)
		return 1. / 4. / fx, 1. / 4. / fy, 0., 0., 0., 0.

	if foc.ndim == 1:
		if len(foc) == num_hst:
			fx = fy = foc[idx]
			return 1. / 4. / fx, 1. / 4. / fy, 0., 0., 0., 0.
		if num_hst == 1 and len(foc) == 2:
			fx, fy = foc
			return 1. / 4. / fx, 1. / 4. / fy, 0., 0., 0., 0.
		if num_hst == 1 and len(foc) == 6:
			return tuple(foc)

	if foc.ndim == 2:
		if foc.shape[0] == 1 and num_hst > 1:
			row = foc[0]
		elif foc.shape[0] == num_hst:
			row = foc[idx]
		else:
			raise ValueError("quadricxy hst_foc must have one row or one row per heliostat")

		if foc.shape[1] == 2:
			fx, fy = row
			return 1. / 4. / fx, 1. / 4. / fy, 0., 0., 0., 0.
		if foc.shape[1] == 6:
			return tuple(row)

	raise ValueError("quadricxy hst_foc must be focal, [fx, fy], or [ax2, ay2, axy, ax, ay, ac]")

def _inline_geometry(geometry, indent):
    """Indent a geometry value beneath an entity's geometry key."""
    return ' ' * indent + 'geometry:\n' + ''.join(
        ' ' * (indent + 2) + line + '\n' for line in geometry.splitlines())


def _heliostat_values(values, name, count, total):
    """Broadcast a scalar or validate values in heliostat input order."""
    values = np.asarray(values, dtype=float)
    if values.ndim == 0:
        return np.full(count, float(values))
    if values.ndim != 1 or len(values) not in (count, total):
        raise ValueError('%s must be a scalar or one value per heliostat' % name)
    return values[:count]


def _mirror_material(reflectivity, slope_error, spectral=False):
    """Material mapping embedded in each mirror geometry."""
    rho = '*ref_mirror' if spectral else '%15.8e' % reflectivity
    return ('front:\n'
            '  mirror: {reflectivity: %s, slope_error: %15.8e}\n'
            'back:\n'
            '  matte: {reflectivity: 0.}\n') % (rho, slope_error)


def _mirror_geometry(shape, width, height, material, focal=None, coefficients=None, slices=4):
    """Write an unanchored geometry using the supplied physical curvature."""
    text = '- material:\n' + ''.join('    ' + line + '\n' for line in material.splitlines())
    surface = {'flat': 'plane', 'sphere': 'hemisphere',
               'parabolic-cylinder': 'parabolic-cylinder', 'quadricxy': 'quadricxy'}.get(shape, 'parabol')
    text += '  %s:\n' % surface
    if shape == 'quadricxy':
        for name, value in zip(('ax2', 'ay2', 'axy', 'ax', 'ay', 'ac'), coefficients):
            text += '    %s: %e\n' % (name, value)
    elif shape == 'sphere':
        text += '    radius: %e\n' % (2. * focal)
    elif shape != 'flat':
        text += '    focal: %e\n' % focal
    text += '    clip:\n'
    text += '    - operation: AND\n'
    text += '      vertices: [ [%e, %e], [%e, %e], [%e, %e], [%e, %e] ]\n' % (
        -width/2., -height/2., -width/2., height/2., width/2., height/2., width/2., -height/2.)
    text += '    slices: %d\n' % slices
    return text


def _pylon_geometry():
    return ('- material: *material_black\n  '
            + yamltransform([0, 0, -0.003], [0, 90, 0])
            + '\n  cylinder: {height: 0.001, radius: 0.200, slices: 4}\n')


def _focal_value(values, index):
    values = np.asarray(values, dtype=float)
    return float(values) if values.ndim == 0 else float(values.reshape(-1)[index])


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


def gen_yaml(sun, hst_pos, hst_foc, hst_aims, hst_w, hst_h
		, rho_refl, slope_error, receiver, rec_param, rec_abs
		, outfile_yaml, outfile_recv
		, hemisphere='North', tower_h=0.01, tower_r=0.01,  spectral=False
		, medium=0, one_heliostat=False, cant=False
		, fct_w=0, fct_h=0, fct_gap=0, n_row=0, n_col=0, shape='paraboloid'
		, virtual_target=None, cant_foc=None):

	"""Generate the heliostat field and receiver YAML input files for Solstice ray-tracing simulation.

	1. the sun
	  * `sun` (`Sun` object): parameters relating to the solar source

	2. the field
	  * `hst_pos` (nx3 numpy array): heliostat positions (x, y, z) (first of the 'field' parameters)
	  * `hst_foc` (nx1 numpy array): heliostat focal length or facet focal length
	  * `hst_aims` (nx3 numpy array): heliostat aiming point (ax, ay, az)
	  * `hst_w` (float): heliostat mirror width  (in x direction)
	  * `hst_h` (float): heliostat mirror height (in y direction)
	  * `hst_z` (float): heliostat center height (in z direction)
	  * `rho_refl` (scalar or length n): reflector reflectivity in heliostat input order.
		With spectral=True, the existing mirror spectrum supplies reflectivity instead.
	  * `slope_error` (scalar or length n): reflector surface slope error rms, radians, in heliostat input order
	  * `cant` (bool): True, multi-facets canted heliostats, False, ideally shaped single facet heliostat
	  * `cant_foc` (scalar or length n, optional): actual canting focus for each heliostat.
		Defaults to `hst_foc` for scalar focal lengths. For directional/explicit quadricxy
		curvatures, defaults to the distance from the heliostat to its aiming point.
	  * `fct_w` (float): facet width (if cant==True) 
	  * `fct_h` (float): facet height (if cant==True) 
	  * `fct_gap` (float): facet gaps (if cant==True)  
	  * `n_row` (int): number of rows for the facet arrangement (if cant==True)  
	  * `n_col` (int): number of cols for the facet arrangement (if cant==True)  
	  * `shape` (str): "paraboloid", "quadricxy", "flat", "parabolic-cylinder" or 'sphere' shaped heliostats/facets
		For single-facet "quadricxy", `hst_foc` can be a focal length, [fx, fy], or [ax2, ay2, axy, ax, ay, ac].
		Canted quadricxy facets use the same `hst_foc` formats; use `cant_foc` to set their canting focus.
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

	Heliostat children and geometries are emitted directly inside each entity.
	Mirror curvature uses `hst_foc` directly; mirror materials are embedded per heliostat.

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
	#iyaml+='   ' + yamltransform(pos=[0,0,tower_h*0.5],rot=[0,90,0]) + '\n'  
	iyaml+='    cylinder: {height: %7.3f, radius: %7.3f, slices: %d }\n' % (tower_h, tower_r, slices) 
	iyaml+='\n'
	#
	# Receiver Geometry
	#
	if receiver=='flat':
		geom, rec_entt, rcv = flat_receiver(rec_param, hemisphere, virtual_target=virtual_target)
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
	else:
		hst_x=hst_pos[:,0]
		hst_y=hst_pos[:,1]
		hst_z=hst_pos[:,2]
		aim_x=hst_aims[:,0]
		aim_y=hst_aims[:,1]
		aim_z=hst_aims[:,2]
		num_hst=len(hst_x)
	slices = 4 # slices for the envelop circle
	

	# Validate against the full input field, then select the emitted heliostats.
	reflectivities = _heliostat_values(rho_refl, 'rho_refl', num_hst, len(hst_pos))
	slope_errors = _heliostat_values(slope_error, 'slope_error', num_hst, len(hst_pos))

	# Generate each heliostat using its actual optical parameters. Geometry
	summary = [['x', 'y', 'z', 'ax2', 'ay2', 'axy', 'ax', 'ay', 'ac']
		if shape == 'quadricxy' else ['x', 'y', 'z', 'focal']]
	for i in range(num_hst):
		position = [hst_x[i], hst_y[i], hst_z[i]]
		aim = [aim_x[i], aim_y[i], aim_z[i]]
		material = _mirror_material(reflectivities[i], slope_errors[i], spectral)
		coefficients = _quadricxy_coefficients(hst_foc, i, len(hst_pos)) if shape == 'quadricxy' else None
		focal = None if shape == 'quadricxy' else _focal_value(hst_foc, i)
		geometry = _mirror_geometry(shape, fct_w if cant else hst_w,
			fct_h if cant else hst_h, material, focal, coefficients,
			slices=1 if shape == 'flat' and not cant else 4)
		if cant:
			if cant_foc is not None:
				canting_focus = _focal_value(cant_foc, i)
			elif shape != 'quadricxy' or np.asarray(hst_foc).ndim == 0 or (
				np.asarray(hst_foc).ndim == 1 and len(np.asarray(hst_foc)) == len(hst_pos)):
				canting_focus = _focal_value(hst_foc, i)
			else:
				canting_focus = np.linalg.norm(np.asarray(position) - aim)
			data = heliostat_canted_facets(hst_w, hst_h, fct_w, fct_h,
				fct_gap, n_row, n_col, canting_focus).reshape(n_row, n_col, 4)
			iyaml += '- entity:\n'
			iyaml += '    name: H_%d\n' % i
			iyaml += '    ' + yamltransform([hst_x[i], hst_y[i], hst_z[i]], [0,0,0]) + '\n'
			iyaml += '    children:\n'
			iyaml += '      - name: hst_%d\n' % i
			iyaml += '        primary: 0\n'
			iyaml += _inline_geometry(_pylon_geometry(), 8)
			iyaml += '        children:\n'
			iyaml += '          - name: pivot\n'
			iyaml += '            ' + yamltransform([0,0,0], [0,0,0]) + '\n'
			iyaml += '            zx_pivot: \n'
			iyaml += '              spacing: 0\n'
			iyaml += '              ref_point: [0,0,0]\n'
			iyaml += '              target: {position: [%e, %e, %e]}\n' % (aim_x[i], aim_y[i], aim_z[i])
			iyaml += '            children:\n'
			iyaml += '              - name: facets\n'
			iyaml += '                primary: 0\n'
			iyaml += '                ' + yamltransform([0,0,0], [0,0,0]) + '\n'
			iyaml += _inline_geometry(_mirror_geometry('flat', 2e-5, 2e-5, material, slices=1), 16)
			iyaml += '                children:\n'
			for col in range(n_col):
				for row in range(n_row):
					x, z, rotx, roty = data[row, col]
					iyaml += '                    - name: facet_%d_c%d_r%d\n' % (i, col, row)
					iyaml += '                      primary: 1\n'
					iyaml += '                      ' + yamltransform([x,0,z], [rotx,roty,0]) + '\n'
					iyaml += _inline_geometry(geometry, 22)
		else:
			iyaml += _heliostat_entity(i, position, aim, geometry)
			if shape != 'flat':
				summary.append(position + (list(coefficients) if coefficients is not None else [focal]))
	if not cant and shape != 'flat':
		np.savetxt('heliostat-info-summary.csv', np.asarray(summary), fmt='%s', delimiter=',')


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


	with open(outfile_yaml,'w') as f:
		f.write(iyaml)

	with open(outfile_recv,'w') as f:
		f.write(rcv) 

def gen_yaml_target_aligned(sun, hst_pos, hst_aims, Fx, Fy, hst_w, hst_h
		, rho_refl, slope_error, receiver, rec_param, rec_abs
		, outfile_yaml, outfile_recv
		, hemisphere='North', tower_h=0.01, tower_r=0.01,  spectral=False
		, medium=0, one_heliostat=False):

	"""Generate the heliostat field and receiver YAML input files for Solstice ray-tracing simulation.
	Target-aligned heliostat tracking + paraboloid shape with 2 focal lengthes

	1. the sun
	  * `sun` (`Sun` object): parameters relating to the solar source

	2. the field
	  * `hst_pos` (nx3 numpy array): heliostat positions (x, y, z) (first of the 'field' parameters)
	  * `hst_aims` (nx3 numpy array): heliostat aiming point (ax, ay, az)
      * `Fx` (numpy array, n): tangential focal length (x direction)
      * `Fy` (numpy array, n): sagittal focal length (y direction)
	  * `hst_w` (float): heliostat mirror width  (in x direction)
	  * `hst_h` (float): heliostat mirror height (in y direction)
	  * `hst_z` (float): heliostat center height (in z direction)
	  * `rho_refl` (scalar or length n): reflector reflectivity in heliostat input order.
		With spectral=True, the existing mirror spectrum supplies reflectivity instead.
	  * `slope_error` (scalar or length n): reflector surface slope error rms, radians, in heliostat input order
	  * `tower_h` (float): tower height (m)
	  * `tower_r` (float): tower radius (a cylindrical shape tower) (m)
	  * `target_aligned' (bool): use target_aligned tracking or not (False is using azi-ele tracking), only avaiable for single facet heliostats at the moment
	3. the receiver
	  * `receiver` (str): ``'flat'``, ``'cylinder'``, or ``'stl' or 'multi-aperture'`` (first of the 'receiver' parameters)
	  * `rec_abs` (float): receiver absorptivity
	  * `rec_param` (numpy array or str): each element contains the geometrical parameter of the corresponding receiver.
	4. others
	  * `spectral` (bool): True - simulate the spectral dependent performance (first of the 'other' parameters)
	  * `medium` (float): if the atmosphere is surrounded by non-participant medium, medium=0; otherwise it is the extinction coefficient in m-1
	  * `one_heliosat` (boolean): if `True`, implements ray tracing from just one heliostat.

	Heliostat children and geometries are emitted directly inside each entity.
	Mirror curvature uses `hst_foc` directly; mirror materials are embedded per heliostat.

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
	#iyaml+='   ' + yamltransform(pos=[0,0,tower_h*0.5],rot=[0,90,0]) + '\n'  
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
		slant=((hst_x-aim_x)**2+(hst_y-aim_y)**2+(hst_z-aim_z)**2)**0.5

		slants=np.r_[slant]
	else:
		hst_x=hst_pos[:,0]
		hst_y=hst_pos[:,1]
		hst_z=hst_pos[:,2]
		aim_x=hst_aims[:,0]
		aim_y=hst_aims[:,1]
		aim_z=hst_aims[:,2]
		num_hst=len(hst_x)
		slants=((hst_x-aim_x)**2+(hst_y-aim_y)**2+(hst_z-aim_z)**2)**0.5
	slices = 4 # slices for the envelop circle
	
	# 
	### Section (5)			
	reflectivities = _heliostat_values(rho_refl, 'rho_refl', num_hst, len(hst_pos))
	slope_errors = _heliostat_values(slope_error, 'slope_error', num_hst, len(hst_pos))
	for i in range(num_hst):
		material = _mirror_material(reflectivities[i], slope_errors[i], spectral)
		coefficients = (1. / 4. / Fx[i], 1. / 4. / Fy[i], 0., 0., 0., 0.)
		geometry = _mirror_geometry('quadricxy', hst_w, hst_h, material, coefficients=coefficients)
		iyaml += _heliostat_entity(i, [hst_x[i], hst_y[i], hst_z[i]],
			[aim_x[i], aim_y[i], aim_z[i]], geometry, target_aligned=True)

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


	with open(outfile_yaml,'w') as f:
		f.write(iyaml)

	with open(outfile_recv,'w') as f:
		f.write(rcv) 



def heliostat_canted_facets(hst_w, hst_h, fct_w, fct_h, gap, n_row, n_col, foc):
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

			n0=np.r_[-x/2./foc, 1, -z/2./foc] # facet pointing direction           
			n0=n0/np.linalg.norm(n0)	
			rotx, roty=get_rotation(O, n0)
			
			data=np.append(data, (x, z, rotx, roty))
			#print(x, 0, z, rotx, roty)
	
	data=data.reshape(int(len(data)/4), 4)		
	
	return data


def flat_receiver(rec_param, hemisphere='North', virtual_target=None):
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
	if len(rec_param)>8:
		tilt_y=rec_param[8]
	else:
		tilt_y=0

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

		entt+='    ' + yamltransform(pos=[x, y , z],rot=[-90.-tilt, tilt_y, 0]) + '\n'
	else:
		entt+='    ' + yamltransform(pos=[x, y , z],rot=[90.+tilt, -tilt_y, 0]) + '\n'

	entt+='    geometry: *%s\n' % 'target_g'


	# CREATE a virtual target entity (primary = 0).
	entt+='\n- entity:\n'
	entt+='    name: virtual_target_e\n'
	entt+='    primary: 0\n'
	if virtual_target is None:
		radius=np.sqrt(rec_w**2+rec_h**2)/2.*1.5
		entt+='    transform: { translation: [%e,%e,%e]}\n' % (x, y, z)
		entt+='    geometry: \n'
		entt+='      - material: *%s\n' % 'material_virtual'
		entt+='        sphere: \n'
		entt+='          radius: %s\n' % radius
		entt+='          slices: %d\n' % 20
	else:
		virtual_w=float(virtual_target['width'])
		virtual_h=float(virtual_target['height'])
		virtual_slices=int(virtual_target.get('slices', 200))
		virtual_pos=virtual_target.get('position', (x, y, z))
		virtual_tilt=float(virtual_target.get('tilt', tilt))
		virtual_tilt_y=float(virtual_target.get('tilt_y', tilt_y))
		if hemisphere=='North':
			virtual_rot=[-90.-virtual_tilt, virtual_tilt_y, 0]
		else:
			virtual_rot=[90.+virtual_tilt, -virtual_tilt_y, 0]
		entt+='    ' + yamltransform(pos=virtual_pos, rot=virtual_rot) + '\n'
		entt+='    geometry: \n'
		entt+='      - material: *%s\n' % 'material_virtual'
		entt+='        plane: \n'
		entt+='          clip: \n'
		entt+='          - operation: AND \n'
		entt+='            vertices: [ [%e, %e], [%e, %e], [%e, %e], [%e, %e] ] \n' % (-virtual_w*0.5, -virtual_h*0.5, -virtual_w*0.5, virtual_h*0.5, virtual_w*0.5, virtual_h*0.5, virtual_w*0.5, -virtual_h*0.5)
		entt+='          slices: %d\n' % virtual_slices


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

	radius=np.sqrt(rec_r**2+rec_h**2/4.)*1.5
	entt+='\n- entity:\n'
	entt+='    name: virtual_target_e\n'
	entt+='    primary: 0\n'
	entt+='    transform: { translation: [%e,%e,%e]}\n' % (x, y, z)
	entt+='    geometry: \n' 
	entt+='      - material: *%s\n' % 'material_virtual' 
	entt+='        sphere: \n'
	entt+='          radius: %s\n' % radius   
	entt+='          slices: %d\n' % 20 

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
	entt+='    transform: { translation: [0,0,%e]}\n' % (vir_z)
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


