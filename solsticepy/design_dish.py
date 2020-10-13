import os
import sys
import time
import numpy as np
import matplotlib.pyplot as plt
from uncertainties import ufloat
from scipy.interpolate import interp1d,interp2d
import matplotlib.cm as cm
from scipy.optimize import curve_fit

from .process_raw import *
from .cal_sun import *
from .gen_yaml import gen_yaml, Sun
from .output_motab import output_motab
from .master import *


class Dish:
	'''
	The Solar Dish model includes three parts: 
	the sun, the parabolic concentrator and the receiver.
	'''

	def __init__(self, latitude, casedir):
		'''
		Arguements:
			casedir : str, the directory of the case 
		'''
		self.casedir=casedir

		if not os.path.exists(casedir):
			os.makedirs(casedir)
		self.latitude=latitude
		self.sun=SunPosition()
		self.master=Master(casedir)

	def receiversystem(self, rec_r=0., rec_h=0., rec_x=0., rec_y=0., rec_z=100., rec_grid_w=10, rec_grid_h=10, rec_abs=1.):

		'''
		It is assumed that the dish system has a cylindrical cavity receiver 
		that is defined by the radius and height of the cylinder

		Arguements:
		    (1) rec_r     : float, radius of the cylinder (m) 
		    (2) rec_h     : float, height of the cylinder (m)
		    (3) rec_x     : float, x location of the center of the aperture (m)
		    (4) rec_y     : float, y location of the center of the aperture (m)
		    (5) rec_z     : float, z location of the center of the aperture (m)
		    (6) rec_grid_w  :   int, number of elements in the horizontal(x)/circumferential direction
		    (7) rec_grid_h  :   int, number of elements in the vertical(z) direction
		    (8) rec_abs   : float, receiver surface absorptivity, e.g. 0.9
		'''

		geom=''
		geom+='- geometry: &%s\n' % 'target_g'
		geom+='  - material: *%s\n' % 'material_target'
		geom+='    cylinder: \n'
		geom+='      height: %s\n'%rec_h 
		geom+='      radius: %s\n'%rec_r 
		geom+='      slices: %d\n' % rec_grid_w
		geom+='      stacks: %d\n' % rec_grid_h 
		geom+='\n'

		# CREATE a receiver entity from "target_g" geometry (primary = 0)
		entt=''
		entt+='\n- entity:\n'
		entt+='    name: target_e\n'
		entt+='    primary: 0\n'
		entt+='    transform: { translation: %s, rotation: %s }\n' % ([rec_x, rec_y, rec_z], [0., 0., 0.]) 
		entt+='    geometry: *%s\n' % 'target_g'

		# CREATE a virtual target entity from "target_g" geometry (primary = 0)
		Vsize=100.
		pts = [ [-rec_h*Vsize, -rec_h*Vsize], [-rec_h*Vsize, rec_h*Vsize], [rec_h*Vsize, rec_h*Vsize], [rec_h*Vsize,-rec_h*Vsize] ]
		slices = 4
		entt+='\n- entity:\n'
		entt+='    name: virtual_target_e\n'
		entt+='    primary: 0\n'

		entt+='    transform: { translation: %s, rotation: %s }\n' % ([rec_x, rec_y, rec_z+rec_h/2.+1], [-180., 0, 0])

		entt+='    geometry: \n' 
		entt+='      - material: *%s\n' % 'material_virtual' 
		entt+='        plane: \n'
		entt+='          clip: \n'    
		entt+='          - operation: AND \n'
		entt+='            vertices: %s\n' % pts
		entt+='          slices: %d\n' % 4 

		rcv=''
		rcv+='- name: target_e \n' 
		rcv+='  side: %s \n' % 'FRONT_AND_BACK'
		rcv+='  per_primitive: %s \n' % 'INCOMING_AND_ABSORBED'
		rcv+='- name: virtual_target_e \n'
		rcv+='  side: %s \n' % 'FRONT'
		rcv+='  per_primitive: %s \n' % 'INCOMING'

		return geom, entt, rcv



	def yaml(self, dish_radius, dish_foc, rho_refl, slope_error, rec_r, rec_h, rec_x, rec_y, rec_z, rec_grid_w, rec_grid_h, rec_abs, dni=1000, sunshape=None, csr=0.01, half_angle_deg=0.2664, std_dev=0.2):
		'''
		Generate YAML files for the Solstice simulation	
		The vertice of the dish concentrator is (0,0,0)

		Arguements:
			(1) dish_radius: float, radius of the dish concentrator (m)
			(2) dish_foc   : float, focal length of the dish concentrator (m)
		    (3) rho_refl   : float, reflectivity of the mirror surface
			(4) slope_error: float, surfac slope error of the mirrors (rad)
		    (5) rec_r     : float, radius of the cylinder (m) 
		    (6) rec_h     : float, height of the cylinder (m)
		    (7) rec_x     : float, x location of the center of the aperture (m)
		    (8) rec_y     : float, y location of the center of the aperture (m)
		    (9) rec_z     : float, z location of the center of the aperture (m)
		    (10) rec_grid_w  :   int, number of elements in the horizontal(x)/circumferential direction
		    (11) rec_grid_h  :   int, number of elements in the vertical(z) direction
		    (12) rec_abs   : float, receiver surface absorptivity, e.g. 0.9
			(13) `dni`: Direct normal irradance (W/m2)
			(14) `sunshape`: Sunshape: can be None, ``'pillbox'``,``'gaussian'`` or ``'buie'``
			(15) `half_angle_deg`: sun angular size (in DEGREES, half-angle) (ONLY in case of ``'pillbox'``)
			(16) `csr`: circumsolar ratio (ONLY in case of ``'buie'``)
			(17) `std_dev`: standard deviation of the angular dsn ratio (ONLY in case of ``'gaussian'``

		'''

		self.rho_refl=rho_refl
		outfile_yaml = self.master.in_case(self.casedir, 'input.yaml')
		outfile_recv = self.master.in_case(self.casedir, 'input-rcv.yaml')
		'''
		geom, rec_entt, rcv=self.receiversystem(rec_r, rec_h, rec_x, rec_y, rec_z, rec_grid_w, rec_grid_h, rec_abs)

		iyaml=''
		sun = Sun(sunshape=sunshape, csr=csr, half_angle_deg=half_angle_deg, std_dev=std_dev)

		iyaml += "- sun: &sun\n"		
		iyaml += "    dni: %s\n"%dni

		if sunshape is not None:
			if sunshape=='pillbox':
				iyaml += "    pillbox: {half_angle: %6.4f}\n" % (half_angle_deg)   
			elif self.sunshape=='buie':
				iyaml += "    buie: {csr: %6.4f}\n" % (csr)  
			elif self.sunshape=='gaussian':
				iyaml += "    gaussian: {std_dev: %6.4f}\n" % (std_dev) 
		iyaml+='\n'

		#
		#    Materials
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
		
		Vsize=100.
		pts = [ [-rec_h*Vsize, -rec_h*Vsize], [-rec_h*Vsize, rec_h*Vsize], [rec_h*Vsize, rec_h*Vsize], [rec_h*Vsize,-rec_h*Vsize] ]
		slices = 4
		#
		#    Template
		#
		iyaml+='- template: &self_oriented_dish\n'
		iyaml+='    name: so_parabol\n'
		iyaml+='    transform: { translation: [0,0,50], rotation: [0,0,90] }\n'   
		iyaml+='    x_pivot:\n'
		iyaml+='      ref_point: [0,0,0]\n' 
		iyaml+='      target: { sun: *sun }\n'
		iyaml+='    children:\n'
		iyaml+='      - name: parabol\n'
		iyaml+='        primary: 1\n'
		iyaml+='        geometry:\n'
		iyaml+='        - material: *%s\n' % 'material_mirror' 
		iyaml+='          parabol: \n'
		iyaml+='            focal: %s\n' % dish_foc
		iyaml+='            clip: \n'  
		iyaml+='            - operation: AND \n'
		iyaml+='              circle:\n'
		iyaml+='                radius: %s\n'% dish_radius 
		iyaml+='                center: [0,0]\n' 
		iyaml+='            slices: %d\n' % 50 
		iyaml+='      - name: target_surface\n'
		iyaml+='        transform: { translation: [0, 0, 50] }\n'
		iyaml+='        primary: 0\n'
		iyaml+='        geometry:\n'
		iyaml+='        - material: *%s\n' % 'material_target'
		iyaml+='          cylinder: \n'
		iyaml+='            height: %s\n'%rec_h 
		iyaml+='            radius: %s\n'%rec_r 
		iyaml+='            slices: %d\n' % rec_grid_w
		iyaml+='            stacks: %d\n' % rec_grid_h 
		iyaml+='      - name: virtual_target\n'
		iyaml+='        transform: { translation: [0, 0, 52] }\n'
		iyaml+='        primary: 0\n'
		iyaml+='        geometry:\n'
		iyaml+='        - material: *%s\n' % 'material_virtual'
		iyaml+='          plane: \n'
		iyaml+='            clip: \n'    
		iyaml+='            - operation: AND \n'
		iyaml+='              vertices: %s\n' % pts
		iyaml+='            slices: %d\n' % 4 
		iyaml+='\n'



		#
		#    Entities
		#
		# Dish concentrator
		iyaml+='- entity:\n'
		iyaml+='    name: %s\n' % 'dish_e'
		iyaml+='    transform: { translation: %s, rotation: %s }\n' % ([0, 0, 0], [0, 0, 0]) 
		iyaml+='    children: [*%s ]\n' % 'self_oriented_dish'    


		rcv=''
		rcv+='- name: dish_e.so_parabol.target_surface \n' 
		rcv+='  side: %s \n' % 'FRONT_AND_BACK'
		rcv+='  per_primitive: %s \n' % 'INCOMING_AND_ABSORBED'
		rcv+='- name: dish_e.so_parabol.virtual_target \n'
		rcv+='  side: %s \n' % 'FRONT'
		rcv+='  per_primitive: %s \n' % 'INCOMING'
		'''
		iyaml='''- sun: &sun
    dni: 1
    spectrum: [{wavelength: 1, data: 1}]
    
- material: &specular
      mirror: { reflectivity: 1, slope_error: 0 }
      
- material: &black
      matte: { reflectivity: 0 }
                  
- template: &self_oriented_parabol
    name: "so_parabol"
    transform: { translation: [0, 0, 0], rotation: [0, 0, 90] }
    x_pivot:
      ref_point: [0, 0, 0]
      target: { sun: *sun }
    children:
      - name: "parabol"
        primary: 1
        geometry:
        - material: *specular
          parabol:
            focal: 4
            clip:
            - operation: AND
              vertices: [[-5.0, -5.0], [-5.0, 5.0],
               [5.0, 5.0], [5.0, -5.0]]
      - name: "small_square"
        transform: { translation: [0, 0, 4] }
        primary: 0
        geometry:
        - material:  *black
          plane:
            clip:
              - operation: AND
                vertices: [[-0.50, -0.50], [-0.50, 0.50],
                 [0.50, 0.50], [0.50, -0.50]]
      
- entity:
    name: "reflector"
    transform: { rotation: [0, 0, 0], translation: [0, 0, 0] }
    children: [ *self_oriented_parabol ]'''

		rcv='''- { name: "reflector.so_parabol.small_square", side: FRONT_AND_BACK }'''

		with open(outfile_yaml,'w') as f:
			f.write(iyaml)

		with open(outfile_recv,'w') as f:
			f.write(rcv) 



	def annual_oelt(self, dni_des, num_rays, nd, nh, zipfiles=False, gen_vtk=False, plot=False):
		'''
		Annual performance of a known field
		'''  

		oelt=self.master.run_annual_dish(nd=nd, nh=nh, latitude=self.latitude, num_rays=num_rays, rho_mirror=1., dni=dni_des, gen_vtk=gen_vtk)


		np.savetxt(self.casedir+'/lookup_table.csv', oelt, fmt='%s', delimiter=',')


		designfolder=self.casedir+'/des_point'
		day=self.sun.days(21, 'Mar')
		dec=self.sun.declination(day)
		hra=0. # solar noon
		zen=self.sun.zenith(self.latitude, dec, hra)
		azi=self.sun.azimuth(self.latitude, zen, dec, hra)     
		azi_des, ele_des=self.sun.convert_convention('solstice', azi, zen) 
		azi_des-=90.

		sys.stderr.write("\n"+green('Design Point: \n'))		
		efficiency_total=self.master.run(azi_des, ele_des, num_rays, self.rho_refl, dni_des, folder=designfolder, gen_vtk=gen_vtk, printresult=False, system='dish')
		self.eff_des=efficiency_total.n

		return oelt
		

	def dni_TMY(self, weafile, nd, nh):
		'''
		Argument:

		weafile : str, directory of the weather file .motab
		# col0 - time (s)
		# col1 - GHI
		# col2 -DNI (W/m2)
		# col3 -DHI 
			...
		'''
		with open(weafile) as f:
			content=f.read().splitlines()
		f.close()
	
		lines=len(content)
		seconds=np.array([])
		dni=np.array([])
		for i in range(2, lines-1):
			l=content[i].split(",") 
			seconds=np.append(seconds, l[0])         
			dni=np.append(dni, l[2])
		seconds=seconds.astype(float)
		days=seconds/3600/24
		wea_dec=np.array([])
		for d in days:
			wea_dec=np.append(wea_dec, self.sun.declination(d)) #deg

		wea_hra=((seconds/3600.)%24-12.)*15. #deg
		wea_dni=dni.astype(float)

		hra_lim=180.*(float(nh)/float(nh-1))
		dec_lim=23.45*(float(nd)/float(nd-1))
		hra_bin=np.linspace(-hra_lim, hra_lim, nh+1) 
		dec_bin=np.linspace(-dec_lim, dec_lim, nd+1)
		bins=np.array([hra_bin, dec_bin])


		dni_weight, xbins, ybins=np.histogram2d(wea_hra, wea_dec, bins, weights=wea_dni)

		return dni_weight.T



if __name__=='__main__':
	start=time.time()
	latitude=32.4
	casedir='./test-dish-design'
	tablefile=casedir+'/OELT_Solstice.motab'
	dish_radius=10.
	dish_foc=50.
	rho_refl=1.
	slope_error=2e-3
	rec_r=3.
	rec_h=3.
	rec_x=0.
	rec_y=0.
	rec_z=dish_foc
	rec_grid_w=10
	rec_grid_h=10
	rec_abs=1.
	dni=1000
	sunshape='pillbox'
	half_angle_deg=0.2664

	dni_des=950.
	num_rays=int(1e6)
	nd=5
	nh=5

	dish=Dish(latitude, casedir)
	dish.yaml(dish_radius, dish_foc, rho_refl, slope_error, rec_r, rec_h, rec_x, rec_y, rec_z, rec_grid_w, rec_grid_h, rec_abs, dni=dni, sunshape=sunshape, half_angle_deg=half_angle_deg)
	dish.annual_oelt(dni_des, num_rays, nd, nh, zipfiles=False, gen_vtk=False, plot=False)

	end=time.time()
	print('total time %.2f'%((end-start)/60.), 'min')



