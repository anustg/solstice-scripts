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

	def __init__(self, casedir):
		'''
		Arguements:
			casedir : str, the directory of the case 
		'''
		self.casedir=casedir

		if not os.path.exists(casedir):
			os.makedirs(casedir)
		self.master=Master(casedir)

	def yaml(self, dish_radius, dish_foc, rho_refl, slope_error, rec_r, rec_x, rec_y, rec_z, rec_grid_r, rec_abs, dni=1000, sunshape=None, csr=0.01, half_angle_deg=0.2664, std_dev=0.2):
		'''
		Generate YAML files for the Solstice simulation	
		The vertice of the dish concentrator is (0,0,0)

		Arguements:
			(1) dish_radius: float, radius of the dish concentrator (m)
			(2) dish_foc   : float, focal length of the dish concentrator (m)
		    (3) rho_refl   : float, reflectivity of the mirror surface
			(4) slope_error: float, surfac slope error of the mirrors (rad)
		    (5) rec_r     : float, radius of the receiver aperature (m) 
		    (6) rec_x     : float, x location of the center of the aperture (m)
		    (7) rec_y     : float, y location of the center of the aperture (m)
		    (8) rec_z     : float, z location of the center of the aperture (m)
		    (9) rec_grid_r  :   int, number of elements in the horizontal(x)/circumferential direction
		    (10) rec_abs   : float, receiver surface absorptivity, e.g. 0.9
			(11) `dni`: Direct normal irradance (W/m2)
			(12) `sunshape`: Sunshape: can be None, ``'pillbox'``,``'gaussian'`` or ``'buie'``
			(13) `half_angle_deg`: sun angular size (in DEGREES, half-angle) (ONLY in case of ``'pillbox'``)
			(14) `csr`: circumsolar ratio (ONLY in case of ``'buie'``)
			(15) `std_dev`: standard deviation of the angular dsn ratio (ONLY in case of ``'gaussian'``
		'''
		
		self.rho_refl=rho_refl
		outfile_yaml = self.master.in_case(self.casedir, 'input.yaml')
		outfile_recv = self.master.in_case(self.casedir, 'input-rcv.yaml')

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
		
		#
		#    Template
		#
		iyaml+='- template: &self_oriented_dish\n'
		iyaml+='    name: so_parabol\n'
		iyaml+='    transform: { translation: [0,0, 0], rotation: [0,0,0] }\n'   
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
		iyaml+='      - name: receiver_aperture\n'
		iyaml+='        transform: { translation: [0, 0, %s] }\n'%rec_z
		iyaml+='        primary: 0\n'
		iyaml+='        geometry:\n'
		iyaml+='        - material: *%s\n' % 'material_target'
		iyaml+='          plane: \n'
		iyaml+='            clip: \n'  
		iyaml+='            - operation: AND \n'
		iyaml+='              circle:\n'
		iyaml+='                radius: %s\n'% rec_r
		iyaml+='            slices: %d\n' % rec_grid_r
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
		rcv+='- name: dish_e.so_parabol.receiver_aperture \n' 
		rcv+='  side: %s \n' % 'FRONT_AND_BACK'
		rcv+='  per_primitive: %s \n' % 'INCOMING_AND_ABSORBED'

		with open(outfile_yaml,'w') as f:
			f.write(iyaml)

		with open(outfile_recv,'w') as f:
			f.write(rcv) 


	def get_opt_eff(self, dni_des, num_rays, zipfiles=False, gen_vtk=False, plot=False):
		'''
		Assume the optical efficiency of a dish system is constant 
		at different sun positions
		'''  
		azi_des=0
		ele_des=90

		sys.stderr.write("\n"+green('Optical efficiency: \n'))		
		efficiency_total=self.master.run(azi_des, ele_des, num_rays, self.rho_refl, dni_des, folder=self.casedir, gen_vtk=gen_vtk, printresult=True, system='dish')
		eff_des=efficiency_total.n

		return eff_des
		

if __name__=='__main__':
	start=time.time()
	casedir='./test-dish-design'

	dish_radius=10.
	dish_foc=10.
	rho_refl=0.9
	slope_error=2e-3
	rec_r=1.
	rec_x=0.
	rec_y=0.
	rec_z=dish_foc
	rec_grid_r=50
	rec_abs=0.9
	dni=1000
	sunshape='pillbox'
	half_angle_deg=0.2664

	dni_des=950.
	num_rays=int(1e6)
	nd=5
	nh=9

	dish=Dish(casedir)
	dish.yaml(dish_radius, dish_foc, rho_refl, slope_error, rec_r, rec_x, rec_y, rec_z, rec_grid_r, rec_abs, dni=dni, sunshape=sunshape, half_angle_deg=half_angle_deg)
	eta=dish.get_opt_eff(dni_des, num_rays, zipfiles=False, gen_vtk=True, plot=False)
	print('total efficiency:', eta)

	end=time.time()
	print('total time %.2f'%((end-start)/60.), 'min')



