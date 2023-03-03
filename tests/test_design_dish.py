#! /bin/env python3

from __future__ import division
import unittest

import solsticepy
from solsticepy.design_dish import Dish
from solsticepy.output_motab import output_metadata_motab, output_motab, read_motab
import os
import numpy as np
import time

class TestDesignDish(unittest.TestCase):


	def test_dish_parabolia(self):

		start=time.time()
		casedir='./test-dish-parabolia'
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
		dish=Dish(casedir, nproc=1, verbose=True)
		dish.yaml(dish_radius, dish_foc, rho_refl, slope_error, rec_r, rec_x, rec_y, rec_z, rec_grid_r, rec_abs, dni=dni, sunshape=sunshape, half_angle_deg=half_angle_deg)
		eta=dish.get_opt_eff(dni_des, num_rays, zipfiles=False, gen_vtk=True, plot=False)

		end=time.time()
		print('total time %.2f'%((end-start)/60.), 'min')
		self.assertEqual(round(eta, 3), 0.757)		
		#os.system('rm -rf %s'%casedir)

	@unittest.skip("check whether the dish data is shareable")
	def test_dish_multi_facets(self):

		start=time.time()
		casedir='./test-dish-multi-facets'
		dish_radius=None
		dish_foc=13.4
		rho_refl=0.9
		slope_error=1e-3
		rec_r=1.
		rec_x=0.
		rec_y=0.
		rec_z=dish_foc
		rec_grid_r=200
		rec_abs=0.9
		dni=1000
		sunshape='pillbox'
		half_angle_deg=0.2664

		dni_des=950.
		num_rays=int(1e6)



		# provide the vertices and facet surface to the program
		import re 
		fn='dish.obj' #TODO check if this data can be shared or not
		vertices=np.array([])
		faces=np.array([])
		with open(fn) as f:
			for r in f.readlines():
				if 'v' in r:
					v=re.findall(r"[-+]?\d*\.\d+|\d+", r)
					vertices=np.append(vertices, np.r_[float(v[0]), float(v[1]), float(v[2])])
				elif 'f' in r:
					f=re.findall(r"[-+]?\d*\.\d+|\d+", r)
					faces=np.append(faces, np.r_[int(f[0]), int(f[1]), int(f[2]), int(f[3])])
		vertices=vertices.reshape(int(len(vertices)/3),3)
		faces=faces.reshape(int(len(faces)/4),4)
		faces=faces.astype(int)
		multifacets=True
		fct_w=1.17
		fct_h=1.17
		#print(len(faces))

		dish=Dish(casedir, nproc=1, verbose=True)
		dish.yaml(dish_radius, dish_foc, rho_refl, slope_error, rec_r, rec_x, rec_y, rec_z, rec_grid_r, rec_abs, multifacets, vertices, faces, fct_w, fct_h,	dni=dni, sunshape=sunshape, half_angle_deg=half_angle_deg)
		eta=dish.get_opt_eff(dni_des, num_rays, zipfiles=False, gen_vtk=True, plot=False)

		end=time.time()
		print('total time %.2f'%((end-start)/60.), 'min')
		self.assertEqual(round(eta, 3), 0.766)		
		#os.system('rm -rf %s'%casedir)

if __name__ == '__main__':
	unittest.main()

