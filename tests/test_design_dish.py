#! /bin/env python3

from __future__ import division
import unittest

import solsticepy
from solsticepy.design_dish import Dish
from solsticepy.output_motab import output_matadata_motab, output_motab, read_motab
import os
import numpy as np
import time

class TestDesignCRS(unittest.TestCase):
	def setUp(self):

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
		self.eta=dish.get_opt_eff(dni_des, num_rays, zipfiles=False, gen_vtk=True, plot=False)

		end=time.time()
		print('total time %.2f'%((end-start)/60.), 'min')

	def test_touching(self):

		self.assertEqual(round(self.eta, 5), 0.75717)		

		#os.system('rm -rf %s'%self.casedir)

if __name__ == '__main__':
	unittest.main()

