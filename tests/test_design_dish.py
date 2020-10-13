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
		self.casedir='./test-dish-design'
		self.tablefile=self.casedir+'/OELT_Solstice.motab'

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
		nh=9

		dish=Dish(latitude, casedir)
		dish.yaml(dish_radius, dish_foc, rho_refl, slope_error, rec_r, rec_h, rec_x, rec_y, rec_z, rec_grid_w, rec_grid_h, rec_abs, dni=dni, sunshape=sunshape, half_angle_deg=half_angle_deg)
		dish.annual_oelt(dni_des, num_rays, nd, nh, zipfiles=False, gen_vtk=True, plot=False)

		end=time.time()
		print('total time %.2f'%((end-start)/60.), 'min')

	def test_touching(self):

		if os.path.exists(self.tablefile):
			oelt_generated='successful'

		#os.system('rm -rf %s'%self.casedir)

if __name__ == '__main__':
	unittest.main()

