#! /bin/env python3

from __future__ import division
import unittest

import solsticepy
from solsticepy.design_crs import CRS
from solsticepy.input import Parameters
from solsticepy.output_motab import output_matadata_motab, output_motab
import os
import numpy as np
import time

class TestDesignCRS(unittest.TestCase):
	def setUp(self):

		start=time.time()
		self.casedir='./test-crs-design'
		self.tablefile=self.casedir+'/OELT_Solstice.motab'
		if os.path.exists(self.tablefile):    
		    print('')
		    print('Load exsiting OELT')

		else:

			pm=Parameters()
			pm.Q_in_rcv=565e6
			pm.nd=5
			pm.nh=22
			pm.H_tower=250.
			pm.H_rcv=24.
			pm.W_rcv=24.
			pm.dependent_par()
			pm.saveparam(self.casedir)
			print(pm.fb)
			print(pm.H_tower)
			crs=CRS(latitude=pm.lat, casedir=self.casedir)   
			weafile='/home/yewang/solartherm-master/SolarTherm/Data/Weather/gen3p3_Daggett_TMY3.motab'
			crs.heliostatfield(field=pm.field_type, hst_rho=pm.rho_helio, slope=pm.slope_error, hst_w=pm.W_helio, hst_h=pm.H_helio, tower_h=pm.H_tower, tower_r=pm.R_tower, hst_z=pm.Z_helio, num_hst=pm.n_helios, R1=pm.R1, fb=pm.fb, dsep=pm.dsep)

			crs.receiversystem(receiver=pm.rcv_type, rec_w=float(pm.W_rcv), rec_h=float(pm.H_rcv), rec_x=float(pm.X_rcv), rec_y=float(pm.Y_rcv), rec_z=float(pm.Z_rcv), rec_tilt=float(pm.tilt_rcv), rec_grid=int(pm.n_H_rcv), rec_abs=float(pm.alpha_rcv))

			crs.yaml(dni=900,sunshape=pm.sunshape,csr=pm.crs,half_angle_deg=pm.half_angle_deg,std_dev=pm.std_dev)

			oelt, A_land=crs.field_design_annual(dni_des=900., num_rays=int(5e6), nd=pm.nd, nh=pm.nh, weafile=weafile, method=1, Q_in_des=pm.Q_in_rcv, n_helios=None, zipfiles=False, gen_vtk=False, plot=False)


			if (A_land==0):    
				self.tablefile=None
			else:                                                
				A_helio=pm.H_helio*pm.W_helio
				output_matadata_motab(table=oelt, field_type=pm.field_type, aiming='single', n_helios=crs.n_helios, A_helio=A_helio, eff_design=crs.eff_des, H_rcv=pm.H_rcv, W_rcv=pm.W_rcv, H_tower=pm.H_tower, Q_in_rcv=pm.Q_in_rcv, A_land=A_land, savedir=self.tablefile)


		end=time.time()
		print('total time %.2f'%((end-start)/60.), 'min')

	def test_touching(self):
		if os.path.exists(self.tablefile):
			test='successful'
		self.assertEqual(test,'successful')
		#os.system('rm -rf %s'%self.casedir)


if __name__ == '__main__':
	unittest.main()

