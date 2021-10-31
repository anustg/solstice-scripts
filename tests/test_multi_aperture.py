#! /bin/env python3

from __future__ import division
import unittest

import solsticepy
from solsticepy.design_crs import CRS
from solsticepy.input import Parameters
from solsticepy.output_motab import output_metadata_motab_multi_aperture, read_motab
from solsticepy.design_multi_aperture import MultiApertureConfiguration
import os
import numpy as np
import time

class TestMultiAperture(unittest.TestCase):
	def setUp(self):

		start=time.time()
		self.casedir='./test-multi-aperture'
		self.tablefile=self.casedir+'/OELT_Solstice.motab'
		if os.path.exists(self.tablefile):    
			print('')
			print('Load exsiting OELT')
			self.newcase=False

		else:
			self.newcase=True

			pm=Parameters()
			pm.Q_in_rcv=56e6
			pm.n_row_oelt=3
			pm.n_col_oelt=6
			pm.H_tower=120.
			pm.H_rcv=[12.,12.,12.]
			pm.W_rcv=[12.,12.,12.]
			pm.Z_rcv=[pm.H_tower,pm.H_tower,pm.H_tower]
			pm.fb=0.5
			pm.R1=50.
			pm.field_type='multi-aperture'
			pm.rcv_type='multi-aperture'	
			pm.num_aperture=3
			pm.gamma=180. # deg
			pm.n_W_rcv=10
			pm.n_H_rcv=10	
			pm.helio_rho=0.9
			pm.helio_soil=1.
			pm.helio_sf_ratio=1.

			pm.dependent_par()
			pm.saveparam(self.casedir)

			crs=CRS(latitude=pm.lat, casedir=self.casedir, nproc=1, verbose=True)   
			weafile='../example/demo_TMY3_weather.motab'
			mac=MultiApertureConfiguration(n=pm.num_aperture, gamma=pm.gamma, H_tower=pm.H_tower, R_tower=pm.R_tower, W_rcv=pm.W_rcv, H_rcv=pm.H_rcv, parallel=False)
						
			crs.receiversystem(receiver=pm.rcv_type, rec_w=pm.W_rcv, rec_h=pm.H_rcv, rec_x=pm.X_rcv, rec_y=pm.Y_rcv, rec_z=pm.Z_rcv, rec_tilt=pm.tilt_rcv, rec_grid_w=int(pm.n_W_rcv), rec_grid_h=int(pm.n_H_rcv), rec_abs=pm.alpha_rcv, mac=mac)

			crs.heliostatfield(field=pm.field_type, hst_rho=pm.helio_refl, slope=pm.slope_error, hst_w=pm.W_helio, hst_h=pm.H_helio, tower_h=pm.H_tower, tower_r=pm.R_tower, hst_z=pm.Z_helio, num_hst=pm.n_helios, R1=pm.R1, fb=pm.fb, dsep=pm.dsep, mac=mac)

			crs.yaml(dni=900,sunshape=pm.sunshape,csr=pm.csr,half_angle_deg=pm.half_angle_deg,std_dev=pm.std_dev)

			self.oelt, A_land=crs.field_design_annual(dni_des=900., num_rays=int(1e6), nd=pm.n_row_oelt, nh=pm.n_col_oelt, weafile=weafile, method=1, Q_in_des=pm.Q_in_rcv, n_helios=None, zipfiles=False, gen_vtk=False, plot=False)

			self.n_helios=crs.n_helios
			self.eff_des=crs.eff_des
			self.eff_annual=crs.eff_annual
			self.num_aperture=pm.num_aperture

			if (A_land==0):    
				self.tablefile=None
			else:                                                
				A_helio=pm.H_helio*pm.W_helio

				output_metadata_motab_multi_aperture(
							TABLE=self.oelt, 
							eff_design = crs.eff_des,
							eff_annual = crs.eff_annual,
							A_land     = A_land, 
							H_tower    = pm.H_tower, 
							A_helio    = A_helio, 
							n_helios_total = crs.n_helios, 
							Q_in_rcv_total = crs.Q_in_rcv, 
							num_aperture   = pm.num_aperture, 
							Q_in_rcv   = crs.Q_in_rcv_i, 
							n_helios   = crs.n_helios_i, 
							H_rcv      = pm.H_rcv, 
							W_rcv      = pm.W_rcv,  
							Z_rcv      = pm.Z_rcv,  
							mac        = mac,
							savedir    = self.tablefile)

		end=time.time()
		print('total time %.2f'%((end-start)/60.), 'min')

	def test_touching(self):

		if self.newcase:
			eta_max=np.max(self.oelt[self.num_aperture][3:,3:].astype(float))
		else:

			self.eff_des, self.eff_annual, A_land, H_tower, A_helio, self.n_helios, Q_in_rcv, num_aperture, radial_distance, Q_in_rcv_i, n_helios_i, H_rcv_i, W_rcv_i, Z_rcv_i, solar_hour, declination, OELT=read_motab(self.tablefile, multi_aperture=True)

			eta_max=np.max(OELT[num_aperture])


		print(self.n_helios, self.eff_des, self.eff_annual)
		if os.path.exists(self.tablefile):
			oelt_generated='successful'
		self.assertEqual(oelt_generated,'successful')
		self.assertTrue(abs(self.n_helios-756)/756. < 0.01)
		self.assertTrue(abs(self.eff_des- 0.741)/0.741 < 0.01)
		self.assertTrue(abs(self.eff_annual-0.672)/0.672 < 0.01)
		#os.system('rm -rf %s'%self.casedir)

if __name__ == '__main__':
	unittest.main()

