#! /bin/env python3

"""Single-heliostat quadricxy simulation at equinox solar noon.

Select the heliostat configuration with HELIOSTAT_MODE=single, canted, or both.
"""

from __future__ import division

import os
import unittest

import numpy as np

import solsticepy
from solsticepy.master import Master


# Select "single", "canted", or "both".
HELIOSTAT_MODE = os.environ.get("HELIOSTAT_MODE", "canted").lower()


class TestHelioshape(unittest.TestCase):
	def test_single_heliostat_at_equinox_noon(self):
		if HELIOSTAT_MODE == "both":
			modes = ("single", "canted")
		else:
			self.assertIn(HELIOSTAT_MODE, ("single", "canted"))
			modes = (HELIOSTAT_MODE,)

		for mode in modes:
			with self.subTest(mode=mode):
				eta, performance_hst = self.run_case(mode)       
				self.assertTrue(np.isfinite(eta.nominal_value))
				#self.assertEqual(len(performance_hst), 1)

	def run_case(self, mode):
		latitude = 34.0
		azimuth = 270.0  # Solstice convention: south at solar noon.
		elevation = 90.0 - latitude  # Equinox solar noon: declination = 0 deg.
		num_rays = 200000

		dni = 1000.0
		sun = solsticepy.Sun(
			dni=dni,
			sunshape="pillbox",
			half_angle_deg=0.2664,
		)

		hst_w = 10.0
		hst_h = 10.0
		rho_refl = 0.95
		slope_error = 2.0e-3
		tower_h = 80.0
		tower_r = 0.01

		# One heliostat north of the tower, aimed at the receiver centre.
		hst_pos = np.array([[0.0, 100.0, 3.0]])
		hst_aims = np.array([[0.0, 0.0, tower_h]])
		slant = np.linalg.norm(hst_aims[0] - hst_pos[0])

		receiver = "flat"
		rec_param = np.r_[8.0, 6.0, 100, 100, 0.0, 0.0, tower_h, 0.0]
		rec_abs = 0.9

		casedir = "./test_helioshape_%s" % mode
		master = Master(casedir)
		outfile_yaml = master.in_case(casedir, "input.yaml")
		outfile_recv = master.in_case(casedir, "input-rcv.yaml")

		common = dict(
			sun=sun,
			hst_pos=hst_pos,
			hst_aims=hst_aims,
			hst_w=hst_w,
			hst_h=hst_h,
			rho_refl=rho_refl,
			slope_error=slope_error,
			receiver=receiver,
			rec_param=rec_param,
			rec_abs=rec_abs,
			outfile_yaml=outfile_yaml,
			outfile_recv=outfile_recv,
			hemisphere="North",
			tower_h=tower_h,
			tower_r=tower_r,
			spectral=False,
			medium=0,
			one_heliostat=True,
			shape="quadricxy",
		)

		if mode == "single":
			# One surface with independent focal lengths in x and y.
			solsticepy.gen_yaml(
				hst_foc=np.array([slant, slant]),
				cant=False,
				**common
			)
		else:
			# A 2 x 2 canted assembly. Each facet uses independent fx and fy.
			facet_gap = 0.1
			facet_w = (hst_w - facet_gap) / 2.0
			facet_h = (hst_h - facet_gap) / 2.0
			bands = np.array([[slant, slant, slant, slant]])
			solsticepy.gen_yaml(
				hst_foc=np.array([slant]),
				cant=True,
				bands=bands,
				fct_w=facet_w,
				fct_h=facet_h,
				fct_gap=facet_gap,
				n_row=2,
				n_col=2,
				**common
			)

		return master.run(
			azimuth,
			elevation,
			num_rays,
			rho_refl,
			dni,
			folder=master.in_case(casedir, "sunpos_equinox_noon"),
			gen_vtk=False,
			verbose=False,
		)


if __name__ == "__main__":
	unittest.main()
