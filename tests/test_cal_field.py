#! /bin/env python3

from __future__ import division
import unittest

from solsticepy.cal_field import *
from solsticepy.cal_layout import radial_stagger
import os
import numpy as np

class TestFieldPF(unittest.TestCase):
	def setUp(self):

		# The cal_field script calculates preliminary field performance using cosine
		# factor and outputs the calculated cosine factors to a vtk file for viewing
		# (heliostats can be coloured according to cosine factor)
		towerheight=250.
		width=10.
		height=10.

		pos_and_aim=radial_stagger(latitude=34., num_hst=22640., width=width, height=height, hst_z=5., towerheight=towerheight, R1=80., fb=0.6, dsep=0., field='polar', savedir='.', plot=False)

		pos=pos_and_aim[2:,:3].astype(float)
		aim=pos_and_aim[2:,4:].astype(float)
		azimuth=np.r_[0.]
		zenith=np.r_[12.]
		field=FieldPF(np.r_[0,1,0])
		sun_vec=field.get_solar_vector(azimuth, zenith)
		norms=field.get_normals(towerheight=towerheight, hstpos=pos, sun_vec=sun_vec)
		#field.heliostat(10, 8)
		COORD, TRI, ele, nc=field.mesh_heliostat_field(width=width, height=height, normals=norms, hstpos=pos)
		cos=field.get_cosine(hst_norms=norms, sun_vec=sun_vec)
		self.savedir='./field.vtk'
		COS=np.repeat(cos, ele)
		DATA={'cos':COS}
		NORMS=np.repeat(norms, ele, axis=0)
		gen_vtk(self.savedir, COORD.T, TRI, NORMS, True, DATA)


	def test_touching(self):
		if os.path.exists(self.savedir):
			successed=1
		self.assertEqual(successed,1)
		os.system('rm *.vtk')
		os.system('rm *.csv')


if __name__ == '__main__':
	unittest.main()

