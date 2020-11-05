#! /bin/env python3

from __future__ import division
import unittest

from solsticepy.cal_layout import *
import os
import numpy as np

class TestLayout(unittest.TestCase):
	def setUp(self):

		latitude=34.
		self.num_hst=22640
		width=10.
		height=10.
		hst_z=5., 
		towerheight=250.
		R1=80.
		fb=0.6
		dsep=0.
		field='polar'
		savedir='.'	
		plot=False

		self.pos_and_aim, Nzones, Nrows_zone=radial_stagger(latitude, self.num_hst, width, height, hst_z, towerheight, R1, fb, dsep, field, savedir, plot)

	def test_touching(self):
		num=len(self.pos_and_aim)-2
		self.assertEqual(num, self.num_hst)
		#os.system('rm *.csv')


if __name__ == '__main__':
	unittest.main()

