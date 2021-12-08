#! /bin/env python3

from __future__ import division
import unittest

from solsticepy.cal_layout import *
from solsticepy.design_multi_aperture import MultiApertureConfiguration
import os
import numpy as np

class TestLayout(unittest.TestCase):
	def setUp(self):
		self.latitude=34.
		self.num_hst=22640
		self.width=10.
		self.height=10.
		self.hst_z=5., 
		self.towerheight=250.
		self.towerdiameter=30.
		self.R1=80.
		self.fb=0.6
		self.dsep=0.
		self.savedir='.'	
		self.plot=False


	def test_polarfield(self):
		field='polar'
		pos_and_aim, Nzones, Nrows_zone, TTRow=radial_stagger(self.latitude, self.num_hst, self.width, self.height, self.hst_z, self.towerheight, self.R1, self.fb, self.dsep, field, savedir=self.savedir, plot=self.plot)
		num=len(pos_and_aim)-2
		self.assertEqual(num, self.num_hst)
		#os.system('rm *.csv')

	def test_multiaperture(self):
		field='multi-aperture'
		num_aperture=3
		ang_rang=180.
		rec_w=[12., 12., 12.]
		mac=MultiApertureConfiguration(n=num_aperture, gamma=ang_rang, H_tower=self.towerheight, R_tower=self.towerdiameter/2., W_rcv=rec_w, H_rcv=rec_w, parallel=False)
		
		rec_z=[self.towerheight, self.towerheight, self.towerheight]
		pos_and_aim, Nzones, Nrows_zone, TTRow=radial_stagger(
				latitude=self.latitude, 
				num_hst=self.num_hst, 
				width=self.width, 
				height=self.height, 
				hst_z=self.hst_z, 
				towerheight=self.towerheight, 
				R1=self.R1, 
				fb=self.fb, 
				dsep=self.dsep, 
				field=field, 
				mac=mac,
				savedir=self.savedir, 
				verbose=False,
				plot=self.plot,
				plt_aiming=None)

		num=len(pos_and_aim)-2
		print(num)
		self.assertEqual(num, self.num_hst)


if __name__ == '__main__':
	unittest.main()

