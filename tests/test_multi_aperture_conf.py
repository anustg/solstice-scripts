#! /bin/env python3

from __future__ import division, print_function
import unittest

import solsticepy
from solsticepy.design_multi_aperture import *
import os
import numpy as np
import time

import matplotlib.pyplot as plt

class TestMultiApertureConfiguration(unittest.TestCase):

	def setup(self):
		
		self.w_rcv=10.
		self.h_rcv=10.
		self.H_tower=200.
		

	def test_3_apertures(self):
		n=3
		gamma=np.r_[60., 90., 135., 180., 360.]
		W_rcv=np.array([])
		H_rcv=np.array([]
		for i in range(n):
			W_rcv=np.append(W_rcv, self.w_rcv)
			H_rcv=np.append(H_rcv, self.h_rcv)

		OMEGA=np.array([[60.,  90., 120.],
						[45.,  90., 135.],
						[22.5, 90., 157.5 ],
						[0.,   90., 180.],
						[-90., 30., 150.]])
		LV=np.r_[1, 3, 2]
		Z=np.r_[200., 190., 180.]

		for j in range(len(gamma)):
			g=gamma[j]
			MAC=MultiApertureConfiguration(n, g, self.H_tower, W_rcv, H_rcv)
			for i in range(n):
				omega_i=MAC.get_angular_pos(i)
				xi, yi, zi=MAC.get_cood_pos(i)
				lv_i = MAC.get_lv_idx(i)
				print('')
				print('Aperture ', i)
				print('Omega ', omega_i)
				print('Level ', lv_i)
				print('Pos   ', xi, yi, zi)

				self.assertEqual(omega_i, OMEGA[j,i])
				self.assertEqual(lv_i, LV[i])
				self.assertEqual(zi, Z[i])
				

if __name__ == '__main__':
	unittest.main()

