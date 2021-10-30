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

	def setUp(self):
		
		self.w_rcv=10.
		self.h_rcv=10.
		self.H_tower=200.
		self.R_tower=15.

	def test_3_apertures(self):
		n=3
		gamma=np.r_[60., 90., 135., 180., 360.]
		W_rcv=np.array([])
		H_rcv=np.array([])
		for i in range(n):
			W_rcv=np.append(W_rcv, self.w_rcv)
			H_rcv=np.append(H_rcv, self.h_rcv)

		OMEGA=np.array([[60.,  90., 120.],
						[45.,  90., 135.],
						[22.5, 90., 157.5 ],
						[0.,   90., 180.],
						[-90., 30., 150.]])
		LV=np.r_[1, 3, 2]
		Z=np.r_[200., 180., 190.]

		for j in range(len(gamma)):
			g=gamma[j]
			MAC=MultiApertureConfiguration(n, g, self.H_tower, self.R_tower, W_rcv, H_rcv)
			for i in range(n):
				omega_i=MAC.get_angular_pos(i)
				xi, yi, zi=MAC.get_cood_pos(i)
				lv_i = MAC.get_lv_index(i)
				print('')
				print('Aperture %s/%s'%(i,n))
				print('Omega ', omega_i)
				print('Level ', lv_i)
				print('Pos   ', xi, yi, zi)

				self.assertEqual(omega_i, OMEGA[j,i])
				self.assertEqual(lv_i, LV[i])
				self.assertEqual(zi, Z[i])

			IDX=np.r_[0, 2, 1]
			for lv in range(1, n+1):
				idx=MAC.get_i_index(lv)
				print('lv', lv, ': aperture', idx)
				self.assertEqual(idx, IDX[lv-1])


	def test_6_apertures(self):
		n=6
		gamma=np.r_[60., 90., 135., 180., 360.]
		W_rcv=np.array([])
		#H_rcv=np.array([])
		for i in range(n):
			W_rcv=np.append(W_rcv, self.w_rcv)
			#H_rcv=np.append(H_rcv, self.h_rcv)
		H_rcv=np.r_[10., 12., 15., 15., 12.,10.]


		OMEGA=np.array([[60.,  72.,  84.,  96,    108.,  120.],
						[45.,  63.,  81.,  99.,   117.,  135.],
						[22.5, 49.5, 76.5, 103.5, 130.5, 157.5 ],
						[0.,   36.,  72.,  108.,  144.,  180.],
						[-90., -30., 30.,  90.,   150.,  210., ]])

		LV=np.r_[1, 3, 5, 6, 4, 2]
		Z=np.r_[200., 179., 153.5, 138.5, 167., 190.]

		for j in range(len(gamma)):
			g=gamma[j]
			MAC=MultiApertureConfiguration(n, g, self.H_tower, self.R_tower, W_rcv, H_rcv)
			for i in range(n):
				omega_i=MAC.get_angular_pos(i)
				xi, yi, zi=MAC.get_cood_pos(i)
				lv_i = MAC.get_lv_index(i)
				print('')
				print('Aperture %s/%s'%(i,n))
				print('Omega ', omega_i)
				print('Level ', lv_i)
				print('Pos   ', xi, yi, zi)

				self.assertEqual(omega_i, OMEGA[j,i])
				self.assertEqual(lv_i, LV[i])
				self.assertEqual(zi, Z[i])

			IDX=np.r_[0, 5, 1, 4, 2, 3]
			for lv in range(1, n+1):
				idx=MAC.get_i_index(lv)
				print('lv', lv, ': aperture', idx)
				self.assertEqual(idx, IDX[lv-1])

				
				

if __name__ == '__main__':
	unittest.main()

