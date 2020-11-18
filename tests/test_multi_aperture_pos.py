#! /bin/env python3

from __future__ import division
import unittest

import solsticepy
from solsticepy.cal_layout import multi_aperture_pos
import os
import numpy as np
import time

import matplotlib.pyplot as plt

class TestMultiAperturePos(unittest.TestCase):

	def test_touching(self):

		omega, xc, yc=multi_aperture_pos(rec_w=10., gamma=180., n=3, i=2)
		self.assertEqual(omega, 180)

		omega, xc, yc=multi_aperture_pos(rec_w=10., gamma=240., n=3, i=2)
		self.assertEqual(omega, 210)

		# to plot some figures
		'''
		GAMMA=np.r_[60., 120., 180., 270., 360.]
		N=np.r_[2,3,4,5,6]
		
		for g in GAMMA:
			for n in N:
				OMEGA=np.array([])
				for i in range(n):
					omega=multi_aperture_angular_pos(gamma=g, n=n, i=i)
					OMEGA=np.append(OMEGA, omega)

				ax = plt.subplot(111, projection='polar')
				ax.plot(OMEGA*np.pi/180., np.ones(n), 'o')
				plt.title('$\\gamma$=%s, n=%s'%(g, n))
				plt.show()
				plt.close()

		'''

if __name__ == '__main__':
	unittest.main()

