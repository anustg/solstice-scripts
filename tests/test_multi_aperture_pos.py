#! /bin/env python3

from __future__ import division
import unittest

import solsticepy
from solsticepy.cal_layout import multi_aperture_pos, radial_stagger
import os
import numpy as np
import time

import matplotlib.pyplot as plt

class TestMultiAperturePos(unittest.TestCase):

	def test_pos(self):

		omega, xc, yc=multi_aperture_pos(rec_w=10., gamma=180., n=3, i=2)
		self.assertEqual(omega, 180)

		omega, xc, yc=multi_aperture_pos(rec_w=10., gamma=240., n=3, i=2)
		self.assertEqual(omega, 210)

		# to plot some figures
		'''
		GAMMA=np.r_[60., 120., 180., 270., 360.]
		N=np.r_[2,3,4,5,6]
		
		for g in GAMMA[-1:]:
			for n in N[:1]:
				OMEGA=np.array([])
				for i in range(n):
					omega, xc, yc=multi_aperture_pos(rec_w=10., gamma=g, n=n, i=i)
					OMEGA=np.append(OMEGA, omega)

				ax = plt.subplot(111, projection='polar')
				ax.plot(OMEGA*np.pi/180., np.ones(n), 'o')
				plt.title('$\\gamma$=%s, n=%s'%(g, n))
				plt.show()
				plt.close()

	def test_layout_and_aiming(self):
		GAMMA=np.r_[60., 120., 180., 270., 360.]
		N=np.r_[2,3,4,5,6]
		
		for g in GAMMA:
			for n in N:
				rec_z=[]	
				for i in range(n):	
					rec_z.append(200.)	
				print('\n n:%s g:%s'%(n,g))
				pos_and_aiming, Nzones, Nrows_zone=radial_stagger(
				latitude=34.,
				num_hst=10000,
				width=12., 
				height=12., 
				hst_z=7., 
				towerheight=200., 
				R1=80., 
				fb=0.8, 
				dsep=0., 
				field='multi-aperture', 
				num_aperture=n, 
				gamma=g, 
				rec_w=20., 
				rec_z=rec_z, 
				savedir='.', 
				verbose=False, 
				plot=False, 
				plt_aiming='n%s_gamma%s'%(n,g))

		'''	

if __name__ == '__main__':
	unittest.main()

