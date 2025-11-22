#! /bin/env python3

from __future__ import division
#import unittest

from solsticepy.incident_angles import get_annual_incident
import numpy as np

def test_1():
	
	hst_pos=np.array([[100,200,200], [200,200,200]])  #np.array([[100,200,200]])
	hst_aims=np.array([[0,0,100], [0,0,100]]) #np.array([[0,0,100]])
	latitude=30

	inc_avg=get_annual_incident(hst_pos, hst_aims, latitude, casename='', verbose=False)
	print(inc_avg)
	assert(abs(inc_avg[0]-45.267)/inc_avg[0]<1e-4)
	assert(abs(inc_avg[1]-45.344)/inc_avg[1]<1e-4)
