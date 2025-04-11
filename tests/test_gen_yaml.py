#! /bin/env python3

from __future__ import division
import unittest

from solsticepy.gen_yaml import yamltransform

def test_1():
	s = yamltransform([1,2,3],[4,5,6])
	#print(s)
	assert(s=='transform: { translation: [1.000000e+00,2.000000e+00,3.000000e+00], rotation: [4.000000e+00,5.000000e+00,6.000000e+00] }')


