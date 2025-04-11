#! /bin/env python3

from __future__ import division
import unittest

from solsticepy.gen_yaml import yamltransform

class TestGenYaml(unittest.TestCase):

	def test_1(self):
		s = yamltransform([1,2,3],[4,5,6])
		print(s)
		self.assertEqual(s,"xxx")

