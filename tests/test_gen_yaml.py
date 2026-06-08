#! /bin/env python3

from __future__ import division
import unittest

import tempfile
from pathlib import Path

import numpy as np

from solsticepy.gen_yaml import Sun, gen_yaml, yamltransform

def test_1():
	s = yamltransform([1,2,3],[4,5,6])
	#print(s)
	assert(s=='transform: { translation: [1.000000e+00,2.000000e+00,3.000000e+00], rotation: [4.000000e+00,5.000000e+00,6.000000e+00] }')




def test_quadricxy_uses_azimuth_elevation_tracking():
	sun = Sun(dni=1000, sunshape='pillbox')
	hst_pos = np.array([[0., 100., 8.54], [20., 120., 8.54]])
	hst_aims = np.array([[0., 0., 62.], [0., 0., 62.]])
	hst_foc = np.array([[100., 110.], [120., 130.]])
	rec_param = np.r_[20., 20., 10, 10, 0., 0., 62., 0.]
	with tempfile.TemporaryDirectory() as tmpdir:
		tmp_path = Path(tmpdir)
		outfile_yaml = tmp_path / 'input.yaml'
		outfile_recv = tmp_path / 'input-rcv.yaml'

		gen_yaml(sun, hst_pos, hst_foc, hst_aims, 12.2, 12.2, 1., 1e-3, 'flat',
			rec_param, 1., outfile_yaml=outfile_yaml, outfile_recv=outfile_recv,
			shape='quadricxy')

		text = outfile_yaml.read_text()
		assert text.count('quadricxy:') == 2
		assert 'target_aligned' not in text
		assert 'ax2: 2.500000e-03' in text
		assert 'ay2: 1.923077e-03' in text
