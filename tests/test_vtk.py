#! /bin/env python3

from __future__ import division
import unittest
import solsticepy
from solsticepy.gen_vtk import read_vtk
import os
import numpy as np

class TestVTK(unittest.TestCase):

    def setUP(self):
        vtkfile='./data/example_rec.vtk'
        savedir='./data'
        dataname='Front_faces_Absorbed_flux'
        read_vtk(vtkfile, savedir, dataname)        

    def test_read_vtk(self):

        self.setUP()
        fn_data='./data/example_rec_mesh_data.csv'
        fn_points='./data/example_rec_points.csv'
        data=np.loadtxt(fn_data, delimiter=',', skiprows=1)
        pnts=np.loadtxt(fn_points, delimiter=',', skiprows=1)        
        
        self.assertTrue(os.path.exists(fn_data))
        self.assertTrue(os.path.exists(fn_points))
        self.assertTrue(np.shape(data)[0]==70)
        self.assertTrue(np.shape(data)[1]==4)        
        self.assertTrue(len(pnts)==48)        
        os.unlink(fn_data)
        os.unlink(fn_points)

if __name__ == '__main__':
	unittest.main()

