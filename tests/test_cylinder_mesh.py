#! /bin/env python3

from __future__ import division
import unittest

import solsticepy
from solsticepy.design_crs import CRS
import solsticepy
from solsticepy.master import Master
import os
import numpy as np
import time
from solsticepy.gen_vtk import read_vtk, flux_reader

import re
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
from matplotlib import cm

class TestCylinder(unittest.TestCase):
    def setUp(self):
        self.hemisphere='South'

        self.azimuth=90
        self.elevation=62
        self.DNI = 1000 # W/m2
        self.sunshape = 'pillbox'
        self.half_angle_deg = 0.2664   
        self.num_rays=int(100e6)#1000e7

        # Heliostat
        self.rho_refl=0.9 # mirror reflectivity
        self.slope_error=2.e-3 # radians
        self.hst_w=10.38 # m
        self.hst_h=9.73 # m
        self.H_pedestal=4.49


        # Tower
        self.tower_h=0.01 # tower height
        self.tower_r=0.01 # tower radius
        
        # Receiver
        self.receiver='cylinder' # 'flat' or 'stl'
        self.rec_r=7.725 # width, m
        self.rec_h=25.05609 # height, m
        self.mesh_h=3 #31
        self.mesh_circ= 5 #60
        tilt=0.  # deg
        loc_x=0. # m
        loc_y=0. # m
        self.loc_z=180.33 #171.035 # m
        self.rec_abs=1.
        self.rec_param=np.r_[self.rec_r*2., self.rec_h, self.mesh_circ, self.mesh_h, loc_x, loc_y, self.loc_z, tilt]


    def test(self):
        """ 
        Heliostat ID 8993
        Single facet, focused to slant range, single heliostat, no blocking and shading
        """

        azi=self.azimuth
        ele=self.elevation

        casedir='./results/test-cylinder-mesh'
        if not os.path.exists(casedir):
            os.makedirs(casedir)

        shape='curved'#, 'curved' #'flat'
        cant=False
        case_id=8993
        hst_x=np.r_[582.49]
        hst_y=np.r_[-490.35]

        hst_z=np.ones(len(hst_x))*self.H_pedestal
        hst_pos=np.append(hst_x, (hst_y, hst_z))
        hst_pos=hst_pos.reshape(3, len(hst_x))
        hst_pos=hst_pos.T

        # aim at the receiver center
        hst_azimuth=np.arccos(hst_y/np.sqrt(hst_x**2+hst_y**2))
        hst_azimuth[hst_x>0]=np.pi*2.-hst_azimuth[hst_x>0]
        hst_aims=np.zeros((len(hst_x),3))
        hst_aims[:,0]=-self.rec_r*np.sin(hst_azimuth)
        hst_aims[:,1]=self.rec_r*np.cos(hst_azimuth)
        hst_aims[:,2]=self.loc_z

        # slant range focus
        hst_foc=np.sqrt((hst_x-hst_aims[:,0])**2+(hst_y-hst_aims[:,1])**2+(hst_z-hst_aims[:,2])**2)
        bands=np.array([[None, None]]) 

        master=Master(casedir=casedir)
        outfile_yaml = master.in_case(folder=casedir, fn='input.yaml')
        outfile_recv = master.in_case(folder=casedir, fn='input-rcv.yaml')

        SUN = solsticepy.Sun(dni=self.DNI, sunshape=self.sunshape, half_angle_deg=1e-4)#self.half_angle_deg) 
      
        solsticepy.gen_yaml(sun=SUN, 
                            hst_pos=hst_pos, 
                            hst_foc=hst_foc, 
                            hst_aims=hst_aims, 
                            hst_w=self.hst_w, 
                            hst_h=self.hst_h,
                            rho_refl=self.rho_refl, 
                            slope_error=0., 
                            cant=cant, 
                            bands=bands, 
                            receiver=self.receiver, 
                            rec_param=self.rec_param, 
                            rec_abs=self.rec_abs,
                            outfile_yaml=outfile_yaml, 
                            outfile_recv=outfile_recv,
                            hemisphere=self.hemisphere, 
                            tower_h=self.tower_h, 
                            tower_r=self.tower_r,  
                            spectral=False,
                            medium=0, 
                            one_heliostat=True, 
                            shape=shape)
        
        #master.run(azi, ele, int(self.num_rays/10), self.rho_refl, self.DNI, folder=casedir, gen_vtk=True,  printresult=True, verbose=True, system='crs')
        casename='test-cylinder-mesh'
        vtkfile=casedir+'/%s-%s-target_e.vtk'%(azi, ele)
        points, tri, flux, flux_abs, flux_back, flux_abs_back=flux_reader(vtkfile, casedir, check=True)
	
        num_points=len(points)
        num_cells=len(tri)

        self.assertEqual(num_points, (self.mesh_h+1)*self.mesh_circ+2)	
        self.assertEqual(num_cells,  self.mesh_h*self.mesh_circ*2+self.mesh_circ*2)	
        #
        # Note
        # the mesh of a cylinder in Solstice includes the top and bottom surfaces
        # the total number of cells on the cylindrical side part is self.mesh_h*self.mesh_circ*2  
        # the number of cells on the bottom and top surfaces is self.mesh_circ*2, (bottom first and top later)


if __name__ == '__main__':
	unittest.main()

