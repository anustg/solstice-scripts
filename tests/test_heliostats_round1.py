#! /bin/env python3

from __future__ import division
import unittest

from solsticepy.design_crs import CRS
import solsticepy
from solsticepy.master import Master
import os
import numpy as np
import time
from solsticepy.gen_vtk import read_vtk, flux_reader
import matplotlib.pyplot as plt

class TestHeliostats(unittest.TestCase):
    def setUp(self):
        '''
        Round Robin Tests: with collaboration with NREL and Tietronix

        Phase 1 Tests
        '''
        self.hemisphere='North'
        self.azimuth=270
        self.elevation=59.96
        self.DNI = 1000 # W/m2
        self.sunshape = 'pillbox'
        self.half_angle_deg = 0.2664   
        self.num_rays=int(20e6)#300e6,  1000e7

        # Heliostat
        self.rho_refl=0.9 # mirror reflectivity
        self.slope_error=2.e-3 # radians
        self.hst_w=11.415 # m
        self.hst_h=10.42 # m
        self.H_pedestal=5.65
        self.fct_w=1.605
        self.fct_h=2.06
        self.fct_row=5
        self.fct_col=7
        gap1=(self.hst_w-self.fct_col*self.fct_w)/(self.fct_col-1)
        gap2=(self.hst_h-self.fct_row*self.fct_h)/(self.fct_row-1)
        self.gap=(gap1+gap2)/2.
        print('gap = %.2f m'%self.gap)


        # Tower
        self.tower_h=0.01 # tower height
        self.tower_r=0.01 # tower radius
        
        # Receiver
        self.receiver='flat' # 'flat', 'cylinder' or 'stl'
        self.rec_w=12 # width, m
        self.rec_h=18 # height, m
        self.mesh_h=150 
        self.mesh_w= 100 

        self.loc_x=0. # m
        self.loc_y=0. # m
        self.loc_z=160 # m
        self.rec_abs=1.


    #@unittest.skip(" ")
    def test_task1(self):
        """ 
        Single-facet heliostat
        a) flat facet
        b) curved to slant range

        solstice_PhaseI_Task_1a_N
        solstice_PhaseI_Task_1a_SE
        solstice_PhaseI_Task_1b_N
        solstice_PhaseI_Task_1b_SE

        """
        cases=['a', 'b']
        facet_shapes=['flat', 'curved']
        locations=['N', 'SE']

        total_energy=np.array([[84.9447266666667, 74.5207333333333], [98.9551866666667, 80.0825733333333]])
        peak_flux=np.array([[0.851282, 0.821611], [2.670995, 4.2935 ]])
        

        for i,c in enumerate(cases):
            shape=facet_shapes[i]
            for j, loc in enumerate(locations):   
                start=time.time()    
                casename='solstice_PhaseI_Task_1%s_%s'%(c, loc)
                casedir='./'+casename
                if not os.path.exists(casedir):
                    os.makedirs(casedir)
                vtkfile=casedir+'/%s-%s-target_e.vtk'%(self.azimuth, self.elevation)
                 
                if not os.path.exists(vtkfile):
                    if loc=='N':
                        hst_pos=np.array([[0, 500, self.H_pedestal]])
                        tilt_x=0
                        tilt_y=0
                    elif loc=='SE':
                        hst_pos=np.array([[200, -200, self.H_pedestal]]) 
                        tilt_x=0
                        tilt_y=-225

                    hst_aims=np.array([[0, 0, self.loc_z]])
                    hst_foc=np.sqrt((hst_pos[0,0]-hst_aims[0,0])**2+(hst_pos[0,1]-hst_aims[0,1])**2+(hst_pos[0,2]-hst_aims[0,2])**2) 

                    rec_param=np.r_[self.rec_w, self.rec_h, self.mesh_w, self.mesh_h, self.loc_x, self.loc_y, self.loc_z, tilt_x, tilt_y]

                    master=Master(casedir=casedir)
                    outfile_yaml = master.in_case(folder=casedir, fn='input.yaml')
                    outfile_recv = master.in_case(folder=casedir, fn='input-rcv.yaml')

                    SUN = solsticepy.Sun(dni=self.DNI, sunshape=self.sunshape, half_angle_deg=self.half_angle_deg) 
                    solsticepy.gen_yaml(sun=SUN, 
                                        hst_pos=hst_pos, 
                                        hst_foc=hst_foc, 
                                        hst_aims=hst_aims, 
                                        hst_w=self.hst_w, 
                                        hst_h=self.hst_h,
                                        rho_refl=self.rho_refl, 
                                        slope_error=self.slope_error, 
                                        cant=False, 
                                        bands=np.array([None, None]), 
                                        receiver=self.receiver, 
                                        rec_param=rec_param, 
                                        rec_abs=self.rec_abs,
                                        outfile_yaml=outfile_yaml, 
                                        outfile_recv=outfile_recv,
                                        hemisphere=self.hemisphere, 
                                        tower_h=self.tower_h, 
                                        tower_r=self.tower_r,  
                                        spectral=False,
                                        medium=0, 
                                        one_heliostat=True, 
                                        fct_w=self.fct_w, 
                                        fct_h=self.fct_h, 
                                        fct_gap=self.gap, 
                                        n_row=self.fct_row, 
                                        n_col=self.fct_col, 
                                        shape=shape)
                   
                    master.run(self.azimuth, self.elevation, self.num_rays, self.rho_refl, self.DNI, folder=casedir, gen_vtk=True,  printresult=True, verbose=True, system='crs')
                end=time.time()
                print('Time: %.1f min'%((end-start)/60))
                points, tri, flux, flux_abs, flux_back, flux_abs_back=flux_reader(vtkfile, casedir)
                peak=convert_tri_to_rect(points, tri, flux, casedir)
                data=np.loadtxt(casedir+'/result-formatted.csv', dtype=str, delimiter=',')
                Qabs=data[-2,1].astype(float)
                print('Peak=%.2f, Qabs=%.2f'%(peak, Qabs))
                self.assertTrue(abs(Qabs-total_energy[i,j])/Qabs<0.02)
                self.assertTrue(abs(peak-peak_flux[i,j])/peak<0.02)
                

    #@unittest.skip(" ")
    def test_task3(self):
        """ 
        Multi-facet heliostat, flat facets, canted to slant range
 
        solstice_PhaseI_Task_3_N
        solstice_PhaseI_Task_3_SE


        """

        locations=['N', 'SE']

        total_energy=np.array([95.98694, 77.32753 ]) #96.11243, 77.9103166666667])
        peak_flux=np.array([2.512975, 3.873275])
        

        for j, loc in enumerate(locations):   
            start=time.time()    
            casename='solstice_PhaseI_Task_3_%s'%(loc)
            casedir='./'+casename
            if not os.path.exists(casedir):
                os.makedirs(casedir)
            vtkfile=casedir+'/%s-%s-target_e.vtk'%(self.azimuth, self.elevation)
             
            if not os.path.exists(vtkfile):
                if loc=='N':
                    hst_pos=np.array([[0, 500, self.H_pedestal]])
                    tilt_x=0
                    tilt_y=0
                elif loc=='SE':
                    hst_pos=np.array([[200, -200, self.H_pedestal]]) 
                    tilt_x=0
                    tilt_y=-225

                hst_aims=np.array([[0, 0, self.loc_z]])
                hst_foc=np.sqrt((hst_pos[0,0]-hst_aims[0,0])**2+(hst_pos[0,1]-hst_aims[0,1])**2+(hst_pos[0,2]-hst_aims[0,2])**2) 
                #hst_foc=2
                rec_param=np.r_[self.rec_w, self.rec_h, self.mesh_w, self.mesh_h, self.loc_x, self.loc_y, self.loc_z, tilt_x, tilt_y]

                master=Master(casedir=casedir)
                outfile_yaml = master.in_case(folder=casedir, fn='input.yaml')
                outfile_recv = master.in_case(folder=casedir, fn='input-rcv.yaml')

                SUN = solsticepy.Sun(dni=self.DNI, sunshape=self.sunshape, half_angle_deg=self.half_angle_deg) 
            
                solsticepy.gen_yaml(sun=SUN, 
                                    hst_pos=hst_pos, 
                                    hst_foc=hst_foc, 
                                    hst_aims=hst_aims, 
                                    hst_w=self.hst_w, 
                                    hst_h=self.hst_h,
                                    rho_refl=self.rho_refl, 
                                    slope_error=self.slope_error, 
                                    cant=True, 
                                    bands=np.array([None, None]), 
                                    receiver=self.receiver, 
                                    rec_param=rec_param, 
                                    rec_abs=self.rec_abs,
                                    outfile_yaml=outfile_yaml, 
                                    outfile_recv=outfile_recv,
                                    hemisphere=self.hemisphere, 
                                    tower_h=self.tower_h, 
                                    tower_r=self.tower_r,  
                                    spectral=False,
                                    medium=0, 
                                    one_heliostat=True, 
                                    fct_w=self.fct_w, 
                                    fct_h=self.fct_h, 
                                    fct_gap=self.gap, 
                                    n_row=self.fct_row, 
                                    n_col=self.fct_col, 
                                    shape='flat')
                         
                master.run(self.azimuth, self.elevation, self.num_rays, self.rho_refl, self.DNI, folder=casedir, gen_vtk=True,  printresult=True, verbose=True, system='crs')
            end=time.time()
            print('Time: %.1f min'%((end-start)/60))
            points, tri, flux, flux_abs, flux_back, flux_abs_back=flux_reader(vtkfile, casedir)
            peak=convert_tri_to_rect(points, tri, flux, casedir)
            data=np.loadtxt(casedir+'/result-formatted.csv', dtype=str, delimiter=',')
            Qabs=data[-2,1].astype(float)
            print('Peak=%.2f, Qabs=%.2f'%(peak, Qabs))

            self.assertTrue(abs(Qabs-total_energy[j])/Qabs<0.02)
            self.assertTrue(abs(peak-peak_flux[j])/peak<0.02)


    #@unittest.skip(" ")
    def test_task4(self):
        """ 
        Multi-facet heliostat, paraboloid facet focused to slant range, canted to slant range
 
        solstice_PhaseI_Task_4_N
        solstice_PhaseI_Task_4_SE


        """

        locations=['N', 'SE']

        total_energy=np.array([96.14, 77.33])
        peak_flux=np.array([2.593, 4.245 ])
        

        for j, loc in enumerate(locations):   
            start=time.time()    
            casename='solstice_PhaseI_Task_4_%s'%(loc)
            casedir='./'+casename
            if not os.path.exists(casedir):
                os.makedirs(casedir)
            vtkfile=casedir+'/%s-%s-target_e.vtk'%(self.azimuth, self.elevation)
             
            if not os.path.exists(vtkfile):
                if loc=='N':
                    hst_pos=np.array([[0, 500, self.H_pedestal]])
                    tilt_x=0
                    tilt_y=0
                elif loc=='SE':
                    hst_pos=np.array([[200, -200, self.H_pedestal]]) 
                    tilt_x=0
                    tilt_y=-225

                hst_aims=np.array([[0, 0, self.loc_z]])
                hst_foc=np.sqrt((hst_pos[0,0]-hst_aims[0,0])**2+(hst_pos[0,1]-hst_aims[0,1])**2+(hst_pos[0,2]-hst_aims[0,2])**2) 

                rec_param=np.r_[self.rec_w, self.rec_h, self.mesh_w, self.mesh_h, self.loc_x, self.loc_y, self.loc_z, tilt_x, tilt_y]

                master=Master(casedir=casedir)
                outfile_yaml = master.in_case(folder=casedir, fn='input.yaml')
                outfile_recv = master.in_case(folder=casedir, fn='input-rcv.yaml')

                SUN = solsticepy.Sun(dni=self.DNI, sunshape=self.sunshape, half_angle_deg=self.half_angle_deg) 
                solsticepy.gen_yaml(sun=SUN, 
                                    hst_pos=hst_pos, 
                                    hst_foc=hst_foc, 
                                    hst_aims=hst_aims, 
                                    hst_w=self.hst_w, 
                                    hst_h=self.hst_h,
                                    rho_refl=self.rho_refl, 
                                    slope_error=self.slope_error, 
                                    cant=True, 
                                    bands=np.array([None, None]), 
                                    receiver=self.receiver, 
                                    rec_param=rec_param, 
                                    rec_abs=self.rec_abs,
                                    outfile_yaml=outfile_yaml, 
                                    outfile_recv=outfile_recv,
                                    hemisphere=self.hemisphere, 
                                    tower_h=self.tower_h, 
                                    tower_r=self.tower_r,  
                                    spectral=False,
                                    medium=0, 
                                    one_heliostat=True, 
                                    fct_w=self.fct_w, 
                                    fct_h=self.fct_h, 
                                    fct_gap=self.gap, 
                                    n_row=self.fct_row, 
                                    n_col=self.fct_col, 
                                    shape='paraboloid')
               
                master.run(self.azimuth, self.elevation, self.num_rays, self.rho_refl, self.DNI, folder=casedir, gen_vtk=True,  printresult=True, verbose=True, system='crs')
            end=time.time()
            print('Time: %.1f min'%((end-start)/60))
            points, tri, flux, flux_abs, flux_back, flux_abs_back=flux_reader(vtkfile, casedir)
            peak=convert_tri_to_rect(points, tri, flux, casedir)
            data=np.loadtxt(casedir+'/result-formatted.csv', dtype=str, delimiter=',')
            Qabs=data[-2,1].astype(float)
            print('Peak=%.2f, Qabs=%.2f'%(peak, Qabs))
            self.assertTrue(abs(Qabs-total_energy[j])/Qabs<0.02)
            self.assertTrue(abs(peak-peak_flux[j])/peak<0.02)
                


def convert_tri_to_rect(points, tri, flux, casedir):
    '''
    Convert fluxmap in triangular mesh to rectangular mesh

    Arguments:
        flux, numpy array, flux data in triangular mesh from solstice
        points, points of the triangular mesh
        tri, indices of the triangular mesh

    Return:
        peak, float, peak flux
    '''
    print(np.shape(flux))	
    #test=np.arange(10)
    #print('even', test[::2])
    #print('odd', test[1::2])
    flux=(flux[::2]+flux[1::2])/2.

    m=150
    n=100

    width=12.
    height=18.
    flux=flux.reshape(n,m)
    flux=flux.T
    flux=np.flipud(flux)
    flux=np.fliplr(flux)
    print(np.shape(flux))
    #xx,yy=np.meshgrid(np.linspace(-9, 9, m+1),np.linspace(-6, 6, n+1)) 
    xx=np.linspace(-width/2., width/2., n+1)
    yy=np.linspace(-height/2., height/2., m+1)

    if ('b' in casedir) or 'Task_3' in casedir or 'Task_4' in casedir:
	    fmax=4.5 #2.8 #4.5
    else:
	    fmax=0.9 
    plt.pcolormesh(xx, yy, flux, vmax=fmax, vmin=0)

    plt.colorbar()
    plt.xlim([-width/2., width/2.])
    plt.ylim([-height/2.,height/2.])
    plt.yticks(np.linspace(-9,9,10))
    plt.gca().set_aspect('equal', adjustable='box')
    plt.savefig(open(casedir+'/flux_rect_%sx%s.png'%(m,n), 'wb'), bbox_inches='tight')
    #plt.show()
    plt.close()

    dx=float(width/n)
    dy=float(height/m)
    X=np.linspace(-width/2.+dx/2., width/2.-dx/2., n)
    Y=np.linspace(-height/2.+dy/2., height/2.-dy/2., m)
    XX,YY=np.meshgrid(X,Y)
    np.savetxt(casedir+'/%s_fluxmap.csv'%casedir, flux, fmt='%.6f', delimiter=',')
    np.savetxt(casedir+'/%s_xx.csv'%casedir, XX, fmt='%.2f', delimiter=',')
    np.savetxt(casedir+'/%s_yy.csv'%casedir, YY, fmt='%.2f', delimiter=',')

    peak=np.max(flux)
    return peak

if __name__ == '__main__':
	unittest.main()

