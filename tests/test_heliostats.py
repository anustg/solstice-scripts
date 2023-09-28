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
import matplotlib.pyplot as plt

class TestHeliostats(unittest.TestCase):
    def setUp(self):
        self.hemisphere='South'
        #dd=22
        #mm='Sep'
        #lat=-28.29945   # latitude of the crs plant
        #log=23.364761   # South Africa
        #sun=solsticepy.cal_sun.SunPosition()
        #day=sun.days(dd, mm)
        #delta=sun.declination(day)
        #omega=0. #-15.*4 (8am)  # solar noon is 0
        #theta=sun.zenith(lat, delta, omega)
        #azi=sun.azimuth(lat, theta, delta, omega)
        #azimuth, elevation=sun.convert_convention('solstice', azi, theta)
        self.azimuth=90
        self.elevation=62
        self.DNI = 1000 # W/m2
        self.sunshape = 'pillbox'
        self.half_angle_deg = 0.2664   
        self.num_rays=int(20e6)#300e6,  1000e7

        # Heliostat
        self.rho_refl=0.9 # mirror reflectivity
        self.slope_error=2.e-3 #1e-9 #2.e-3 # radians
        self.hst_w=10.38 # m
        self.hst_h=9.73 # m
        self.H_pedestal=4.49
        self.fct_w=2.06
        self.fct_h=1.605
        self.fct_row=6
        self.fct_col=5
        gap1=(self.hst_w-self.fct_col*self.fct_w)/(self.fct_col-1)
        gap2=(self.hst_h-self.fct_row*self.fct_h)/(self.fct_row-1)
        self.gap=(gap1+gap2)/2.

        # Tower
        self.tower_h=0.01 # tower height
        self.tower_r=0.01 # tower radius
        
        # Receiver
        self.receiver='cylinder' # 'flat' or 'stl'
        self.rec_r=7.725 # width, m
        self.rec_h=25.05609 # height, m
        self.mesh_h=31 #31
        self.mesh_circ= 60 #60
        tilt=0.  # deg
        loc_x=0. # m
        loc_y=0. # m
        self.loc_z=180.33 #171.035 # m
        self.rec_abs=1.
        self.rec_param=np.r_[self.rec_r*2., self.rec_h, self.mesh_circ, self.mesh_h, loc_x, loc_y, self.loc_z, tilt]


    @unittest.skip(" ")
    def test_1(self):
        """ 
        Whole field
        single curved facet heliostats
        a) perfect cant, focus
        b) canting bands
        solar noon
        morning 8 am
        central aiming
        designed aiming

        solstice_Task_1a_AimStrat_1_12_fluxmap
        solstice_Task_1a_AimStrat_1_8_fluxmap
        solstice_Task_1a_AimStrat_2_12_fluxmap
        solstice_Task_1a_AimStrat_2_8_fluxmap
        solstice_Task_1b_AimStrat_1_12_fluxmap
        solstice_Task_1b_AimStrat_1_8_fluxmap
        solstice_Task_1b_AimStrat_2_12_fluxmap
        solstice_Task_1b_AimStrat_2_8_fluxmap

        """
        
        shape='curved'#, 'curved' #'flat'
        cant=False
        
        time=[12, 8]
        focuses=['a', 'b']
        aims=[1, 2]
        for t in time:
            for f in focuses:
                for a in aims:
                    casename='solstice_Task_1%s_AimStrat_%s_%s'%(f, a, t)
                    casedir='./'+casename
                    if not os.path.exists(casedir):
                        os.makedirs(casedir)

                    if t==12:
                        azi=self.azimuth
                        ele=self.elevation
                        vtkfile=casedir+'/%.0f-%.0f-target_e.vtk'%(azi, ele)
                    elif t==8:
                        azi=14.7333
                        ele=26.4378
                        vtkfile=casedir+'/%s-%s-target_e.vtk'%(azi, ele)

                    if not os.path.exists(casedir+'/flux_tri.png'):

                        if '1a' in casedir:
                            bands=np.array([[None, None]])
                        else:
                            bands=np.array([[502, 516],  # band range (<=), focal length
                                [885, 668],
                                [1267, 959],
                                [1650, 1500]])

                        layout=np.loadtxt('./data/heliostats_pos_ID.csv', delimiter=',', skiprows=1)
                        if a==1:
                            offset_z=0
                        elif a==2:
                            offset_z=layout[:,4]

                        hst_x=layout[:,1]
                        hst_y=layout[:,2]
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
                        hst_aims[:,2]=self.loc_z+offset_z

                        # slant range focus
                        #hst_foc=np.sqrt((hst_x-hst_aims[:,0])**2+(hst_y-hst_aims[:,1])**2+(hst_z-hst_aims[:,2])**2)
                        #hst_foc=np.sqrt((hst_x)**2+(hst_y)**2+(hst_z)**2) # slant range is the heliostat to the bottom of the tower (0,0,0)
                        hst_foc=np.sqrt((hst_x)**2+(hst_y)**2+(hst_z-self.loc_z)**2) #slant range is the centre point of the receiver cylinder

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
                                            one_heliostat=False, 
                                            fct_w=self.fct_w, 
                                            fct_h=self.fct_h, 
                                            fct_gap=self.gap, 
                                            n_row=self.fct_row, 
                                            n_col=self.fct_col, 
                                            shape=shape)
                       
                        master.run(azi, ele, self.num_rays, self.rho_refl, self.DNI, folder=casedir, gen_vtk=True,  printresult=True, verbose=True, system='crs')
                        
                    points, tri, flux, flux_abs, flux_back, flux_abs_back=flux_reader(vtkfile, casedir)
                    plot_fluxmap(points, tri, flux, casedir, casename=casename, loc_z_rec=self.loc_z, rec_r=self.rec_r, rec_h=self.rec_h, m=self.mesh_h, n=self.mesh_circ)

    #@unittest.skip(" ")
    def test_2(self):
        """ 
        Whole field
        canted, multi facets heliostats, flat facets
        a) perfect cant, focus
        b) canting bands
        solar noon
        morning 8 am
        central aiming
        designed aiming

        solstice_Task_2a_AimStrat_1_12_fluxmap
        solstice_Task_2a_AimStrat_1_8_fluxmap
        solstice_Task_2a_AimStrat_2_12_fluxmap
        solstice_Task_2a_AimStrat_2_8_fluxmap
        solstice_Task_2b_AimStrat_1_12_fluxmap
        solstice_Task_2b_AimStrat_1_8_fluxmap
        solstice_Task_2b_AimStrat_2_12_fluxmap
        solstice_Task_2b_AimStrat_2_8_fluxmap

        """
        
        shape='flat'#, 'curved' #'flat'
        cant=True
        
        time=[12, 8]
        focuses=['a', 'b']
        aims=[1]#, 2]
        for a in aims:
            for t in time:
                for f in focuses:


                    layout=np.loadtxt('./data/heliostats_pos_ID.csv', delimiter=',', skiprows=1)
                    if a==1:
                        offset_z=0
                    elif a==2:
                        offset_z=layout[:,4]

                    num=len(layout)
                    m=500
                    casename_0='solstice_Task_2%s_AimStrat_%s_%s'%(f, a, t)
                    casedir_0='./'+casename_0
                    for i in range(int(num/m)+1):
                        casename='solstice_Task_2%s_AimStrat_%s_%s-%s'%(f, a, t, i)
                        casedir=casedir_0+'/%s'%casename
                        if not os.path.exists(casedir):
                            os.makedirs(casedir)

                        if not os.path.exists(casedir+'/flux_tri.png'):

                            if 'a' in casedir:
                                bands=np.array([[None, None]])
                            else:
                                bands=np.array([[502, 516],  # band range (<=), focal length
                                    [885, 668],
                                    [1267, 959],
                                    [1650, 1500]])

                            if t==12:
                                azi=self.azimuth
                                ele=self.elevation
                                vtkfile=casedir+'/%.0f-%.0f-target_e.vtk'%(azi, ele)
                            elif t==8:
                                azi=14.7333
                                ele=26.4378
                                vtkfile=casedir+'/%s-%s-target_e.vtk'%(azi, ele)

                            if i<=int(num/m)-1:                                         
                                hst_x=layout[i*m:(i+1)*m,1]
                                hst_y=layout[i*m:(i+1)*m,2]
                                if a==2:
                                    offset_z=offset_z[i*m:(i+1)*m]
                            else:
                                hst_x=layout[i*m:,1]
                                hst_y=layout[i*m:,2]
                                if a==2:
                                    offset_z=offset_z[i*m:]
                            print(len(hst_x), i)
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
                            hst_aims[:,2]=self.loc_z+offset_z

                            # slant range focus
                            #hst_foc=np.sqrt((hst_x-hst_aims[:,0])**2+(hst_y-hst_aims[:,1])**2+(hst_z-hst_aims[:,2])**2)
                            hst_foc=np.sqrt((hst_x)**2+(hst_y)**2+(hst_z-self.loc_z)**2) #slant range is the centre point of the receiver cylinder

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
                                                one_heliostat=False, 
                                                fct_w=self.fct_w, 
                                                fct_h=self.fct_h, 
                                                fct_gap=self.gap, 
                                                n_row=self.fct_row, 
                                                n_col=self.fct_col, 
                                                shape=shape)
                           
                            master.run(azi, ele, self.num_rays, self.rho_refl, self.DNI, folder=casedir, gen_vtk=True,  printresult=True, verbose=True, system='crs')
                            points, tri, flux, flux_abs, flux_back, flux_abs_back=flux_reader(vtkfile, casedir)
                            plot_fluxmap(points, tri, flux, casedir, casename=casename, loc_z_rec=self.loc_z, rec_r=self.rec_r, rec_h=self.rec_h, m=self.mesh_h, n=self.mesh_circ)
                            

                    width=2.*np.pi*self.rec_r
                    height=self.rec_h
                    FLUX=np.zeros((31,60))
                    RESULTS=np.zeros(10)
                    for i in range(int(num/m)+1):
                        casename='solstice_Task_2%s_AimStrat_%s_%s-%s'%(f, a, t, i)
                        casedir='/'+casename

                        flux=np.loadtxt(casedir_0+'%s/%s_fluxmap.csv'%(casedir, casename), delimiter=',')
                        FLUX+=flux
                        XX=np.loadtxt(casedir_0+'%s/%s_xx.csv'%(casedir, casename), delimiter=',')
                        YY=np.loadtxt(casedir_0+'%s/%s_yy.csv'%(casedir, casename), delimiter=',')
                    
                        data=np.loadtxt(casedir_0+'%s/result-formatted.csv'%(casedir), delimiter=',', dtype=str)
                        res=data[1:,1].astype(float)
                        RESULTS+=res

                    data[1:,1]=RESULTS          
                    np.savetxt(casedir_0+'/solstice_Task_2%s_AimStrat_%s_%s_results.csv'%(f, a, t), data[:,:2], fmt='%s', delimiter=',')	  

                    np.savetxt(casedir_0+'/solstice_Task_2%s_AimStrat_%s_%s_fluxmap.csv'%(f, a, t), FLUX, fmt='%.6f', delimiter=',')
                    np.savetxt(casedir_0+'/solstice_Task_2%s_AimStrat_%s_%s_xx.csv'%(f, a, t), XX, fmt='%.2f', delimiter=',')
                    np.savetxt(casedir_0+'/solstice_Task_2%s_AimStrat_%s_%s_yy.csv'%(f, a, t), YY, fmt='%.2f', delimiter=',')	                    
                    
                    plt.pcolormesh(XX[0], YY[:,0], FLUX, cmap='jet')#, vmax=2400, vmin=0)
                    plt.colorbar()
                    plt.xlim([-width/2., width/2.])
                    plt.ylim([-height/2.,height/2.])
                    plt.gca().set_aspect('equal', adjustable='box')
                    plt.savefig(open(casedir_0+'/solstice_Task_2%s_AimStrat_%s_%s_flux_map.png'%(f, a, t), 'wb'), bbox_inches='tight')
                    #plt.show()
                    plt.close()
                                    
              

    @unittest.skip(" ")
    def test_3(self):
        """ 
        Whole field
        canted, multi facets heliostats, curved facets
        a) perfect cant, focus
        b) canting bands
        solar noon
        morning 8 am
        central aiming
        designed aiming

        solstice_Task_3a_AimStrat_1_12_fluxmap
        solstice_Task_3a_AimStrat_1_8_fluxmap
        solstice_Task_3a_AimStrat_2_12_fluxmap
        solstice_Task_3a_AimStrat_2_8_fluxmap
        solstice_Task_3b_AimStrat_1_12_fluxmap
        solstice_Task_3b_AimStrat_1_8_fluxmap
        solstice_Task_3b_AimStrat_2_12_fluxmap
        solstice_Task_3b_AimStrat_2_8_fluxmap

        """
        
        shape='curved'#, 'curved' #'flat'
        cant=True
        
        time=[12, 8]
        focuses=['a', 'b']
        aims=[1, 2]
        for t in time:
            for f in focuses:
                for a in aims:

                    layout=np.loadtxt('./data/heliostats_pos_ID.csv', delimiter=',', skiprows=1)
                    if a==1:
                        offset_z=0
                    elif a==2:
                        offset_z=layout[:,4]

                    num=len(layout)
                    m=500

                    for i in range(int(num/m)):
                        casename='solstice_Task_2%s_AimStrat_%s_%s-%s'%(f, a, t, i)
                        casedir='./test-heliostats-'+casename
                        if not os.path.exists(casedir):
                            os.makedirs(casedir)

                        if not os.path.exists(casedir+'/flux_tri.png'):

                            if 'a' in casedir:
                                bands=np.array([[None, None]])
                            else:
                                bands=np.array([[502, 516],  # band range (<=), focal length
                                    [885, 668],
                                    [1267, 959],
                                    [1650, 1500]])
                                #TODO facets are curved in different ranges and focal lengthes  
                                  

                            if t==12:
                                azi=self.azimuth
                                ele=self.elevation
                                vtkfile=casedir+'/%.0f-%.0f-target_e.vtk'%(azi, ele)
                            elif t==8:
                                azi=14.7333
                                ele=26.4378
                                vtkfile=casedir+'/%s-%s-target_e.vtk'%(azi, ele)

                            if i<int(num/m)-1:                                         
                                hst_x=layout[i*m:(i+1)*m,1]
                                hst_y=layout[i*m:(i+1)*m,2]
                            else:
                                hst_x=layout[i*m:,1]
                                hst_y=layout[i*m:,2]
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
                            hst_aims[:,2]=self.loc_z+offset_z

                            # slant range focus
                            hst_foc=np.sqrt((hst_x-hst_aims[:,0])**2+(hst_y-hst_aims[:,1])**2+(hst_z-hst_aims[:,2])**2)


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
                                                one_heliostat=False, 
                                                fct_w=self.fct_w, 
                                                fct_h=self.fct_h, 
                                                fct_gap=self.gap, 
                                                n_row=self.fct_row, 
                                                n_col=self.fct_col, 
                                                shape=shape)
                           
                            master.run(azi, ele, self.num_rays, self.rho_refl, self.DNI, folder=casedir, gen_vtk=True,  printresult=True, verbose=True, system='crs')
                            points, tri, flux, flux_abs, flux_back, flux_abs_back=flux_reader(vtkfile, casedir)
                            plot_fluxmap(points, tri, flux, casedir, casename=casename, loc_z_rec=self.loc_z, rec_r=self.rec_r, rec_h=self.rec_h, m=self.mesh_h, n=self.mesh_circ)

    @unittest.skip(" ")
    def test_4c(self):
        """ 
        Heliostat ID 5473
        Multi-facets, canting and focusing bands, single heliostat with surrounding helios for blocking and shading
        """

        casedir='./test-heliostats-Task_4c'
        if not os.path.exists(casedir):
            os.makedirs(casedir)

        shape='curved'#, 'curved' #'flat'
        cant=True
        case_id=5473

        hst_x, hst_y, hst_z, hst_pos=heliostat_selections(case_id, casedir,self.H_pedestal)

        # aim at the receiver center
        hst_azimuth=np.arccos(hst_y/np.sqrt(hst_x**2+hst_y**2))
        hst_azimuth[hst_x>0]=np.pi*2.-hst_azimuth[hst_x>0]
        hst_aims=np.zeros((len(hst_x),3))
        hst_aims[:,0]=-self.rec_r*np.sin(hst_azimuth)
        hst_aims[:,1]=self.rec_r*np.cos(hst_azimuth)
        hst_aims[:,2]=self.loc_z

        # slant range focus
        hst_foc=np.sqrt((hst_x-hst_aims[:,0])**2+(hst_y-hst_aims[:,1])**2+(hst_z-hst_aims[:,2])**2)
        bands=np.array([[502, 516],  # band range (<=), focal length
                [885, 668],
                [1267, 959],
                [1650, 1500]])

        master=Master(casedir=casedir)
        outfile_yaml = master.in_case(folder=casedir, fn='input.yaml')
        outfile_recv = master.in_case(folder=casedir, fn='input-rcv.yaml')

        SUN = solsticepy.Sun(dni=self.DNI, sunshape=self.sunshape, half_angle_deg=self.half_angle_deg) 
        '''
        solsticepy.gen_yaml(sun=SUN, 
                            hst_pos=hst_pos, 
                            hst_foc=hst_foc, 
                            hst_aims=hst_aims, 
                            hst_w=self.hst_w, 
                            hst_h=self.hst_h,
	                        rho_refl=self.rho_refl, 
                            slope_error=self.slope_error, 
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
                            one_heliostat=False, 
                            fct_w=self.fct_w, 
                            fct_h=self.fct_h, 
                            fct_gap=self.gap, 
                            n_row=self.fct_row, 
                            n_col=self.fct_col, 
                            shape=shape)
        '''
        master.run(self.azimuth, self.elevation, self.num_rays*50, self.rho_refl, self.DNI, folder=casedir, gen_vtk=True,  printresult=True, verbose=True, system='crs')
        vtkfile=casedir+'/%.0f-%.0f-target_e.vtk'%(self.azimuth, self.elevation)
        points, tri, flux, flux_abs, flux_back, flux_abs_back=flux_reader(vtkfile, casedir)
        plot_fluxmap(points, tri, flux, casedir, loc_z_rec=self.loc_z, rec_r=self.rec_r, rec_h=self.rec_h, m=self.mesh_h, n=self.mesh_circ)


    @unittest.skip(" ")
    def test_4d(self):
        """ 
        Heliostat ID 8993
        Single facet, focused to slant range, single heliostat, no blocking and shading
        """
        time=[12, 8]
        for t in time:
            if t==12:
                azi=self.azimuth
                ele=self.elevation
            elif t==8:
                azi=14.7333
                ele=26.4378

            casedir='./solstice_Task_4d_AimStrat_1_%.0f-no slope-0 sunshape'%t
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
                                fct_w=self.fct_w, 
                                fct_h=self.fct_h, 
                                fct_gap=self.gap, 
                                n_row=self.fct_row, 
                                n_col=self.fct_col, 
                                shape=shape)
            
            master.run(azi, ele, int(self.num_rays/10), self.rho_refl, self.DNI, folder=casedir, gen_vtk=True,  printresult=True, verbose=True, system='crs')
            casename='solstice_Task_4d_AimStrat_1_%.0f-no slope-0 sunshape'%t
            vtkfile=casedir+'/%s-%s-target_e.vtk'%(azi, ele)
            points, tri, flux, flux_abs, flux_back, flux_abs_back=flux_reader(vtkfile, casedir)
            print('\n***********')
            plot_fluxmap(points, tri, flux, casedir, casename, loc_z_rec=self.loc_z, rec_r=self.rec_r, rec_h=self.rec_h, m=self.mesh_h, n=self.mesh_circ)
            print('***********\n\n')

    @unittest.skip(" ")
    def test_4efgh(self):
        """ 
        Heliostat ID 8993
        Single facet, flat mirror, single heliostat, no blocking and shading, no slope error, no sunshape
        """
        time=[12, 8]
        labels=['e', 'f', 'g', 'h', 'i']
        shapes=['curved', 'curved', 'curved', 'curved', 'flat']
        
        aim_offsets_AZI=np.r_[1./self.rec_r, -1./self.rec_r, 0., 0., 0.]
        aim_offsets_Z=np.r_[0., 0., 1., -1., 0.]

        for i in range(5):
            label=labels[i]
            aimoff_azi=aim_offsets_AZI[i]
            aimoff_Z=aim_offsets_Z[i]

            if label=='i':
                fct_w=self.hst_w
                fct_h=self.hst_h
            else:
                fct_w=self.fct_w
                fct_h=self.fct_h                

            for t in time:
                if t==12:
                    azi=self.azimuth
                    ele=self.elevation
                elif t==8:
                    azi=14.7333
                    ele=26.4378

                casedir='./results/solstice_Task_4%s_time%.0f'%(label, t)
                if not os.path.exists(casedir):
                    os.makedirs(casedir)

                shape=shapes[i] #'curved'#, 'curved' #'flat'
                cant=False
                case_id=8993
                hst_x=np.r_[582.489]
                hst_y=np.r_[-490.346]

                hst_z=np.ones(len(hst_x))*self.H_pedestal
                hst_pos=np.append(hst_x, (hst_y, hst_z))
                hst_pos=hst_pos.reshape(3, len(hst_x))
                hst_pos=hst_pos.T

                # aim at the receiver center
                hst_azimuth=np.arccos(hst_y/np.sqrt(hst_x**2+hst_y**2))
                hst_azimuth[hst_x>0]=np.pi*2.-hst_azimuth[hst_x>0]
                hst_aims=np.zeros((len(hst_x),3))
                hst_aims[:,0]=-self.rec_r*np.sin(hst_azimuth+aimoff_azi)
                hst_aims[:,1]=self.rec_r*np.cos(hst_azimuth+aimoff_azi)
                hst_aims[:,2]=self.loc_z+aimoff_Z
                #print(hst_aims[:,0], hst_aims[:,1],hst_aims[:,2])    


                sx=-self.rec_r*np.sin(hst_azimuth)
                sy=self.rec_r*np.cos(hst_azimuth)
                sz=self.loc_z
                # slant range focus
                hst_foc=np.sqrt((hst_x-sx)**2+(hst_y-sy)**2+(hst_z-sz)**2)   #np.sqrt((hst_x-hst_aims[:,0])**2+(hst_y-hst_aims[:,1])**2+(hst_z-hst_aims[:,2])**2) # None #
                bands=np.array([[None, None]]) 
                #print(hst_foc)


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
                                    fct_w=fct_w, 
                                    fct_h=fct_h, 
                                    fct_gap=self.gap, 
                                    n_row=self.fct_row, 
                                    n_col=self.fct_col, 
                                    shape=shape)
                
                master.run(azi, ele, int(self.num_rays/10), self.rho_refl, self.DNI, folder=casedir, gen_vtk=True,  printresult=True, verbose=True, system='crs')
                casename='solstice_Task_4%s_AimStrat_1_%.0f'%(label, t)
                vtkfile=casedir+'/%s-%s-target_e.vtk'%(azi, ele)
                points, tri, flux, flux_abs, flux_back, flux_abs_back=flux_reader(vtkfile, casedir)
                print('\n***********')
                plot_fluxmap(points, tri, flux, casedir, casename, loc_z_rec=self.loc_z, rec_r=self.rec_r, rec_h=self.rec_h, m=self.mesh_h, n=self.mesh_circ)
                print('***********\n\n')



    @unittest.skip(" ")
    def test_7a(self):
        """ 
        Heliostat ID 5473
        Multi-facets, slant range canting, slant range focusing, single heliostat, 0 slope error
        """

        casedir='./test-heliostats-Task_7a-1'
        if not os.path.exists(casedir):
            os.makedirs(casedir)

        hst_x=np.r_[46.51]
        hst_y=np.r_[-1580.01]
        shape='curved'#, 'curved' #'flat'
        cant=True

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

        SUN = solsticepy.Sun(dni=self.DNI, sunshape=self.sunshape, half_angle_deg=self.half_angle_deg) 
        solsticepy.gen_yaml(sun=SUN, 
                            hst_pos=hst_pos, 
                            hst_foc=hst_foc, 
                            hst_aims=hst_aims, 
                            hst_w=self.hst_w, 
                            hst_h=self.hst_h,
	                        rho_refl=self.rho_refl, 
                            slope_error=1e-9, 
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
                            one_heliostat=False, 
                            fct_w=self.fct_w, 
                            fct_h=self.fct_h, 
                            fct_gap=self.gap, 
                            n_row=self.fct_row, 
                            n_col=self.fct_col, 
                            shape=shape)

        master.run(self.azimuth, self.elevation, self.num_rays, self.rho_refl, self.DNI, folder=casedir, gen_vtk=True,  printresult=True, verbose=True, system='crs')
        vtkfile=casedir+'/%.0f-%.0f-target_e.vtk'%(self.azimuth, self.elevation)
        points, tri, flux, flux_abs, flux_back, flux_abs_back=flux_reader(vtkfile, casedir)
        plot_fluxmap(points, tri, flux, casedir, loc_z_rec=self.loc_z, rec_r=self.rec_r, rec_h=self.rec_h, m=self.mesh_h, n=self.mesh_circ)

    @unittest.skip(" ")
    def test_7c(self):
        """ 
        Heliostat ID 5473, noon
        Single facet, perfect focusing, 0 slope error, with surrounding heliostats for blcoking and shading
        """

        casedir='./test-heliostats-Task_7c'
        if not os.path.exists(casedir):
            os.makedirs(casedir)

        shape='curved'#, 'curved' #'flat'
        cant=False
        case_id=5473

        hst_x, hst_y, hst_z, hst_pos=heliostat_selection(case_id, casedir,self.H_pedestal)

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

        SUN = solsticepy.Sun(dni=self.DNI, sunshape=self.sunshape, half_angle_deg=self.half_angle_deg) 
        '''
        solsticepy.gen_yaml(sun=SUN, 
                            hst_pos=hst_pos, 
                            hst_foc=hst_foc, 
                            hst_aims=hst_aims, 
                            hst_w=self.hst_w, 
                            hst_h=self.hst_h,
	                        rho_refl=self.rho_refl, 
                            slope_error=1e-9, 
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
                            one_heliostat=False, 
                            fct_w=self.fct_w, 
                            fct_h=self.fct_h, 
                            fct_gap=self.gap, 
                            n_row=self.fct_row, 
                            n_col=self.fct_col, 
                            shape=shape)
        '''

        master.run(self.azimuth, self.elevation, self.num_rays*50, self.rho_refl, self.DNI, folder=casedir, gen_vtk=True,  printresult=True, verbose=True, system='crs')
        vtkfile=casedir+'/%.0f-%.0f-target_e.vtk'%(self.azimuth, self.elevation)
        points, tri, flux, flux_abs, flux_back, flux_abs_back=flux_reader(vtkfile, casedir)
        plot_fluxmap(points, tri, flux, casedir, loc_z_rec=self.loc_z, rec_r=self.rec_r, rec_h=self.rec_h, m=self.mesh_h, n=self.mesh_circ)

    @unittest.skip(" ")
    def test_7d(self):
        """ 
        Heliostat ID 5473, 8am
        single facet, slant range canting, slant range focusing, single heliostat, 0 slope error
        """

        casedir='./test-heliostats-Task_7d'
        if not os.path.exists(casedir):
            os.makedirs(casedir)

        hst_x=np.r_[46.51]
        hst_y=np.r_[-1580.01]
        shape='curved'#, 'curved' #'flat'
        cant=False

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

        SUN = solsticepy.Sun(dni=self.DNI, sunshape=self.sunshape, half_angle_deg=self.half_angle_deg) 
        solsticepy.gen_yaml(sun=SUN, 
                            hst_pos=hst_pos, 
                            hst_foc=hst_foc, 
                            hst_aims=hst_aims, 
                            hst_w=self.hst_w, 
                            hst_h=self.hst_h,
	                        rho_refl=self.rho_refl, 
                            slope_error=1e-9, 
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
                            fct_w=self.fct_w, 
                            fct_h=self.fct_h, 
                            fct_gap=self.gap, 
                            n_row=self.fct_row, 
                            n_col=self.fct_col, 
                            shape=shape)

        azi=14.7333
        ele=26.4378
        casename='solstice_Task_7d'
        #master.run(azi, ele, self.num_rays, self.rho_refl, self.DNI, folder=casedir, gen_vtk=True,  printresult=True, verbose=True, system='crs')
        vtkfile=casedir+'/%.4f-%.4f-target_e.vtk'%(azi, ele)
        points, tri, flux, flux_abs, flux_back, flux_abs_back=flux_reader(vtkfile, casedir)
        plot_fluxmap(points, tri, flux, casedir, casename, loc_z_rec=self.loc_z, rec_r=self.rec_r, rec_h=self.rec_h, m=self.mesh_h, n=self.mesh_circ)

    @unittest.skip(" ")
    def test_8a(self):
        """ 
        Each single facet of heliostat ID 5473
        slant range canting, slant range focusing, 0 slope error
        With the sun as a point source
        """

        casedir='./test-heliostats-Task_8'
        if not os.path.exists(casedir):
            os.makedirs(casedir)

        hst_x=np.r_[46.51]
        hst_y=np.r_[-1580.01]
        shape='curved'#, 'curved' #'flat'
        cant=True

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

        SUN = solsticepy.Sun(dni=self.DNI, sunshape=self.sunshape, half_angle_deg=1e-4) 
		
        solsticepy.gen_yaml(sun=SUN, 
                            hst_pos=hst_pos, 
                            hst_foc=hst_foc, 
                            hst_aims=hst_aims, 
                            hst_w=self.hst_w, 
                            hst_h=self.hst_h,
	                        rho_refl=self.rho_refl, 
                            slope_error=1e-9, 
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
                            one_heliostat=False, 
                            fct_w=self.fct_w, 
                            fct_h=self.fct_h, 
                            fct_gap=self.gap, 
                            n_row=self.fct_row, 
                            n_col=self.fct_col, 
                            shape=shape)

        master.run(self.azimuth, self.elevation, self.num_rays, self.rho_refl, self.DNI, folder=casedir, gen_vtk=True,  printresult=True, verbose=True, system='crs')
        vtkfile=casedir+'/%.0f-%.0f-target_e.vtk'%(self.azimuth, self.elevation)
        points, tri, flux, flux_abs, flux_back, flux_abs_back=flux_reader(vtkfile, casedir)
        plot_fluxmap(points, tri, flux, casedir, loc_z_rec=self.loc_z, rec_r=self.rec_r, rec_h=self.rec_h, m=self.mesh_h, n=self.mesh_circ)

    '''   
    def test_corners(self):
        datadir='/media/yewang/Data/Work/Research/Topics/yewang/NREL-raytracing/5473 details'
        cx=np.loadtxt(datadir+'/corners_5473_theoretical_global_x.csv', delimiter=',')
        cy=np.loadtxt(datadir+'/corners_5473_theoretical_global_y.csv', delimiter=',')
        cz=np.loadtxt(datadir+'/corners_5473_theoretical_global_z.csv', delimiter=',')

        #hst=np.loadtxt('./test-heliostat-facets/corners.csv', delimiter=',')
        fn='./test-heliostat-facets/90-62-primaries.vtk'
        points, poly=read_vtk(fn)
        print(np.shape(points))
        X=points[:,0]
        Y=points[:,1]
        Z=points[:,2]
        tri=poly[:,1:]
        c=np.sum(Z[tri.astype(int)], axis=1)/3.

        idx=np.array([])
        dX=np.array([])
        dY=np.array([])
        dZ=np.array([])
        diff=0
        for i in range(len(cx)):
            for j in range(4):
                A=np.r_[cx[i,j], cy[i,j], cz[i,j]]  
                dist=np.sqrt((X-A[0])**2+(Y-A[1])**2+(Z-A[2])**2)
                idx1=np.argmin(dist)
                idx=np.append(idx, idx1)
                diff+=np.min(dist)
                dx=X[idx1]-cx[i,j]
                dy=Y[idx1]-cy[i,j]
                dz=Z[idx1]-cz[i,j]
                dX=np.append(dX, dx)
                dY=np.append(dY, dy)
                dZ=np.append(dZ, dz)
                print(np.min(dist), Z[idx1]-cz[i,j])
        print('Total diff', diff)
        idx=idx.astype(int)
        plt.scatter(cx, cy, marker='o', c=cz)
        plt.scatter(X[idx], Y[idx], marker='x', c=Z[idx])
        plt.legend()
        #plt.show()
        plt.close()

        plt.scatter(X[idx], Y[idx], c=dZ, cmap='seismic')
        cbr=plt.colorbar()
        cbr.ax.tick_params(labelsize=16)
        plt.xticks(fontsize=16)
        plt.yticks(fontsize=16)
        plt.xlabel('X', fontsize=16)
        plt.ylabel('Y', fontsize=16)
        #plt.show()
        plt.close()
  
            
        #plt.scatter(hst[:,0], hst[:,1], marker='.', c=hst[:,2])
        plt.triplot(X, Y, tri)
        plt.tripcolor(X, Y, tri, facecolors=c)
        plt.scatter(cx[:,0], cy[:,0], marker='x', c=cz[:,0])
        plt.scatter(cx[:,1], cy[:,1], marker='x', c=cz[:,1])
        plt.scatter(cx[:,2], cy[:,2], marker='x', c=cz[:,2])
        plt.scatter(cx[:,3], cy[:,3], marker='x', c=cz[:,3])
        cbr=plt.colorbar()
        cbr.ax.tick_params(labelsize=16)
        plt.xticks(fontsize=16)
        plt.yticks(fontsize=16)
        plt.xlabel('X', fontsize=16)
        plt.ylabel('Y', fontsize=16)
        plt.show()
        plt.close()
    '''		
def heliostat_selection(case_id, casedir, H_pedestal):

    layout=np.loadtxt('./data/heliostats_pos_ID.csv', delimiter=',', skiprows=1)
    ID=layout[:,0]
    hst_x=layout[:,1]
    hst_y=layout[:,2]
    hst_z=np.ones(len(hst_x))*H_pedestal
    hst_pos=np.append(hst_x, (hst_y, hst_z))
    hst_pos=hst_pos.reshape(3, len(hst_x))
    hst_pos=hst_pos.T

    idx=ID==case_id
    hst0_x=hst_x[idx]
    hst0_y=hst_y[idx]
    r=np.sqrt(hst0_x**2+hst0_y**2)
    R=np.sqrt(hst_x**2+hst_y**2)

    phi=np.arctan2(hst0_y, hst0_x)
    PHI=np.arctan2(hst_y, hst_x)

    if phi<0:
        phi+=np.pi*2.
    PHI[PHI<0]=PHI[PHI<0]+np.pi*2.

    idx1=((R<1.01*r) *(R>0.8*r))*(PHI>0.95*phi)*(PHI<1.05*phi)
    hst_x=hst_x[idx1]
    hst_y=hst_y[idx1]
    hst_z=hst_z[idx1]
    hst_pos=hst_pos[idx1]
    np.savetxt(casedir+'/hst_select_%s.csv'%case_id, ID[idx1], fmt='%.0f', delimiter=',')
    #plt.plot(hst_x, hst_y, '.')
    #plt.plot(hst0_x, hst0_y, 'x')
    #plt.show()
    #plt.close()
    return hst_x, hst_y, hst_z, hst_pos


def plot_fluxmap(points, tri, flux, casedir, casename, loc_z_rec=171.035, rec_r=7.75, rec_h= 25.05609, m=31, n=60):

	X=points[:,0]
	Y=points[:,1]
	Z=points[:,2]-loc_z_rec

	THETA=np.array([])
	for i in range(len(X)):
		x=X[i]
		y=Y[i]
		if x>=0 and y>=0:
			theta=np.arcsin(x/rec_r)
		elif x>=0 and y<0:
			theta=np.pi/2.+np.arcsin(-y/rec_r)
		elif x<0 and y <0:
			theta=-np.pi/2.-np.arcsin(-y/rec_r)			
		elif x<0 and y>=0:
			theta=np.arcsin(x/rec_r)
		THETA=np.append(THETA, theta)

	idx=(Y+rec_r>0.01)
	circ=THETA*rec_r


	plt.tripcolor(X, Z, tri[:-2*31], facecolors=flux[:-2*31], cmap='jet') 
	plt.colorbar()
	plt.savefig(open(casedir+'/flux_tri.png', 'wb'), bbox_inches='tight')
	plt.close() 

	print('flux shape 0', np.shape(flux))  
	flux=(flux[::2]+flux[1::2])/2. # triangle bins combined to rectangular bins
	print('flux shape', np.shape(flux), n)
	flux=flux[:-n] #TODO why this?
	width=2.*np.pi*rec_r
	height=rec_h

	flux=flux.reshape(n,m)
	flux=flux.T
	flux=np.fliplr(flux)

	FLUX=np.array([])
	tri=tri[::2]
	idx_x=tri[:-n,0].reshape(n, m)
	idx_y=tri[:-n,1].reshape(n, m)
	idx_x=np.fliplr(idx_x.T)
	idx_y=np.fliplr(idx_y.T)


	for i in range(n):
		q=flux[:,i]
		ix=int(idx_x[0,i])
		iy=int(idx_y[0,i])
		x=points[ix,0]
		y=points[iy,1]
		if x>=0 and y>=0:
			FLUX=np.append(FLUX, q)

	for i in range(n):
		q=flux[:,i]
		ix=int(idx_x[0,i])
		iy=int(idx_y[0,i])
		x=points[ix,0]
		y=points[iy,1]
		if x>=0 and y<0:
			FLUX=np.append(FLUX, q)

	for i in range(n):
		q=flux[:,i]
		ix=int(idx_x[0,i])
		iy=int(idx_y[0,i])
		x=points[ix,0]
		y=points[iy,1]
		if x<0 and y<0:
			FLUX=np.append(FLUX, q)

	for i in range(n):
		q=flux[:,i]
		ix=int(idx_x[0,i])
		iy=int(idx_y[0,i])
		x=points[ix,0]
		y=points[iy,1]
		if x<0 and y>=0:
			FLUX=np.append(FLUX, q)
	FLUX=FLUX.reshape(n, m)
	FLUX=FLUX.T
	

	xx=np.linspace(-width/2., width/2., n+1)
	yy=np.linspace(-height/2., height/2., m+1)

	plt.pcolormesh(xx, yy, FLUX, cmap='jet')#, vmax=2400, vmin=0)
	plt.colorbar()
	plt.xlim([-width/2., width/2.])
	plt.ylim([-height/2.,height/2.])
	plt.gca().set_aspect('equal', adjustable='box')
	plt.savefig(open(casedir+'/flux_rect_%sx%s.png'%(m,n), 'wb'), bbox_inches='tight')
	#plt.show()
	plt.close()

	dx=float(width/n)
	dy=float(height/m)
	X=np.linspace(-width/2.+dx/2., width/2.-dx/2., n)
	Y=np.linspace(-height/2.+dy/2., height/2.-dy/2., m)
	XX,YY=np.meshgrid(X,Y)
	np.savetxt(casedir+'/%s_fluxmap.csv'%(casename), FLUX, fmt='%.6f', delimiter=',')
	np.savetxt(casedir+'/%s_xx.csv'%(casename), XX, fmt='%.2f', delimiter=',')
	np.savetxt(casedir+'/%s_yy.csv'%(casename), YY, fmt='%.2f', delimiter=',')	

if __name__ == '__main__':
	unittest.main()

