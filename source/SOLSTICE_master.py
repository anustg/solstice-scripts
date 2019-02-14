from spectral_data import SolarSpectrum, MirrorRhoSpectrum
from get_raw import proces_raw_results

import subprocess
import sys
import matplotlib
import numpy as np
import matplotlib.cm as cm
import matplotlib.mlab as mlab
import matplotlib.pyplot as plt
import os
import math
import time


class SolsticeScene:

    def __init__(self, receiver, folder, casename, pmfile, annual=False):
        '''
        Import and define the parameters
        Arguments:
        receiver - str, 'flat' or 'blade' or 'STL'
        folder - str, the directory of saving the case
        casename - str, define that name of the case, for saving files
        pmfile - str, the directory of the parameter input csv file
        annual - bool, do an annual performance simulation or not
        '''
        self.folder=folder
        self.casename=casename
        self.annual=annual

        self.receiver=receiver

        # Solar Related parameters
        pm=np.loadtxt(pmfile, dtype=str, delimiter=',', skiprows=9)
        
        self.azimuth=pm[0,4].astype(float)
        self.zenith=pm[1,4].astype(float)
        self.num_rays=int(pm[5,4].astype(float))

        self.spectral=pm[6,4].astype(bool)
        self.air=pm[7,4].astype(float)

        self.sunshape=pm[2,4]
        self.sunsize=pm[3,4].astype(float)
        self.dni=pm[4,4].astype(float)


        # heliostat related parameters
        pm=np.loadtxt(pmfile, dtype=str, delimiter=',', skiprows=19)
        self.fieldname=pm[0,4]
        self.mirror_reflectivity=pm[1,4].astype(float)
        self.slope=pm[2,4].astype(float)
        self.hst_dir=pm[3,4]
        self.hst_w=pm[4,4].astype(float)
        self.hst_h=pm[5,4].astype(float)
        self.tower_h=pm[6,4].astype(float)
        self.tower_r=pm[7,4].astype(float)
        self.tower_slice=int(pm[8,4].astype(float))

        # receiver related parameters
        pm=np.loadtxt(pmfile, dtype=str, delimiter=',', skiprows=31)
        self.rec_w=pm[0,4].astype(float)
        self.rec_h=pm[1,4].astype(float)
        self.absorptivity=pm[2,4].astype(float)
        self.rec_x=pm[3,4].astype(float)
        self.rec_y=pm[4,4].astype(float)
        self.rec_z=pm[5,4].astype(float)
        self.rec_tilt=pm[6,4].astype(float)
        self.rec_slice=int(pm[7,4].astype(float))

    def gen_YAML(self):

        '''
        Generate YAML file 
        '''
        #TODO try pyYAML?

        # OPEN the input YAML file
        fd = open('%s/%s.yaml'%(self.folder, self.casename), 'w')


        # -------------------------------- Creat the specturm properties -----------------------
        if self.spectral:

            I_sun=SolarSpectrum()
            # CREATE the spectrum for the sun
            fd.write('- spectrum: &%s  \n' % 'solar_spectrum')
            for i in range(0,len(I_sun)-1):
	            fd.write('  - {wavelength: %15.8e, data: %15.8e }\n' % (float(I_sun[i][0]),float(I_sun[i][1])) ) 
            i = len(I_sun)-1
            fd.write('  - {wavelength: %15.8e, data: %15.8e }\n' % (float(I_sun[i][0]),float(I_sun[i][1])) )
            fd.write('\n')

            # CREATE the spectrum for the reflectivity (mirror)
            mirror_rho= MirrorRhoSpectrum()
            mirror_ref=mirror_rho
            for i in range(0,len(mirror_rho)):
                mirror_ref[i][0] = mirror_rho[len(mirror_rho)-1-i][0]/1000.
                mirror_ref[i][1] = mirror_rho[len(mirror_rho)-1-i][1]/100.
            mirror_ref.append([4,0.9])
            fd.write('- spectrum: &%s  \n' % 'ref_mirror')
            for i in range(0,len(mirror_ref)-1):
	            fd.write('  - {wavelength: %15.8e, data: %15.8e }\n' % (float(mirror_ref[i][0]),float(mirror_ref[i][1])) ) 
            i = len(mirror_ref)-1
            fd.write('  - {wavelength: %15.8e, data: %15.8e }\n' % (float(mirror_ref[i][0]),float(mirror_ref[i][1])) )
            fd.write('\n')

        if self.air>1e-99:
            I_sun=SolarSpectrum()
            # CREATE the spectrum for the air extinction coefficient
            ke_air = I_sun
            for i in range(0,len(I_sun)):
                ke_air[i][1] = self.air 
  
            fd.write('- spectrum: &%s  \n' % 'air_kext')
            for i in range(0,len(ke_air)-1):
	            fd.write('  - {wavelength: %15.8e, data: %15.8e }\n' % (float(ke_air[i][0]),float(ke_air[i][1])) ) 
            i = len(ke_air)-1
            fd.write('  - {wavelength: %15.8e, data: %15.8e }\n' % (float(ke_air[i][0]),float(ke_air[i][1])) )
            fd.write('\n')


        # ----------------------------- CREATION OF THE SUN ------------------------
        #
        # CREATE The sun
        if self.spectral:
            fd.write('- sun: {dni: %15.8e, spectrum: *%s, %s: {half_angle: %6.4f}}\n' % (float(self.dni),'solar_spectrum', self.sunshape,self.sunsize) )
        else:
            fd.write('- sun: {dni: %15.8e, %s: {half_angle: %6.4f}}\n' % (float(self.dni), self.sunshape,self.sunsize) )


        # ----------------------------- CREATION OF THE ATMOSPHERE ------------------
        #
        if self.air>1e-99:
            fd.write('- atmosphere: {extinction: *%s}\n' % 'air_kext' )
        fd.write('\n')

        # ------------------------------CREATION OF MATERIALS ----------------------
        #
        # CREATE an occultant material
        r_f = 0. # front
        r_b = 0. # and back reflectivity
        fd.write('- material: &%s\n' % 'material_black')
        fd.write('   front:\n')
        fd.write('     matte: {reflectivity: %6.4f }\n' % float(r_f))    
        fd.write('   back:\n')
        fd.write('     matte: {reflectivity: %6.4f }\n' % float(r_b))
        fd.write('\n')
        #
        # CREATE a material for the target
        r_f = 1.-self.absorptivity # front
        r_b = 0. # and back reflectivity
        fd.write('- material: &%s\n' % 'material_target')
        fd.write('   front:\n')
        fd.write('     matte: {reflectivity: %6.4f }\n' % float(r_f))    
        fd.write('   back:\n')
        fd.write('     matte: {reflectivity: %6.4f }\n' % float(r_b))
        fd.write('\n')
        #
        # CREATE a specular material
        slope_error = self.slope 
        r_b = 0. # and back reflectivity
        fd.write('- material: &%s\n' % 'material_mirror')
        fd.write('   front:\n')
        if self.spectral:
            fd.write('     mirror: {reflectivity: *%s, slope_error: %15.8e }\n' % ('ref_mirror',float(slope_error) ) ) 
        else:
            fd.write('     mirror: {reflectivity: %6.4f, slope_error: %15.8e }\n' % (self.mirror_reflectivity,float(slope_error) ) ) 
  
        fd.write('   back:\n')
        fd.write('     matte: {reflectivity: %6.4f }\n' % float(r_b)) 
        fd.write('\n')
        #
        # CREATE a material for the Large target to compute spillage
        fd.write('- material: &%s\n' % 'material_virtual')
        fd.write('   virtual:\n')
        fd.write('\n')


        #----------------------- CREATION OF GEOMETRIES -----------------------------
        #
        #   Tower Geometry 
        #- cylindrical shape
        slices = self.tower_slice # slices for the envelop circle
        fd.write('- geometry: &%s\n' % 'tower_g' )
        fd.write('  - material: *%s\n' % 'material_black' )
        #fd.write('    transform: { translation: %s, rotation: %s }\n' % ([0, 0, h_tow*0.5], [0, 90, 0]) )
        fd.write('    cylinder: {height: %7.3f, radius: %7.3f, slices: %d }\n' % (self.tower_h, self.tower_r, slices) )
        fd.write('\n')
        #


        #    Receiver Geometry
        #
        if self.receiver =='flat':
            #
            slices = self.rec_slice # slices for the fluxmap
            #
            pts = [ [-self.rec_w*0.5, -self.rec_h*0.5], [-self.rec_w*0.5, self.rec_h*0.5], [self.rec_w*0.5, self.rec_h*0.5], [self.rec_w*0.5,-self.rec_h*0.5] ]
            #
            fd.write('- geometry: &%s\n' % 'target_g' )
            fd.write('  - material: *%s\n' % 'material_target' )
            fd.write('    plane: \n')
            fd.write('      clip: \n')    
            fd.write('      - operation: AND \n')
            fd.write('        vertices: %s\n' % pts)
            fd.write('      slices: %d\n' % slices ) 
            fd.write('\n')
            #
        elif self.receiver =='blade':

            # the back wall flat plane
            self.rec_h = self.rec_h # target height
            self.rec_w = self.rec_w # target width

            pts = [ [-self.rec_w*0.5, -self.rec_h*0.5], [-self.rec_w*0.5, self.rec_h*0.5], [self.rec_w*0.5, self.rec_h*0.5], [self.rec_w*0.5,-self.rec_h*0.5] ]
            #
            fd.write('- geometry: &%s\n' % 'target_g' )
            fd.write('  - material: *%s\n' % 'material_target' )
            fd.write('    plane: \n')
            fd.write('      clip: \n')    
            fd.write('      - operation: AND \n')
            fd.write('        vertices: %s\n' % pts)
            fd.write('      slices: %d\n' % slices ) 
            fd.write('\n')
        
            nb = 4 # number of blades >= 3
            Tb = 0.05 # blade thickness, m
            tilt = 90. #tilt angle of the blades, degree 
            Db = 0.86*1.157 #Depth of the blade (bottom)
            Wb = 8.5 # Width of the blade (and of the receiver)
            Hb = 8.5 # Height of the bladed receiver

            pi = 4.*np.arctan(1.) 
            tiltr =  tilt*pi/180. #tilt in radians
            # Depth of the top blade
            Lb = Db + Tb/np.tan(tiltr) # Depth of the top blade Lb = Db + Tb/tan(tilt) 
            # space between blades
            space = (Hb-nb*Tb)/(nb-1.) 

            #* Geometry */
            #print_section("Geometry") 

            fd.write("- geometry: &blade_top \n")
            fd.write("    - material: *material_target\n")
            trans_y = Lb*np.sin(tiltr)*0.5 
            trans_z = Tb
            fd.write("      transform: {translation: [0,  %s,  %s], rotation: [ %s, 0, 0]}\n" % (trans_y,trans_z,-90+tilt))  #Tb-Lb*cos(tiltr)*0.5
            fd.write("      plane : {clip: [{operation: AND, vertices: [[ %s,  %s],[ %s,  %s],[ %s,  %s],[ %s,  %s]]}], slices: 64}\n"%( 
            -Wb*0.5,-Lb*0.5, -Wb*0.5,Lb*0.5, Wb*0.5,Lb*0.5, Wb*0.5,-Lb*0.5))

            fd.write("- geometry: &blade_bottom \n") 
            fd.write("    - material: *material_target\n") 
            trans_y = Db*np.sin(tiltr)*0.5  
            fd.write("      transform: {translation: [0,  %s,  %s], rotation: [ %s, 0, 0]}\n"% (trans_y ,0.,90+tilt))  #-Db*cos(tiltr)*0.5
            fd.write("      plane : {clip: [{operation: AND, vertices: [[ %s,  %s],[ %s,  %s],[ %s,  %s],[ %s,  %s]]}], slices: 64}\n"%(
            -Wb*0.5,-Db*0.5, -Wb*0.5,Db*0.5, Wb*0.5,Db*0.5, Wb*0.5,-Db*0.5)) 

            fd.write("- geometry: &blade_front \n") 
            fd.write("    - material: *material_target\n") 
            trans_y = Db*np.sin(tiltr)+Tb*np.cos(tiltr)*0.5  
            trans_z = Tb- 0.5*Tb/np.sin(tiltr) - Db*np.cos(tiltr)*0.5  
            fd.write("      transform: {translation: [0, %s,  %s], rotation: [ %s, 0, 0]}\n"%( trans_y, trans_z, 180+tilt))  #trans_z + Tb*0.5
            fd.write("      plane : {clip: [{operation: AND, vertices: [[ %s,  %s],[ %s,  %s],[ %s,  %s],[ %s,  %s]]}], slices: 32}\n"%( 
            -Wb*0.5,-Tb*0.5, -Wb*0.5,Tb*0.5, Wb*0.5,Tb*0.5, Wb*0.5,-Tb*0.5) )

            fd.write("- geometry: &blade_end_xmin \n") 
            fd.write("    - material: *material_target\n") 
            trans_y = Lb*np.sin(tiltr)*0.5  
            trans_z = Tb  - 0.5*Tb/np.sin(tiltr) 
            fd.write("      transform: {translation: [ %s,  %s,  %s], rotation: [ %s, -90, 90]}\n"%(-Wb*0.5,trans_y,trans_z, -90+tilt) )
            fd.write("      plane : {clip: [{operation: AND, vertices: [[ %s,  %s],[ %s,  %s],[ %s,  %s],[ %s,  %s]]}], slices: 32}\n"%( 
            -Db*0.5-Tb/np.tan(tiltr),-Tb*0.5, -Db*0.5,Tb*0.5, Db*0.5,Tb*0.5, Db*0.5,-Tb*0.5) )

            fd.write("- geometry: &blade_end_xmax \n") 
            fd.write("    - material: *material_target\n") 
            fd.write("      transform: {translation: [ %s,  %s,  %s], rotation: [ %s,90, 90]}\n"%(Wb*0.5,trans_y,trans_z, -90+tilt) )
            fd.write("      plane : {clip: [{operation: AND, vertices: [[ %s,  %s],[ %s,  %s],[ %s,  %s],[ %s,  %s]]}], slices: 32}\n"%( 
            -Db*0.5,-Tb*0.5, -Db*0.5-Tb/np.tan(tiltr),Tb*0.5, Db*0.5,Tb*0.5, Db*0.5,-Tb*0.5) )

            fd.write("- geometry: &space \n") 
            fd.write("    - material: *material_target\n") 
            fd.write("      transform: {translation: [0, 0,  %s], rotation: [-90, 0, 0]}\n"%( Db*np.cos(tiltr)*0.5+Tb/np.sin(tiltr)+space*0.5) )
            fd.write("      plane : {clip: [{operation: AND, vertices: [[ %s,  %s],[ %s,  %s],[ %s,  %s],[ %s,  %s]]}], slices: 128}\n"%( 
            -Wb*0.5,-space*0.5, -Wb*0.5,space*0.5, Wb*0.5,space*0.5, Wb*0.5,-space*0.5) )
  

        # IMPORT CSV file for the heliostat positions
        hst_info=np.loadtxt(self.hst_dir,delimiter=',', skiprows=2)
        hst_x=hst_info[:,0]
        hst_y=hst_info[:,1]
        hst_z=hst_info[:,2]

        foc=hst_info[:,3]

        aim_x=hst_info[:,4]
        aim_y=hst_info[:,5]
        aim_z=hst_info[:,6]

        #foc = np.sqrt((hst_x-rec_x)**2+ (hst_y-rec_y)**2+(hst_z-rec_z)**2) # ideal focus

        h_hst = self.hst_h # heliostat height
        w_hst = self.hst_w # heliostat width
        slices = 4 # slices for the envelop circle
        pts_hst = [ [-w_hst*0.5, -h_hst*0.5], [-w_hst*0.5, h_hst*0.5], [w_hst*0.5, h_hst*0.5], [w_hst*0.5,-h_hst*0.5] ]

        # CREATE a reflective facet (mirror) used in the PS10 heliostat template
        for i in range(0,len(foc)):
            name_hst_g = 'hst_g_'+str(i)
            fd.write('- geometry: &%s\n' % name_hst_g )
            fd.write('  - material: *%s\n' % 'material_mirror' )
        #    fd.write('    transform: { translation: %s, rotation: %s }\n' % ([hst_x[i], hst_y[i], hst_z[i]], [0, 0, 0]) )
            fd.write('    parabol: \n')
            fd.write('      focal: %s\n' % foc[i]) 
            fd.write('      clip: \n')    
            fd.write('      - operation: AND \n')
            fd.write('        vertices: %s\n' % pts_hst)
            fd.write('      slices: %d\n' % slices )  

        # CREATE the pylon "pylon_g" geometry cylindrical shape
        h_pyl = 0.001 # pylon height
        r_pyl = 0.2 # pylon radius
        slices = 4 # slices for the envelop circle
        fd.write('- geometry: &%s\n' % 'pylon_g' )
        fd.write('  - material: *%s\n' % 'material_black' )
        fd.write('    transform: { translation: %s, rotation: %s }\n' % ([0, 0, -h_pyl*3], [0, 90, 0]) )
        fd.write('    cylinder: {height: %7.3f, radius: %7.3f, slices: %d }\n' % (h_pyl,r_pyl,slices) )
        #
            
        # -------------------------------------- CREATE THE TEMPLATES using the geometries

        # CREATE the heliostat templates

        for i in range(0,len(foc)):    
            name_hst_t = 'hst_t_'+str(i)
            fd.write('- template: &%s\n' % name_hst_t )
            name_hst_n = 'hst_'+ str(i)
            fd.write('    name: %s\n' % name_hst_n )
            fd.write('    primary: 0\n' )    
            fd.write('    geometry: *pylon_g\n')
            fd.write('    children: \n' )
            fd.write('    - name: pivot\n')
            fd.write('      zx_pivot: {target: {position: %s}} \n' % ([aim_x[i],aim_y[i],aim_z[i]]) )
            fd.write('      children: \n')
            fd.write('      - name: reflect_surface\n')
            fd.write('        primary: 1\n')
            fd.write('        transform: {rotation: [-90,0,0]} \n')    
            name_hst_g = 'hst_g_'+str(i)
            fd.write('        geometry: *%s\n' % name_hst_g )


        if self.receiver=='blade':

            #/* Template for the bladed receiver */
            fd.write("- template: &bladed_rcv_t\n") 
            fd.write("    name: bladed_receiver_t\n")    # metal structure
            fd.write("    primary: 0\n") 
            fd.write("    geometry:\n") 
            fd.write("    - material: *material_virtual\n") 
            fd.write("      transform: {translation: [0,  %s, 0], rotation: [0, 0, 0]}\n"%0.) #,-Db) 
            # half width of the square virtual surface
            VV = 4.*Hb  
            fd.write("      plane: {clip: [{operation: AND, vertices: [[ %s,  %s],[ %s,  %s],[ %s,  %s],[ %s,  %s]]}], slices: 128}\n"%( 
            -VV,-VV, -VV,VV, VV,VV, VV,-VV) )
            fd.write("    children:\n") 
            i=0
            while i< nb:
                i+=1
                fd.write("    - name: blade_bottom_%d\n"%i)   # First blade at Z minimum
                fd.write("      primary: 0\n") 
                fd.write("      transform: {translation: [0, 0,  %s]}\n"%(-Hb*0.5+ (float(i)*(Tb/np.sin(tiltr)+space))))
                fd.write("      geometry: *blade_bottom\n") 
                fd.write("    - name: blade_top_%d\n"%i)   # First blade at Z minimum
                fd.write("      primary: 0\n") 
                fd.write("      transform: {translation: [0, 0,  %s]}\n"%(-Hb*0.5+ (float(i))*(Tb/np.sin(tiltr)+space)))
                fd.write("      geometry: *blade_top\n") 
                fd.write("    - name: blade_front_%d\n"%i)   # First blade at Z minimum
                fd.write("      primary: 0\n") 
                fd.write("      transform: {translation: [0, 0,  %s]}\n"%(-Hb*0.5+ (float(i))*(Tb/np.sin(tiltr)+space)))
                fd.write("      geometry: *blade_front\n") 
                fd.write("    - name: blade_end_xmin_%d\n"%i)   # First blade at Z minimum
                fd.write("      primary: 0\n") 
                fd.write("      transform: {translation: [0, 0,  %s]}\n"%(-Hb*0.5+ (float(i))*(Tb/np.sin(tiltr)+space)))
                fd.write("      geometry: *blade_end_xmin\n") 
                fd.write("    - name: blade_end_xmax_%d\n"%i)   # First blade at Z minimum
                fd.write("      primary: 0\n") 
                fd.write("      transform: {translation: [0, 0,  %s]}\n"%(-Hb*0.5+ (float(i))*(Tb/np.sin(tiltr)+space)))
                fd.write("      geometry: *blade_end_xmax\n") 
                if i<nb:
                    fd.write("    - name: space_%d\n"%i)   # First blade at Z minimum
                    fd.write("      primary: 0\n") 
                    fd.write("      transform: {translation: [0, 0,  %s]}\n"%(-Hb*0.5+ ((float(i))*(Tb/np.sin(tiltr)+space))))
                    fd.write("      geometry: *space\n") 

            


        #
        #
        # -------------------------------------- CREATE THE ENTITIES using the geometries or the templates
        #
        ### CREATE entities from geometries: specifying if primary reflector
        #
        if self.receiver=='flat':
            # CREATE a target entity from "target_g" geometry (primary = 0)
            fd.write('\n- entity:\n')
            fd.write('    name: target_e\n')
            fd.write('    primary: 0\n')
            fd.write('    transform: { translation: %s, rotation: %s }\n' % ([self.rec_x, self.rec_y, self.rec_z], [-90-self.rec_tilt, 0, 0]) )
            fd.write('    geometry: *%s\n' % 'target_g')


            pts = [ [-self.rec_w*10., -self.rec_h*10.], [-self.rec_w*10., self.rec_h*10.], [self.rec_w*10., self.rec_h*10.], [self.rec_w*10.,-self.rec_h*10.] ]
            slices = 4
            # CREATE a target entity from "target_g" geometry (primary = 0)
            fd.write('\n- entity:\n')
            fd.write('    name: virtual_target_e\n')
            fd.write('    primary: 0\n')
            fd.write('    transform: { translation: %s, rotation: %s }\n' % ([self.rec_x, self.rec_y-5., self.rec_z], [-90-self.rec_tilt, 0, 0]))
            fd.write('    geometry: \n' )
            fd.write('      - material: *%s\n' % 'material_virtual' )

            fd.write('        plane: \n')
            fd.write('          clip: \n')    
            fd.write('          - operation: AND \n')
            fd.write('            vertices: %s\n' % pts)
            fd.write('          slices: %d\n' % slices ) 


        elif self.receiver=='blade':

            # CREATE a target entity from "target_g" geometry (primary = 0)
            fd.write('\n- entity:\n')
            fd.write('    name: target_e\n')
            fd.write('    primary: 0\n')
            fd.write('    transform: { translation: %s, rotation: %s }\n' % ([self.rec_x, self.rec_y, self.rec_z], [-90-self.rec_tilt, 0, 0]) )
            fd.write('    geometry: *%s\n' % 'target_g')


            pts = [ [-self.rec_w*10., -self.rec_h*10.], [-self.rec_w*10., self.rec_h*10.], [self.rec_w*10., self.rec_h*10.], [self.rec_w*10.,-self.rec_h*10.] ]
            slices = 4
            # CREATE a target entity from "target_g" geometry (primary = 0)
            fd.write('\n- entity:\n')
            fd.write('    name: virtual_target_e\n')
            fd.write('    primary: 0\n')
            fd.write('    transform: { translation: %s, rotation: %s }\n' % ([self.rec_x, self.rec_y-5., self.rec_z], [-90-self.rec_tilt, 0, 0]))
            fd.write('    geometry: \n' )
            fd.write('      - material: *%s\n' % 'material_virtual' )
            fd.write('        plane: \n')
            fd.write('          clip: \n')    
            fd.write('          - operation: AND \n')
            fd.write('            vertices: %s\n' % pts)
            fd.write('          slices: %d\n' % slices ) 

            fd.write("- entity: {name: bladed_receiver_e, children: [*bladed_rcv_t], transform: {translation: [%s, %s, %s]}}\n"%(self.rec_x, self.rec_y, self.rec_z))



        elif self.receiver=='STL':
            #/*entity from STL file */

            fd.write("- entity: \n")
            fd.write("    name: STL_receiver_e\n")
            fd.write("    primary: 0\n")
            fd.write("    transform: {translation: [0, 0, %s], rotation: [90, 0, 0]}\n"% (self.rec_z))
            fd.write("    geometry:\n")
            fd.write("    - material: *material_target\n")
            fd.write("      transform: {translation: [0, 0, 0], rotation: [0, 0, 0]}\n")
            fd.write("      stl : {path: %s/plane_10x10_5cm.stl }  \n"%(self.folder))



        #
        # CREATE a tower entity from "tow" geometry (primary = 0)
        fd.write('\n- entity:\n')
        fd.write('    name: tower_e\n')
        fd.write('    primary: 0\n' )
        fd.write('    transform: { translation: %s, rotation: %s }\n' % ([0, -self.tower_r, self.tower_h*0.5], [0, 0, 0]) )
        fd.write('    geometry: *%s\n' % 'tower_g')


        #
        #
        ### CREATE entities from templates
        #
        #
        for i in range(0,len(foc)):
            name_e ='H_'+str(i)
            name_hst_t = 'hst_t_'+str(i)
            fd.write('\n- entity:\n')
            fd.write('    name: %s\n' % name_e)
            fd.write('    transform: { translation: %s, rotation: %s }\n' % ([hst_x[i], hst_y[i], hst_z[i]], [0, 0, 0]) )
            fd.write('    children: [ *%s ]\n' % name_hst_t)

        # END OF THE YAML INPUT FILE
        fd.close() 


        # WRITE THE YAML FILE FOR THE RECEIVERS
        fd = open('%s/%s-rcv.yaml'%(self.folder,self.casename), 'w')

        if self.receiver=='flat':
            fd.write('- name: target_e \n' )
            fd.write('  side: %s \n' % 'FRONT_AND_BACK')
            fd.write('  per_primitive: %s \n' % 'INCOMING_AND_ABSORBED')

            fd.write('- name: virtual_target_e \n' )
            fd.write('  side: %s \n' % 'FRONT')
            fd.write('  per_primitive: %s \n' % 'INCOMING')


        elif self.receiver=='blade':
            fd.write('- name: target_e \n' )
            fd.write('  side: %s \n' % 'FRONT_AND_BACK')
            fd.write('  per_primitive: %s \n' % 'INCOMING_AND_ABSORBED')

            i=0
            while i<nb:
                i+=1
                fd.write("- {name: bladed_receiver_e.bladed_receiver_t.blade_bottom_%d, side: FRONT, per_primitive: ABSORBED}\n"%i)
                fd.write("- {name: bladed_receiver_e.bladed_receiver_t.blade_top_%d, side: FRONT, per_primitive: ABSORBED}\n"%i)
                fd.write("- {name: bladed_receiver_e.bladed_receiver_t.blade_front_%d, side: FRONT, per_primitive: ABSORBED}\n"%i)
                fd.write("- {name: bladed_receiver_e.bladed_receiver_t.blade_end_xmin_%d, side: FRONT, per_primitive: ABSORBED}\n"%i)
                fd.write("- {name: bladed_receiver_e.bladed_receiver_t.blade_end_xmax_%d, side: FRONT, per_primitive: ABSORBED}\n"%i)
                if i<nb:
                    fd.write("- {name: bladed_receiver_e.bladed_receiver_t.space_%d, side: FRONT, per_primitive: ABSORBED}\n"%i)

            fd.write("- {name: bladed_receiver_e.bladed_receiver_t, side: FRONT, per_primitive: ABSORBED}\n")
        

        elif self.receiver=='STL':
            
            #/* Receivers */
            fd.write("- {name: STL_receiver_e, side: FRONT_AND_BACK, per_primitive: ABSORBED}\n")




        # -------------------------------------- END OF THE YAML PARSER WRITTING
        fd.close()

    def runSOLSTICE(self, savefile, azi=0., zenith=0., view=False):
        '''
        run SOLSTICE with the corresponding sun position  
        and postprocessing the result   
        azi: from East to North
        zenith: 0 is the horizontal (in Solstice)
        view - if check it in paraview   
        '''


        if not self.annual:
            azi=self.azimuth
            zenith=self.zenith

       
        #os.system('source /home/yewang/SOLSTICE/Solstice-0.8.1-GNU-Linux64/etc/solstice.profile')    
        os.system('solstice -D%s,%s -v -n %s -R %s/%s-rcv.yaml -fo %s/simul %s/%s.yaml'%(azi,zenith,self.num_rays, self.folder, self.casename,savefile,self.folder,self.casename))


        if view:
            os.system('solstice -D%s,%s -g format=obj:split=geometry -fo %s/geom %s/%s.yaml'%(azi, zenith,savefile, self.folder, self.casename))
            os.system('solstice -D%s,%s -q -n 100 -R %s/%s-rcv.yaml -p default %s/%s.yaml > %s/solpaths'%(azi, zenith, self.folder, self.casename, self.folder, self.casename, savefile ))

            # postprocessing in C (provided by Cyril Caliot)
            #Read "simul" results and produce a text file with the raw results
            os.system('gcc ./postprocessing/solppraw.c -o %s/solppraw'%savefile)
            os.system('%s/solppraw %s/simul'%(savefile, savefile))

            #Read "simul" results and produce receiver files (.vtk) of incoming and/or absorbed solar flux per-primitive
            os.system('gcc ./postprocessing/solmaps.c -o %s/solmaps'%savefile)
            os.system('%s/solmaps %s/simul'%(savefile,savefile))

            #Read "geom" and "simul" file results and produce primaries and receivers files (.vtk), and .obj geometry files
            os.system('gcc ./postprocessing/solpp.c -o %s/solpp'%savefile)
            os.system('%s/solpp %s/geom %s/simul'%(savefile,savefile, savefile))

            #Read "solpaths" file and produce readable file (.vtk) by paraview to visualize the ray paths
            os.system('gcc ./postprocessing/solpaths.c -o %s/solpath'%savefile)
            os.system('%s/solpath %s/solpaths'%(savefile, savefile))

            os.system('mv *vtk %s'%savefile)
            os.system('mv *obj %s'%savefile)
            os.system('mv *txt %s'%savefile)

            rawfile='%s/simul'%savefile
            eta_tot=proces_raw_results(rawfile,savefile)

        else:
            #Read "simul" results and produce a text file with the raw results
            os.system('gcc ./postprocessing/solppraw.c -o %s/solppraw'%savefile)
            os.system('%s/solppraw %s/simul'%(savefile,savefile))
            os.system('mv *txt %s'%savefile)
            rawfile='%s/simul'%savefile
            eta_tot=proces_raw_results(rawfile,savefile)

        return eta_tot


class AnnualPerformance:
    
    def __init__(self, folder, casename, pmfile):
       
        self.folder=folder
        # define the geometry configuration in SOLSTICE
        self.scene=SolsticeScene('flat',self.folder, casename, pmfile, annual=True) 
        self.scene.gen_YAML()
        self.sun_position()

    def sun_position(self):
        # sun_position that defined in SolarTherm
        self.AZI=np.linspace(0., 180., 7)
        self.ELE=np.linspace(0., 90., 4)

        # convert sun position into SOLSTICE convention
        self.AZI_c=90.-self.AZI
        self.ZENITH_c=self.ELE
        
        for i in xrange(len(self.AZI_c)):
            if self.AZI_c[i]<0.:
                self.AZI_c[i]+=360.


    def run(self):
        '''
        SINGLE Process
        '''
        table=np.array([])
        for j in xrange(len(self.ELE)):
            ele=self.ELE[j]
            zenith_c=self.ZENITH_c[j]
               
            for i in xrange(len(self.AZI)):
                azi=self.AZI[i]
                azi_c=self.AZI_c[i]
 
                savefile=self.folder+'/azi%s_ele_%s'%(int(azi), int(ele))
                if not os.path.exists(savefile):
                    os.makedirs(savefile)

                print ''
                print 'CASE (%s/%s)'%((j*len(self.AZI)+i+1), (len(self.AZI))*(len(self.ELE)))
                print 'azi%s_ele%s'%(int(azi), int(ele))

                eta=self.scene.runSOLSTICE(savefile=savefile, azi=azi_c, zenith=zenith_c, view=False)
                table=np.append(table, eta)

                print eta

        table=table.reshape(len(self.ELE), len(self.AZI))
        full_table=np.hstack((table, np.fliplr(table)[:, 1:]))
        print np.shape(table)

        ELE=self.ELE.reshape(len(self.ELE), 1)
        AZI=np.append(9999, self.AZI)
        AZI=AZI.reshape(1, len(AZI))

        table=np.hstack((ELE,table))
        table=np.vstack((AZI, table))

        np.savetxt(self.folder+'/lookup_table_ori.csv', table, fmt='%.4f', delimiter=',')

        # apply symetry to make the full look up table
        print np.shape(full_table)
        full_table=np.hstack((ELE,full_table))
        AZI2=np.append(AZI, AZI[0,1:-1]+360.-AZI[0,-2])
        AZI2=AZI2.reshape(1, len(AZI2))
        print np.shape(AZI2)
        full_table=np.vstack((AZI2, full_table))

        np.savetxt(self.folder+'/lookup_table_full.csv', full_table, fmt='%.4f', delimiter=',')        
                
        
       

if __name__=='__main__':

    #======================================================
    #
    #   Try annual performance or just run a single case?
    #
    case='single' # Case = 'annual' or '1p'
    #
    #=======================================================
    
    if case=='single':
     
        casename='demo'
        casefolder='./demo'
        if not os.path.exists(casefolder):
            os.makedirs(casefolder)
        pmfile='./source/set_parameter.csv'
        receiver='flat'
        scene=SolsticeScene(receiver, casefolder, casename, pmfile, annual=False)
        scene.gen_YAML()
        eta=scene.runSOLSTICE(savefile=casefolder, view=True)
        print eta
    
    elif case=='annual':

        init_time=time.time() 
        casename='demo'
        casefolder='./demo'
        if not os.path.exists(casefolder):
            os.makedirs(casefolder)
        pmfile='./source/set_parameter.csv'
        an=AnnualPerformance(casefolder,casename, pmfile)
        an.run()
        total_time=time.time()-init_time
        print ''
        print 'FINISH'
        print 'Total Time (min)', total_time/60.








