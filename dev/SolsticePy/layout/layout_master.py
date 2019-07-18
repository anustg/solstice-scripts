import numpy as N


class LayoutGeneration:
    '''
    Generate a heliostat layout
    version 1: 
               (1) size is based on the total reflector area 
               (2) shape options: polar or surrounding
               (3) layout options: radial staggered, mueen or corn field
               (4) performance optimisation is not considered
    '''


    def __init__(self,pmfile):
        '''
        read the layout related parameters
        now - they are read from the csv file
        future - they are assigned by SolarTherm

        '''
        # Solar Related parameters
        #pm=N.loadtxt(pmfile, dtype=str, delimiter=',', skiprows=9)        
        #self.azimuth=pm[0,4].astype(float)
        #self.zenith=pm[1,4].astype(float)
        #self.num_rays=int(pm[5,4].astype(float))
        #self.spectral=pm[6,4].astype(bool)
        #self.air=pm[7,4].astype(float)
        #self.sunshape=pm[2,4]
        #self.sunsize=pm[3,4].astype(float)
        #self.dni=pm[4,4].astype(float)


        # heliostat related parameters
        pm=N.loadtxt(pmfile, dtype=str, delimiter=',', skiprows=19)
        for i in xrange(len(pm)):
            for j in xrange(len(pm[i])):
                pm[i,j]=pm[i,j].replace('"','')
        self.field_layout=pm[0,4]
        self.field_shape=pm[1,4]
        self.field_area=pm[2,4].astype(float)
        #self.mirror_reflectivity=pm[3,4].astype(float)
        #self.slope=pm[4,4].astype(float)
        self.hst_dir=pm[5,4]
        self.hst_w=pm[6,4].astype(float)
        self.hst_h=pm[7,4].astype(float)
        self.hst_z=pm[8,4].astype(float)
        self.tower_h=pm[9,4].astype(float)
        self.tower_r=pm[10,4].astype(float)
        self.tower_slice=int(pm[11,4].astype(float))

        # receiver related parameters
        #pm=N.loadtxt(pmfile, dtype=str, delimiter=',', skiprows=33)
        #self.rec_w=pm[0,4].astype(float)
        #self.rec_h=pm[1,4].astype(float)
        #self.absorptivity=pm[2,4].astype(float)
        #self.rec_x=pm[3,4].astype(float)
        #self.rec_y=pm[4,4].astype(float)
        #self.rec_z=pm[5,4].astype(float)
        #self.rec_tilt=pm[6,4].astype(float)
        #self.rec_slice=int(pm[7,4].astype(float))

    def radial_stagger(self):
        pass

    def mueen(self):
        pass

    def corn_field(self):
        pass


