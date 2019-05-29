import numpy as N
#from tracer.models.heliostat_field import solar_vector
import matplotlib.pyplot as plt

class Field:

    def __init__(self, hstpos, aimpos, azimuth, zenith):
        '''

        hstpos: (n,3) array, the heliostats position
        aimpos: (n,3) array, the aiming points of each heliostats
        '''
        self.hstpos=hstpos
        self.aimpos=aimpos
        self.sun_vec=self.get_solar_vector(azimuth, zenith)
        self.get_normal()
        self.cosine_factor=self.get_cosine()

    def get_solar_vector(self, azimuth, zenith):
        """
        Calculate the solar vector using elevation and azimuth.

        Arguments:
        azimuth - the sun's azimuth, in deg, 
	        from South increasing towards to the West
        zenith - angle created between the solar vector and the Z axis, in deg.

        Returns: a 3-component 1D array with the solar vector.
        """
        #copied from tracer.models.heliostat_field import solar_vector
        #TODO change it to Solstice convetion if necessary

        azimuth*=N.pi/180.
        zenith*=N.pi/180.

        sun_z = N.cos(zenith)
        sun_y=-N.sin(zenith)*N.cos(azimuth)
        sun_x=-N.sin(zenith)*N.sin(azimuth)
        sun_vec = N.r_[sun_x, sun_y,sun_z] 

        return sun_vec

    def get_normal(self):
        tower_vec=self.aimpos-self.hstpos 
        tower_vec/=N.sqrt(N.sum(tower_vec**2, axis=1)[:,None])

        self.normals=self.sun_vec+tower_vec
        self.normals/=N.sqrt(N.sum(self.normals**2, axis=1)[:,None])        


    def get_cosine(self):

        cos_factor=N.sum(self.normals*self.sun_vec, axis=1)

        return cos_factor

    def plot_cosine(self):
        x=self.hstpos[:,0]
        y=self.hstpos[:,1]
        z=self.hstpos[:,2]

        plt.figure(1)
        cm = plt.cm.get_cmap('rainbow')
        cs=plt.scatter(x, y, c=self.cosine_factor, cmap=cm,s=60)
        plt.colorbar(cs)
        plt.show()
        
        

        
        
          

if __name__=='__main__':
    #pos_and_aiming=N.loadtxt('/media/yewang/Work/svn_ye/Solstice-tutorial/cases/2-Validation/layout.csv', delimiter=',', skiprows=2)
    pos_and_aiming=N.loadtxt('/home/yewang/Github_repo/anustg/solstice-scripts/src-Linux/srcPy/layout/pos_and_aiming.csv', delimiter=',', skiprows=2)
    pos=pos_and_aiming[:,:3]
    aim=pos_and_aiming[:,4:]
    azimuth=0.
    zenith=12.
    field=Field(pos, aim, azimuth, zenith)
    field.plot_cosine()
        
