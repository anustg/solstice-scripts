import numpy as N
#from tracer.models.heliostat_field import solar_vector
import matplotlib.pyplot as plt

class FieldPF:

    def __init__(self, azimuth, zenith, receiver_norm):
        '''
        Evaluation of field initial performance
        1. cosine factors: sun - heliostat
        2. another cosine factors: receiver -heliostat
        
        azimuth: float, the solar azimuth angle, from South to West
        zenith: float, the solar zenith angle, 0 from vertical
        receiver_normal: (3,) array, normal vector of the receiver aperature 

        '''
        #self.sun_vec=self.get_solar_vector(azimuth, zenith)
        self.rec_norm=receiver_norm.reshape(3,1)
        self.azimuth=azimuth
        self.zenith=zenith


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

    def get_normals(self, towerheight, hstpos, sun_vec):
        
        tower_vec=-hstpos
        tower_vec[:,-1]+=towerheight
        tower_vec/=N.sqrt(N.sum(tower_vec**2, axis=1)[:,None])

        hst_norms=sun_vec+tower_vec
        hst_norms/=N.sqrt(N.sum(hst_norms**2, axis=1)[:,None])        
        return hst_norms

    def get_rec_view(self, towerheight, hstpos):
        # angle between normal of the receiver apterture and the RH
        # RH is the vector from the receiver to the heliostat (i.e. -tower_vec)

        tower_vec=-hstpos
        tower_vec[:,-1]+=towerheight
        tower_vec/=N.sqrt(N.sum(tower_vec**2, axis=1)[:,None])

        view=N.arccos(N.dot(-tower_vec, self.rec_norm))
        view=view.flatten()
        vis_idx=(view<1.5) # the heliostat that can be seen by the receiver

        return vis_idx


    def get_cosine(self,towerheight, hstpos):
        cos_factor=N.zeros(len(hstpos))
        i=0
        for az in self.azimuth:
            for zen in self.zenith:
                sun_vec=self.get_solar_vector(az, zen)
                
                hst_norms=self.get_normals(towerheight, hstpos, sun_vec)
                cos_factor=(cos_factor*float(i)+N.sum(hst_norms*sun_vec, axis=1))/float(i+1)

        return cos_factor


    def plot_cosine(self, savename,pm):
        x=self.hstpos[:,0]
        y=self.hstpos[:,1]
        z=self.hstpos[:,2]
        av=(self.view<N.pi/2.)
        plt.figure(1)
        cm = plt.cm.get_cmap('rainbow')
        cs=plt.scatter(x[av], y[av], c=self.cosine_factor[av], cmap=cm,s=30)
        plt.colorbar(cs)
        plt.title('Cosine factors\n Tower height %s m'%pm )
        plt.savefig(open(savename, 'w'),dpi=500, bbox_inches='tight')
        plt.close()	


    def plot_select(self, savename,pm):
        x=self.hstpos[:,0]
        y=self.hstpos[:,1]
        z=self.hstpos[:,2]

        plt.figure(1)
        fig,ax1=plt.subplots()
        cm = plt.cm.get_cmap('rainbow')
     
        av=(self.view<1.5)
        cs=plt.scatter(x[av], y[av], c=self.view[av], cmap=cm)
        plt.colorbar(cs)
        nv=(self.view>=1.5)
        plt.scatter(x[nv], y[nv], c='gray')
        plt.title('Receiver view \n tilt %s deg'%pm )
        #plt.title('Receiver view \n Tower height %s m'%pm )
        plt.savefig(open(savename, 'w'),dpi=500, bbox_inches='tight')
        plt.close()	


    
            
        
        

        
        
          

if __name__=='__main__':
    #pos_and_aiming=N.loadtxt('/media/yewang/Work/svn_ye/Solstice-tutorial/cases/2-Validation/layout.csv', delimiter=',', skiprows=2)
    pos_and_aiming=N.loadtxt('/home/yewang/Github_repo/anustg/solstice-scripts/src-Linux/srcPy/layout/pos_and_aiming.csv', delimiter=',', skiprows=2)
    pos=pos_and_aiming[:,:3]
    aim=pos_and_aiming[:,4:]
    azimuth=0.
    zenith=12.
    field=Field(pos, aim, azimuth, zenith)
    field.plot_cosine()
        
