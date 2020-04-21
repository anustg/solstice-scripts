import numpy as np
#from tracer.models.heliostat_field import solar_vector
#import matplotlib.pyplot as plt
from scipy.spatial import Delaunay
from .gen_vtk import gen_vtk

class FieldPF:
	"""Preliminary calculation of heliostat field performance.

	1. cosine factors: sun - heliostat
	2. another cosine factors: receiver - heliostat

	Note the angle conventions in this program:

	  * azimuth: the solar azimuth angle, from South to West
	  * zenith: the solar zenith angle, 0 from vertical


	``Example``

		>>> from solsticepy.cal_field import *
		>>> pos_and_aiming=np.loadtxt('./pos_and_aiming.csv', delimiter=',', skiprows=2) #load the field layout file (refer to cal_layout.py)
		>>> pos=pos_and_aiming[:,:3]
		>>> aim=pos_and_aiming[:,4:]
		>>> azimuth=np.r_[0.]
		>>> zenith=np.r_[12.]
		>>> field=FieldPF(np.r_[0,1,0])
		>>> sun_vec=field.get_solar_vector(azimuth, zenith)
		>>> norms=field.get_normals(towerheight=70., hstpos=pos, sun_vec=sun_vec)
		>>> COORD, TRI, ele, nc=field.mesh_heliostat_field(width=10., height=8., normals=norms, hstpos=pos)
		>>> cos=field.get_cosine(hst_norms=norms, sun_vec=sun_vec)
		>>> savedir='./field.vtk'
		>>> COS=np.repeat(cos, ele)
		>>> DATA={'cos':COS}
		>>> NORMS=np.repeat(norms, ele, axis=0)
		>>> gen_vtk(savedir, COORD.T, TRI, NORMS, True, DATA)

		The field.vtk file is saved in the local directory and can be open in ParaView to visualise the cosine factor of each individual heliostat. See more details in the 'gen_vtk.py' in the next section.


	"""

	def __init__(self, receiver_norm=np.r_[0,1,0]):
		"""

		``Argument``

		  * receiver_normal (numpy array): A numpy 3-vector with unit normal direction of the receiver aperature (default: [0,1,0])

		"""
		self.rec_norm=receiver_norm.reshape(3,1)


	def get_solar_vector(self, azimuth, zenith):
		"""Calculate the solar vector using azimuth and zenith angles.

		``Arguments``

		  * azimuth (float): the sun's azimuth (deg), from South increasing towards to the West
		  * zenith (float): angle created between the solar vector and the Z axis (deg)

		``Return``

		  * sun_vec (numpy array): a 3-component 1D array with the solar vector
		"""
		#copied from tracer.models.heliostat_field import solar_vector
		#TODO change it to Solstice convetion if necessary

		azimuth*=np.pi/180.
		zenith*=np.pi/180.

		sun_z = np.cos(zenith)
		sun_y=-np.sin(zenith)*np.cos(azimuth)
		sun_x=-np.sin(zenith)*np.sin(azimuth)
		sun_vec = np.r_[sun_x, sun_y,sun_z] 

		return sun_vec

	def get_normals(self, towerheight, hstpos, sun_vec):
		"""Calculate the normal vectors of each heliostat

		``Arguments``

		  * towerheight (float): tower height
		  * hstpos (nx3 numpy array): heliostat positions
		  * sun_vec (numpy array): the solar vector

		``Return``

		  * hst_norms (numpy array): normal vectors of each heliostat

		"""

		tower_vec=-hstpos
		tower_vec[:,-1]+=towerheight
		tower_vec/=np.sqrt(np.sum(tower_vec**2, axis=1)[:,None])

		hst_norms=sun_vec+tower_vec
		hst_norms/=np.sqrt(np.sum(hst_norms**2, axis=1)[:,None])        
		return hst_norms

	def get_rec_view(self, towerheight, hstpos):
		"""Check the visibility of each heliostat from the view of the receiver
	
		``Arguments``
	
		  * towerheight (float): tower height 
		  * hstpos (numpy array): position of each heliostat

		``Return``

		  * vis_idx (numpy array): the indices of the heliostats that can be seen by the receiver 
		
		"""

		# angle between normal of the receiver apterture and the RH
		# RH is the vector from the receiver to the heliostat (i.e. -tower_vec)

		tower_vec=-hstpos
		tower_vec[:,-1]+=towerheight
		tower_vec/=np.sqrt(np.sum(tower_vec**2, axis=1)[:,None])

		view=np.arccos(np.dot(-tower_vec, self.rec_norm))
		view=view.flatten()
		vis_idx=(view<1.5) # the heliostat that can be seen by the receiver

		return vis_idx


	def get_cosine(self,hst_norms, sun_vec):

		"""Calculate the cosine factor between the sun and the heliostats

		``Arguments``

		  * hst_norms (numpy array): normal vectors of the heliostats
		  * sun_vec (numpy array): solar vector

		``Return``

		  * cos_factor (numpy array): the cosine factor between the sun and the heliostats
		
		"""

		cos_factor=np.sum(hst_norms*sun_vec, axis=1)

		return cos_factor

	def mesh_heliostat(self, width, height):
		"""The local coordinate of the triangular mesh of a heliostat 

		``Arguments``
		  * width (float): width of the heliostat
		  * height (float): height of the heliostat

		``Returns``

		  * coords (numpy array): coordinates of the vertices
		  * tri (numpy array): the indices of the triangular mesh

		"""
		x=np.linspace(-width/2., width/2., 2)
		y=np.linspace(-height/2., height/2., 2)

		xx, yy=np.meshgrid(x, y)
		coords=np.column_stack([xx.ravel(),yy.ravel()])   
		tri=Delaunay(coords).simplices
		#plt.figure(1) 
		#plt.triplot(coords[:,0], coords[:,1], tri)   
		#plt.show()
		return coords, tri

	def mesh_heliostat_field(self, width, height, normals, hstpos):
		"""Generate the necessary elements to create the VTK file to view the heliostat layout 

		``Arguments``
		  * width (float): width of the heliostat
		  * height (float): height of the heliostat
		  * normals (numpy array): normal vectors of the heliostats
		  * hstpos (numpy array): the positions of the heliostats

		``Returns``

		  * COORD (numpy array): the coordinates of the indices of the helisotat field
		  * TRI (numpy array): the indices of the triangular mesh of the heliostat field
		  * ele (int): number of element of the triangular mesh
		  * nc (int): number of indices
		"""		


		coord1, tri1=self.mesh_heliostat(width, height)

		ele=len(tri1)
		nc=len(coord1)

		num_hst=len(hstpos)

		COORD=np.zeros(num_hst*nc*3).reshape(num_hst*nc, 3)
		TRI=np.zeros(num_hst*ele*3).reshape(num_hst*ele, 3)

		norm_x=normals[:,0]
		norm_y=normals[:,1]
		norm_z=normals[:,2]          

		for i in range(num_hst):            
		 
		    TRI[i*ele: (i+1)*ele]=tri1+i*nc

		    trans = np.eye(4)
		    trans[:3,3] = hstpos[i]

		    x=coord1[:,0]
		    y=coord1[:,1]
		    z=float(hstpos[i,2])*np.ones(nc)  

		    cd=np.vstack((x, y, z, np.ones(nc)))   
		    cd_t=np.dot(trans, cd)

		    xx=cd_t[0]
		    yy=cd_t[1]
		    zz=cd_t[2]

		    COORD[i*nc: (i+1)*nc,0]=xx
		    COORD[i*nc: (i+1)*nc,1]=yy
		    COORD[i*nc: (i+1)*nc,2]=zz

		#plt.figure(1) 
		#plt.triplot(COORD[:,0], COORD[:,1], TRI)   
		#plt.show()        
		return COORD, TRI, ele, nc

	def plot_cosine(self, savename):
		"""Plot the cosine factors of heliostats in Matplotlib

		``Argument``

		  * savename (str): the directory to save the figure, with suffix, e.g. '.png' or '.jpg'

		``Return``

		  * No return value (a figure is written in the specified path)

		"""
		x=self.hstpos[:,0]
		y=self.hstpos[:,1]
		z=self.hstpos[:,2]
		av=(self.view<np.pi/2.)
		plt.figure(1)
		cm = plt.cm.get_cmap('rainbow')
		cs=plt.scatter(x[av], y[av], c=self.cosine_factor[av], cmap=cm,s=30)
		plt.colorbar(cs)
		plt.title('Cosine factors')
		plt.savefig(open(savename, 'w'),dpi=500, bbox_inches='tight')
		plt.close()	


	def plot_select(self, savename, tilt):
		"""Plot the heliostats that can be seen by the receiver in Matplotlib

		``Argument``

		  * savename (str): the directory to save the figure, with suffix, e.g. '.png' or '.jpg'
		  * savename (float): receiver tilted angle

		``Return``

		  * No return value (a figure is written in the specified path)

		"""

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
       
def rotx(ang):
    """Generate a homogenous transform for ang radians around the x axis"""
    s = np.sin(ang); c = np.cos(ang)
    return np.array([
        [1., 0, 0, 0],
        [0, c,-s, 0],
        [0, s, c, 0],
        [0, 0, 0, 1.]
    ])

def roty(ang):
    """Generate a homogenous transform for ang radians around the y axis"""
    s = np.sin(ang); c = np.cos(ang)
    return np.array([
        [c, 0, s, 0],
        [0, 1., 0, 0],
        [-s,0, c, 0],
        [0, 0, 0, 1.]
    ])

def rotz(ang):
    """Generate a homogenous trransform for ang radians around the z axis"""
    s = np.sin(ang); c = np.cos(ang)
    return np.array([
        [c,-s, 0, 0],
        [s, c, 0, 0],
        [0, 0, 1., 0],
        [0, 0, 0, 1.]
    ])

def translate(x=0, y=0, z=0):
    """Generate a homogenous transform for translation by x, y, z"""
    return np.array([
        [1., 0, 0, x],
        [0, 1., 0, y],
        [0 ,0, 1., z],
        [0, 0, 0, 1.]
    ])
          

if __name__=='__main__':
    #pos_and_aiming=np.loadtxt('/media/yewang/Work/svn_ye/Solstice-tutorial/cases/2-Validation/layout.csv', delimiter=',', skiprows=2)
    pos_and_aiming=np.loadtxt('./pos_and_aiming.csv', delimiter=',', skiprows=2)
    pos=pos_and_aiming[:,:3]
    aim=pos_and_aiming[:,4:]
    azimuth=np.r_[0.]
    zenith=np.r_[12.]
    field=FieldPF(np.r_[0,1,0])
    sun_vec=field.get_solar_vector(azimuth, zenith)
    norms=field.get_normals(towerheight=70., hstpos=pos, sun_vec=sun_vec)
    #field.heliostat(10, 8)
    COORD, TRI, ele, nc=field.view_heliostats(width=10., height=8., normals=norms, hstpos=pos)
    cos=field.get_cosine(hst_norms=norms, sun_vec=sun_vec)
    savedir='./field.vtk'
    COS=np.repeat(cos, ele)
    DATA={'cos':COS}
    NORMS=np.repeat(norms, ele, axis=0)
    gen_vtk(savedir, COORD.T, TRI, NORMS, True, DATA)
        
