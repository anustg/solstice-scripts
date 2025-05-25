import numpy as np

def get_rotation(O, n0):
    '''
    O the center of the facet
    n0 the facet pointing direction
    '''

    n0=n0/np.linalg.norm(n0)

    xn=n0[0]
    yn=n0[1]
    zn=n0[2]

    cosx=zn/np.sqrt(yn**2+zn**2)
    theta_x=np.arccos(cosx)*180./np.pi

    cosy=np.sqrt(yn**2+zn**2)/np.sqrt(xn**2+yn**2+zn**2)
    theta_y=np.arccos(cosy)*180./np.pi

    x=O[0]
    z=O[2]

    rotx=-theta_x
    if x>=0:
        roty=-theta_y
    elif x<0:
        roty=theta_y

    return rotx, roty
  

def find_centroid(xx, yy, weights):
	
	sw=np.sum(weights)
	x_sum=np.sum(xx*weights)
	y_sum=np.sum(yy*weights)
	
	centroid_x=x_sum/sw
	centroid_y=y_sum/sw

	return centroid_x, centroid_y


