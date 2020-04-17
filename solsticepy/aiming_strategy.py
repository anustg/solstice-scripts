'''
Modified deviation-based aiming strategy for flat receiver.
'''

import numpy as np
from sys import path
#from numpy import *

def aiming_flat(folder,h_rec,l_rec,C_aiming,csv,tower_h,Exp):
	"""Aiming strategy for a flat receiver

	TODO add details here.
	"""

	title=np.array(['x', 'y', 'z', 'foc', 'aim x', 'aim y', 'aim z', 'm', 'm', 'm', 'm', 'm', 'm', 'm'])
	hst_info=np.loadtxt(csv,delimiter=',', skiprows=2) 
	num_hst=hst_info.size/7
	
	pos_and_aiming_new=np.array([])
	foc=hst_info[:,3]
	hst_info_ranked = hst_info[np.argsort(foc)[::1]]
	for j in range(len(hst_info)):
		lmax=np.max(hst_info_ranked[:,3])
		lmin=np.min(hst_info_ranked[:,3])
		li=hst_info_ranked[j,3]
		
		# for z coordinates
		d0=0.5*h_rec*C_aiming[0]*((lmax-li)/(lmax-lmin))**Exp[0]
		r0=round(random.uniform(0.,1),3)
		if r0<0.5:
			random_index0=-1
		elif r0 >=0.5:
			random_index0=1
		# the random index is -1 or +1
		hst_info_ranked[j,6]=tower_h+d0*random_index0 # z coordinates of aiming points
		
		# for x coordinates
		d1=0.5*h_rec*C_aiming[1]*((lmax-li)/(lmax-lmin))**Exp[1]
		r1=round(random.uniform(0.,1),3)
		if r1<0.5:
			random_index1=-1
		elif r1 >=0.5:
			random_index1=1
		hst_info_ranked[j,4]=d1*random_index1 # x coordinates of aiming points
		hst_info_ranked[j,3]=np.sqrt((hst_info_ranked[j,0]-hst_info_ranked[j,4])**2+(hst_info_ranked[j,1]-hst_info_ranked[j,5])**2+(hst_info_ranked[j,2]-hst_info_ranked[j,6])**2) # focal lenght
	pos_and_aiming_new=np.append(pos_and_aiming_new, hst_info_ranked)
	
	pos_and_aiming_new=np.append(title,pos_and_aiming_new)
	pos_and_aiming_new=pos_and_aiming_new.reshape(len(pos_and_aiming_new)/7, 7)
	csv_new=csv='%s/pos_and_aiming_new.csv' % folder # the output field file
	np.savetxt(csv_new, pos_and_aiming_new, fmt='%s', delimiter=',')
	
if __name__=='__main__':
	h_rec=2.75 # receiver height
	l_rec=2.75
	folder=path[0]
	csv='%s/sandia_hstat_coordinates_new.csv' % folder # the input field file
	tower_h=60. # that is the optical tower height

	C_aiming=np.array([0.5,0.5]) # in two directions, z and x, i==0 is z; i==1 is x
	Exp=np.array([1.5,1.5])
	aiming_flat(folder,h_rec,l_rec,C_aiming,csv,tower_h,Exp)
