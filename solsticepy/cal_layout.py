import numpy as np
import sys
import os
import matplotlib
#matplotlib.use("agg")
import matplotlib.pyplot as plt
from .cal_sun import *
from .gen_vtk import gen_vtk

def radial_stagger(latitude, num_hst, width, height, hst_z, towerheight, R1, fb, dsep=0., field='polar', num_aperture=0, gamma=0., rec_w=0., rec_z=[], savedir='.', verbose=False, plot=False, plt_aiming=None):
	'''Generate a radial-stagger heliostat field, ref. Collado and Guallar, 2012, Campo: Generation of regular heliostat field.

	``Arguments``
	  * latitude (float): latitude of the field location (deg)
	  * num_hst (int)   : number of heliostats
	  * width (float)   : mirror width (m)
	  * height (float)  : mirror height (m)
	  * hst_z (float)   : the vertical location of each heliostat (m)
	  * towerheight (float): tower height (m)
	  * R1 (float)      : distance from the first row to the bottom of the tower, i.e. (0, 0, 0)
	  * fb (float)      : the field layout growing factor, in (0, 1)
	  * dsep (float)    : separation distance (m)
	  * field (str)     : 'polar-half' or 'surround-half' or 'polar' or 'surround' field or 'multi-aperture', the 'half' option is for simulation a symmetric field
	  * num_aperture(int): number of apertures, for a multi-aperture configuration
	  * gamma (float)   : the anangular range of the multi-aperture configration (deg)	 
	  * rec_z (list)    : a list of the elevation heights of the apertures 
	  * savedir (str)   : directory of saving the pos_and_aiming.csv
	  * verbose(bool)   : write results to disk or not
	  * plot (bool)     : True - plot the layout by Matplotlib

	``Returns``

	  * pos_and_aiming (nx7 numpy array): position, focal length and aiming point of each generated heliostat
	  * a pos_and_aiming.csv file is created and written to the savedir


	``Example``

		>>> from solsticepy.cal_layout import *
		>>> latitude=34.
		>>> hst_width=13.
		>>> hst_height=10.
		>>> hst_z=5.
		>>> target_area=75715.
		>>> num_hst=target_area/hst_width/hst_height*2 # create a twice large field, so that the inefficient heliostat can be trimmed-off after optical simulations
		>>> tower_height=110.
		>>> R1=40. # the distance from the first row to the bottom of the tower
		>>> fb=0.6
		>>> casefolder='.' # replace it as the directory of your case folder
		>>> plot=False # if you want to plot the field layout, set it as True
		>>> pos_and_aim=radial_stagger(latitude, num_hst, hst_width, hst_height, hst_z, tower_height, R1, fb, savedir=casefolder,plot=plot)

		An nx7 array is returned, and the pos_and_aim.csv file is saved in the local directory

	'''

	# heliostat diagonal distantce
	DH=np.sqrt(height**2+width**2) 

	# distance between contiguous helistat center on the X and Y plane
	DM=DH+dsep

	# minimum radial increment
	delta_Rmin=0.866*DM

	# number of heliostats in the first row
	Nhel1 =int(2.*np.pi*R1/DM)


	# the total number of zones (estimated)
	#Nzones=int(np.log(5.44*3*(num_hst/az_rim*np.pi)/Nhel1**2+1)/np.log(4))+1

	X={}
	Y={}
	Nrows_zone=np.array([])
	Nhel_zone=np.array([])
	delta_az_zone=np.array([])

	num=0
	i=0
	sys.stderr.write('DM '+repr(DM)+'\n')
	sys.stderr.write('dRm'+repr(delta_Rmin)+'\n')

	while num<num_hst*3:
		Nrows= int((2.**(i))*Nhel1/5.44)
		Nhel=(2**(i))*Nhel1
		R=Nhel/2./np.pi*DM
		delta_az=2.*np.pi/Nhel

		Nrows_zone=np.append(Nrows_zone, Nrows)
		Nhel_zone=np.append(Nhel_zone, Nhel)
		delta_az_zone=np.append(delta_az_zone, delta_az)

		nh=np.arange(Nhel)
		azimuth=np.zeros((int(Nrows), int(Nhel)))
		azimuth[0::2, :]=delta_az/2.+nh*delta_az # the odd rows
		azimuth[1::2, :]=nh*delta_az

		row=np.arange(Nrows)
		r=R+row*delta_Rmin

		xx=r[:, None]*np.sin(azimuth)
		yy=r[:, None]*np.cos(azimuth)

		X[i]=xx
		Y[i]=yy
		num+=len(xx.flatten())
		print('Zone', i, 'Nrow', Nrows, 'Nhel', Nhel)
		i+=1
	Nzones=i

	sys.stderr.write("\n")
	sys.stderr.write("Denest field %d\n"%(num))

	# expanding the field
	#
	wr=width/height
	const=(1.-(1.-fb)*wr/(2.*wr-(np.sqrt(1.+wr**2)+dsep/height)))*height

	XX=np.array([])
	YY=np.array([])
	ZONE=np.array([])  # zone index
	ROW=np.array([])   # row index among the rows in a zone
	TTROW=np.array([]) # row index among the total rows
	NHEL=np.array([])  # No. index among the heliostats in a row
	AZIMUTH=np.array([])

	for i in range(Nzones):
		Nrows=int(Nrows_zone[i])
		Nhel=int(Nhel_zone[i])
		delta_az=delta_az_zone[i]

		R=np.zeros((Nrows, Nhel))

		if i==0:
			# first zone
			R[0]=R1 # first row
		else:
			# second zones
			R[0,  ::2]=Rn+1.5*DRn
			R[0,  1::2]=Rn+1.5*DRn            
			#R[0,-1]=0.5*(R[0,0]+Rn[-1])

		xx=X[i]
		yy=Y[i]
		zz=np.ones(np.shape(xx))*hst_z
		cosw, coseT=cal_cosw_coset(latitude, towerheight,xx, yy, zz)
		row=np.arange(Nrows)
		cosw=cosw.reshape(Nrows, Nhel)  
		coseT=coseT.reshape(Nrows, Nhel)

		Delta_R=cosw/coseT*const
		Delta_R[Delta_R<delta_Rmin]=delta_Rmin

		for j in range(1, Nrows):
			R[j]=R[j-1]+Delta_R[j-1]

		Rn=R[-1]
		DRn=Delta_R[-1]

		nh=np.arange(Nhel)
		azimuth=np.zeros((Nrows, Nhel))
		azimuth[0::2, :]=delta_az/2.+nh*delta_az # the odd rows
		azimuth[1::2, :]=nh*delta_az

		azimuth=azimuth.flatten()
		R=R.flatten()
		nhels, rows=np.meshgrid(nh, row)
		nhels=nhels.flatten()
		rows=rows.flatten()

		if field=='polar':
			if i<2:
				idx=(azimuth>1.5*np.pi)+(azimuth<0.5*np.pi)

			else:
				idx=(azimuth>(1.5*np.pi+i*np.pi/40.))+(azimuth<(np.pi/2.-i*np.pi/40.))

			xx=R[idx]*np.sin(azimuth[idx])
			yy=R[idx]*np.cos(azimuth[idx])
			AZIMUTH=np.append(AZIMUTH, azimuth[idx])
			rows=rows[idx]
			ROW=np.append(ROW, rows)
			NHEL=np.append(NHEL, nhels[idx])
			zone=np.ones(np.shape(rows))*i

		else:                       
			xx=R*np.sin(azimuth)
			yy=R*np.cos(azimuth)  
			AZIMUTH=np.append(AZIMUTH, azimuth)
			ROW=np.append(ROW, rows)
			NHEL=np.append(NHEL, nhels)
			zone=np.ones(np.shape(rows))*i	

		XX=np.append(XX, xx)
		YY=np.append(YY, yy)
		ZONE=np.append(ZONE, zone)

		if len(TTROW)==0:
			TTROW=np.append(TTROW, rows)
		else:			
			TTROW=np.append(TTROW, rows+np.max(TTROW)+1)
				
	num_hst=int(num_hst)

	if field=='multi-aperture':

		nt=len(XX)
		hstpos=np.zeros(nt*3).reshape(nt, 3)
		hstpos[:, 0]=XX
		hstpos[:, 1]=YY
		hstpos[:,2]=hst_z

		ANGLE=np.array([])
		NORMRCV=np.array([])
		C=np.array([])

		APOS=np.array([])
		for i in range(num_aperture):
			ang_pos, xc, yc=multi_aperture_pos(rec_w, gamma, num_aperture, i)
			print(ang_pos, xc, yc)

			zc=rec_z[i]
			APOS=np.append(APOS, ang_pos)

			c=np.r_[xc, yc, zc]

			oc=np.r_[xc, yc, 0]
			C=np.append(C, c)

			norm_rcv=oc/np.linalg.norm(oc)	
			NORMRCV=np.append(NORMRCV, norm_rcv)	

			vec_CH=hstpos-c	
			L_CH=np.linalg.norm(vec_CH, axis=1)
			L_CH=L_CH.reshape(len(L_CH), 1)
			norm_CH=vec_CH/L_CH
		
			angle=np.arccos(np.sum(norm_rcv*norm_CH, axis=1))
			ANGLE=np.append(ANGLE, angle)

		ANGLE=ANGLE.reshape(num_aperture, len(angle))
		C=C.reshape(num_aperture, 3)

		idx_aim=np.argmin(ANGLE, axis=0)
		angle_min=np.amin(ANGLE, axis=0)

		idx_hst=((angle_min<80.*np.pi/180.)) # valid heliostats that can be seen from the receiver

		idx_aim=idx_aim[idx_hst]
		XX=XX[idx_hst]
		YY=YY[idx_hst]
		ZONE=ZONE[idx_hst]
		ROW=ROW[idx_hst]
		NHEL=NHEL[idx_hst]
		TTROW=TTROW[idx_hst]
		AZIMUTH=AZIMUTH[idx_hst]

		aiming=C[idx_aim]
		aim_x=aiming[:,0]
		aim_y=aiming[:,1]	

		aim_x=aim_x[:num_hst]
		aim_y=aim_y[:num_hst]	
		idx_aim=idx_aim[:num_hst]

		aim_z=np.array([])
		for i in range(num_hst):
			aim_z=np.append(aim_z, rec_z[idx_aim[i]])


	else:
		aim_x=np.zeros(num_hst)
		aim_y=np.zeros(num_hst)
		aim_z=np.ones(num_hst)*towerheight
		idx_aim=np.zeros(num_hst)


	XX=XX[:num_hst]
	YY=YY[:num_hst]
	ZONE=ZONE[:num_hst]
	ROW=ROW[:num_hst]
	NHEL=NHEL[:num_hst]
	TTROW=TTROW[:num_hst]
	AZIMUTH=AZIMUTH[:num_hst]*180./np.pi


	sys.stderr.write("\nExpanded field %d\n"%(num_hst,))

	hstpos=np.zeros(num_hst*3).reshape(num_hst, 3)
	hstpos[:, 0]=XX
	hstpos[:, 1]=YY
	hstpos[:,2]=hst_z

	# the aiming point is a default point
	# it is required to be revised for receiver performance


	foc=np.sqrt((XX-aim_x)**2+(YY-aim_y)**2+(hstpos[:,2]-aim_z)**2)

	pos_and_aiming=np.append(XX, (YY, hstpos[:,2], foc, aim_x, aim_y, aim_z, idx_aim, AZIMUTH, ZONE, ROW, NHEL, TTROW, np.arange(num_hst)))
	title=np.array(['x', 'y', 'z', 'foc', 'aim x', 'aim y', 'aim z', 'aim-rec-index','Azimuth pos','Zone', 'Row', 'No.', 'row index', 'No. index',  'm', 'm', 'm', 'm', 'm', 'm', 'm', '-', 'deg','-', '-', '-', '-', '-'])
	pos_and_aiming=pos_and_aiming.reshape(14, num_hst)
	pos_and_aiming=np.append(title, pos_and_aiming.T)
	pos_and_aiming=pos_and_aiming.reshape(num_hst+2, 14)

	if verbose:
		if not os.path.exists(savedir):
			os.makedirs(savedir)
		np.savetxt('%s/pos_and_aiming.csv'%savedir, pos_and_aiming, fmt='%s', delimiter=',')

	if plot:
		if not os.path.exists(savedir):
			os.makedirs(savedir)

		fts=24
		plt.figure(dpi=100.,figsize=(12,9))
		plt.plot(XX, YY, '.')
		#plt.xlim(-1000, 1000)
		#plt.ylim(-1000, 1000)
		plt.xticks(fontsize=fts)
		plt.yticks(fontsize=fts)
		plt.xlabel('x (m)', fontsize=fts)
		plt.ylabel('y (m)', fontsize=fts)
		plt.savefig(savedir+'/field_layout.png', bbox_inches='tight')
		plt.close()

	'''
	if plt_aiming!=None:

		NORMRCV=NORMRCV.reshape(num_aperture, 3)
		plt.figure(dpi=100.,figsize=(12,9))
		plt.scatter(XX, YY, c=np.arctan(aim_x/aim_y))
		

		plt.scatter(C[:,0], C[:,1], s=1)

		origin = np.array([[0, 0, 0],[0, 0, 0]]) 	
		plt.quiver(*origin, NORMRCV[:,0], NORMRCV[:,1],scale=10)
		#plt.colorbar()
		#plt.grid()
		plt.savefig(savedir+'/aiming_%s.png'%plt_aiming, bbox_inches='tight')
		plt.close()
	'''
	return pos_and_aiming, Nzones, Nrows_zone

def cal_cosw_coset(latitude, towerheight, xx, yy, zz):
	'''
	The factors to growing the heliostat field, see eq.(2) Francisco J. Collado, Jesus Guallar, Campo: Generation of regular heliostat fields, 2012

	``Arguments``
	  * latitude (float)   : latitude of the field location (deg)
	  *	towerheight (float): tower height (m)
      * xx, yy, zz (float) : coordinates of heliostats

	``Returns``
	  * cosw (array) : cos(omega)
	  * coseT (array): cos(epsilon_T) 

	'''

	hst_pos=np.append(xx, (yy, zz))
	hst_pos=hst_pos.reshape(3, len(xx.flatten())) # 3 x n
	tower_vec=-hst_pos
	tower_vec[-1]+=towerheight
	tower_vec/=np.sqrt(np.sum(tower_vec**2, axis=0)) # 3 x n
	unit=np.array([[0.], [0.], [1.]])
	unit=np.repeat(unit, len(tower_vec[0]), axis=1)
	coseT=np.sum(tower_vec*unit, axis=0)

	sun=SunPosition()
	dd=sun.days(21, 'Mar')
	delta=sun.declination(dd)
	h=np.arange(8, 17)

	omega= -180.+15.*h
	theta=sun.zenith(latitude, delta, omega) # solar zenith angle 

	phi=np.array([]) # solar azimuth angle
	for i in range(len(h)):
		p=sun.azimuth(latitude, theta[i], delta, omega[i])
		phi=np.append(phi, p)

	theta*=np.pi/180.
	phi*=np.pi/180.
	 
	cosw=np.zeros(len(tower_vec[0]))
	sun_z = np.cos(theta)
	sun_y=-np.sin(theta)*np.cos(phi)
	sun_x=-np.sin(theta)*np.sin(phi)
	sun_vec = np.append(sun_x, (sun_y, sun_z))
	sun_vec=sun_vec.reshape(3, len(sun_x)) # 3xs

	for i in range(len(sun_vec[0])):
		sv=np.repeat(sun_vec[:,i], len(tower_vec[0])).reshape(3, len(tower_vec[0]))
		hst_norms=sv+tower_vec
		hst_norms/=np.linalg.norm(hst_norms, axis=0)
		cosw+=np.sum(sv*hst_norms, axis=0)
	cosw/=float(len(sun_vec[0]))
	return cosw, coseT


def aiming_cylinder(r_height,r_diameter, pos_and_aiming, savefolder, c_aiming=0.):
	'''
	The aiming method is developed following the deviation-based multiple aiming, by Shuang Wang. Reference: Augsburger G. Thermo-economic optimisation of large solar tower power plants[R]. EPFL, 2013.

	``Arguments``
	  * r_height (float)   : receiver height (m)
	  * r_diameter (float) : receiver diameter (m)
	  * pos_and_aiming (array): the array returned from the function radial_stagger
	  * savefolder (str)   : directory to save the aiming point results
	  * c_aiming (float)   : an aiming co-efficient

	``Returns``

	  * hst_info_ranked (array) : the heliostats ranked according to focal lenghts	
	  * a pos_and_aiming.csv file is created and written to the savedir, which contains pos_and_aiming (nx7 numpy array): position, focal length and the updated aiming point of each heliostat 

	'''
	r_radius=0.5*r_diameter
	hst_info=pos_and_aiming[2:].astype(float)
	num_hst=hst_info.size/7
	#print num_hst
	foc=hst_info[:,3] # the focal lenghts of all the heliostats

	# ranked the hsts according to focal lenghts
	hst_info_ranked = hst_info[np.argsort(foc)[::1]]

	for i in range(num_hst):
		
		if (i+1)%2==0: # negative
			hst_info_ranked[i,6]=hst_info_ranked[i,6]+0.5*r_height*c_aiming*(float(num_hst)-1-i)/num_hst
		else:
			hst_info_ranked[i,6]=hst_info_ranked[i,6]-0.5*r_height*c_aiming*(float(num_hst)-1-i)/num_hst
		
		hst_info_ranked[i,4]=hst_info_ranked[i,0]*r_radius/np.sqrt(hst_info_ranked[i,0]**2+hst_info_ranked[i,1]**2)
		hst_info_ranked[i,5]=hst_info_ranked[i,1]*r_radius/np.sqrt(hst_info_ranked[i,0]**2+hst_info_ranked[i,1]**2)
		#print hst_info_ranked[i,0],hst_info_ranked[i,1],hst_info_ranked[i,4],hst_info_ranked[i,5]
		hst_info_ranked[i,3]=np.sqrt((hst_info_ranked[i,0]-hst_info_ranked[i,4])**2+(hst_info_ranked[i,1]-hst_info_ranked[i,5])**2+(hst_info_ranked[i,2]-hst_info_ranked[i,6])**2)


	title=np.array(['x', 'y', 'z', 'foc', 'aim x', 'aim y', 'aim z', 'm', 'm', 'm', 'm', 'm', 'm', 'm'])

	pos_and_aiming_new=np.append(title, hst_info_ranked)
	pos_and_aiming_new=pos_and_aiming_new.reshape(len(pos_and_aiming_new)/7, 7)

	csv_new=savefolder+'/pos_and_aiming.csv'# the output field file
	np.savetxt(csv_new, pos_and_aiming_new, fmt='%s', delimiter=',')

	return hst_info_ranked	

def multi_aperture_pos(rec_w, gamma, n, i):
	"""
	This function returens the angular position of each aperture
	in the multi-aperture configration that has n apertures
	in the angular range of gamma

	Arguments:
	rec_w, list, a list of the width of all the apertures
	gamma, float, angular range (deg) of the multi-aperture configration
			      which is defined as the angle from the most right to 
				  the most left aperture
	n, int, number of apertures
	i, int, the i-th aperture (starts from 0 for the most right aperture)

	Return:
	omega_i, float, the angular position of the i-th aperture
					in the coordinate system, the angular position 
					starts from +x and increases counter-clockwise
	"""	
	if gamma%360.==0:
		omega_i=360./float(n)*float(i)-90.
	else:
		omega_i=90.-gamma/2.+gamma/float(n-1)*float(i)

	W=max(rec_w)*1.2 # 20% space
	alpha=gamma/float(n-1)*np.pi/180.
	if np.tan(alpha/2.)<1e-20:
		r=W/2.
	else:
		r=W/2./np.tan(alpha/2.)

	xc=r*np.cos(omega_i*np.pi/180.)
	yc=r*np.sin(omega_i*np.pi/180.)

	return omega_i, xc, yc


if __name__=='__main__':
    
	pos_and_aim=radial_stagger(latitude=34., num_hst=22640, width=10., height=10., hst_z=0., towerheight=250, R1=80, fb=0., dsep=0., field='polar', savedir='.', plot=True)

        
