from __future__ import print_function
import numpy as np


class MultiApertureConfiguration:

	def __init__(self, n, gamma, H_tower, R_tower, W_rcv, H_rcv, parallel=False):
		"""
		The functions in this class configurate a multi-aperture receiver model, including
			(1) angular position of each aperture
			(2) elevation position and level index of each aperture
			(3) index of the paired aperture
			(4) normal vector of each aperture plane

		by given 
			(1) total number of apertures (n)
			(2) total angular range (gamma)
			(3)	tower height (H_tower)
			(4)	tower radius (R_tower)
			(5) width and height of each aperture (W_rcv, H_rcv)
		    (6) if the multiple apertures are in parallel or cascaded (default is cascaded)

		Note:
			* the index of each aperture (i) starts from the most right aperture and increases counter clockwise
			* the index of the elevation level (lv) starts from the highest and increases downwards 


		Arguments:
			n      :   int, number of apertures (n>=2)
			gamma  : float, the total angular range (the angle btw the most right and left aperture) in deg
			H_tower: float, tower height (m)
			R_tower: float, radius of the tower (m)
			w_rcv  : a list of float that contains the width of each aperture (m)
			h_rcv  : a list of float that contains the height of each aperture (m)
			parallel: bool, the multiple apertures are in parallel or cascaded
		"""
		self.n=n
		self.gamma=gamma
		self.H_tower=H_tower
		self.parallel=parallel

		if n%2==0:
			self.mode='even'
		else:
			self.mode='odd'

		self.elevation_level()
		self.angular_pos()

		self.W_rcv=[]
		self.H_rcv=[]
		if 'lv_1' in W_rcv:	
			# the receiver dimension is given in the sequence of receiver level
			for i in range(self.n):
				lv=self.get_lv_index(i)
				self.W_rcv.append(W_rcv['lv_%.0f'%lv])
				self.H_rcv.append(H_rcv['lv_%.0f'%lv])

		else:
			self.W_rcv=W_rcv
			self.H_rcv=H_rcv


		self.elevation_height()

		W=max(self.W_rcv)*1.2 # 20% space
		alpha=gamma/float(n-1)*np.pi/180.

		if np.tan(alpha/2.)<1e-20:
			self.r=W/2.
		else:
			self.r=W/2./np.tan(alpha/2.)
		if self.r<R_tower:
			self.r=R_tower

		print('Aperture radial distance:', self.r)

	def angular_pos(self):
		self.Omega=np.array([])

		if self.gamma%360.==0:
			for i in range(self.n):
				omega=360./float(self.n)*float(i)-90.				
				self.Omega=np.append(self.Omega, omega)
		else:
			for i in range(self.n):
				omega=90.-self.gamma/2.+self.gamma/float(self.n-1)*float(i)
				self.Omega=np.append(self.Omega, omega)

	def elevation_level(self):
		self.LV=np.array([])

		if self.mode=='odd':
			for i in range(self.n):
				if i<=(self.n-1)/2:
					lv=2*i+1
				else:
					lv=-2*(i-self.n)
				self.LV=np.append(self.LV, lv)

		elif self.mode=='even':
			for i in range(self.n):		
				if i<=(self.n/2-1):
					lv=2*i+1
				else:
					lv=-2*(i-self.n)
				self.LV=np.append(self.LV, lv)

	def elevation_height(self):
		self.Elev=np.array([])
		
		self.i_idx=self.LV.argsort()
		
		for lv in range(1, self.n+1):
			if lv==1:
				elev=self.H_tower
			else:
				idx=self.i_idx[lv-1]
				idx1=self.i_idx[lv-2]

				h_lv=self.H_rcv[idx]
				h_lv1=self.H_rcv[idx1]
				elev-=(h_lv+h_lv1)/2.
			self.Elev=np.append(self.Elev, elev)
			

	def get_angular_pos(self, i):
		return self.Omega[i]
	
	def get_elev_height(self, lv):
		return self.Elev[lv-1]

	def get_lv_index(self, i):
		"return the level index of aperture i"
		return int(self.LV[i])

	def get_i_index(self, lv):
		"return the aperture index at level lv"
		return int(self.i_idx[lv-1])

	def get_pair_idx(self, i):
		if self.mode=='even':
			assert i<=(n/2-1), "i must lower than %s"%(n/2-1)
		elif self.mode=='odd':
			assert i<=((n-1)/2), "i must lower than %s"%((n-1)/2)
		return -i-1

	def get_cood_pos(self, i):
		if self.parallel:
			z_i=self.H_tower

		else:
		# cascaded
			lv=self.get_lv_index(i)
			z_i=self.get_elev_height(lv)

		omega_i=self.get_angular_pos(i)
		x_i=self.r*np.cos(omega_i*np.pi/180.)
		y_i=self.r*np.sin(omega_i*np.pi/180.)

		return x_i, y_i, z_i





