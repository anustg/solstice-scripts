import numpy as np
import matplotlib.pyplot as plt

class Case:

	def __init__(self,casedir):
		'''
		casedir: str, the case directory
		'''
		self.casedir=casedir
		self.layout()

	def layout(self):
		'''
		hst_fn: heliostat pos_and_aim.csv file
		pfm_fn: heliostat performance annual_hst.csv file
		idx_fn: the index of the selected heliostats, selected_hst.csv file
		'''
		hst_fn=self.casedir+'/des_point/pos_and_aiming.csv'
		pfm_fn=self.casedir+'/annual_hst.csv'
		idx_fn=self.casedir+'/selected_hst.csv'

		self.pos_and_aim=np.loadtxt(hst_fn, skiprows=2, delimiter=',')
		X=self.pos_and_aim[:,0]
		Y=self.pos_and_aim[:,1]

		self.annual=np.loadtxt(pfm_fn, delimiter=',')
		selected=np.loadtxt(idx_fn, delimiter=',')
		selected=selected.astype(int)

		return X[selected], Y[selected], self.annual[selected]


	def plot_initial_layout(self, perf=True):
		'''
		casedir: str, the case directory
		'''
		x,y,annual=self.layout()
		if perf:
			# the initial large field
			plt.scatter(self.pos_and_aim[:,0], self.pos_and_aim[:,1], c=self.annual)
			plt.plot(x, y, 'r.')
		else:
			plt.plot(self.pos_and_aim[:,0], self.pos_and_aim[:,1], '.')

		plt.colorbar()
		plt.show()
		plt.close()

	def plot_designed_layout(self, perf=True, vmin=None, vmax=None, markersize=5, savefig=None, title=None, xlabel=None, ylabel=None):
		'''
		casedir: str, the case directory
		'''
		x,y,annual=self.layout()
		print(np.max(annual), np.min(annual))

		if perf:
			# the design
			if vmin==None:
				plt.scatter(x, y, c=annual, s=markersize)
			else:
				plt.scatter(x, y, c=annual, s=markersize, vmin=vmin, vmax=vmax)

		else:
			plt.plot(x, y, '.')
		plt.colorbar()
		if title!=None:
			plt.title(title)
		if xlabel!=None:
			plt.xlabel(xlabel)
		if ylabel!=None:
			plt.ylabel(ylabel)
		if savefig==None:
			plt.show()
		else:
			plt.savefig(savefig, bbox_inches='tight')

		plt.close()

	def plot_aiming(self, savefig=None):
		X=self.pos_and_aim[:,0]
		Y=self.pos_and_aim[:,1]
		aim_x=self.pos_and_aim[:,4]
		plt.scatter(X, Y, c=aim_x)	
		if savefig==None:
			plt.show()
		else:
			plt.savefig(savefig, bbox_inches='tight')
		plt.close()	

	def plot_oelt(self, num_aperture=1, vmax=0.85, vmin=0.45,savefig=None):

		if num_aperture==1:
			table=np.loadtxt(self.casedir+'/lookup_table.csv', dtype=str, delimiter=',')

			## comparison
			dec=table[3:,2].astype(float)
			hra=table[2, 3:].astype(float)
			oelt=table[3:,3:].astype(float)

			plt.figure(1)
			plt.pcolormesh(hra, dec, oelt)
			plt.colorbar()
			plt.xlabel('Solar hour angle (deg)')
			plt.ylabel('Declination angle (deg)')
			plt.show()
			plt.close()

		else:
			for i in range(num_aperture):
				table=np.loadtxt(self.casedir+'/lookup_table_%s.csv'%i, dtype=str, delimiter=',')

				## comparison
				dec=table[3:,2].astype(float)
				hra=table[2, 3:].astype(float)
				oelt=table[3:,3:].astype(float)

				plt.figure(1)
				plt.pcolormesh(hra, dec, oelt, vmax=vmax, vmin=vmin)
				plt.colorbar()
				plt.xlabel('Solar hour angle (deg)')
				plt.ylabel('Declination angle (deg)')

				if savefig==None:
					plt.show()
				else:
					plt.savefig(self.casedir+'/oelt_aperture_%s.png'%i, bbox_inches='tight')
				plt.close()

			table=np.loadtxt(self.casedir+'/lookup_table_total.csv', dtype=str, delimiter=',')

			## comparison
			dec=table[3:,2].astype(float)
			hra=table[2, 3:].astype(float)
			oelt=table[3:,3:].astype(float)

			plt.figure(1)
			plt.pcolormesh(hra, dec, oelt, vmax=vmax, vmin=vmin)
			plt.colorbar()
			plt.xlabel('Solar hour angle (deg)')
			plt.ylabel('Declination angle (deg)')
			if savefig==None:
				plt.show()
			else:
				plt.savefig(self.casedir+'/oelt_total.png', bbox_inches='tight')

			plt.close()
		
def plot_fluxmap_cylinder(points, tri, flux, casedir, casename, loc_z_rec=171.035, rec_r=7.75, rec_h= 25.05609, m=31, n=60):

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


	flux=(flux[::2]+flux[1::2])/2.
	flux=flux[:-n]
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


def plot_fluxmap_flat(points, tri, flux, casedir,  m=None, n=None):

	X=points[:,0]
	Y=points[:,1]
	Z=points[:,2]

	print(len(X), len(tri), len(flux))
	plt.tripcolor(X, Y, tri, facecolors=flux, cmap='jet') 
	plt.colorbar()
	plt.savefig(open(casedir+'/flux_tri.png', 'wb'), bbox_inches='tight')
	plt.close() 



if __name__=='__main__':
	casedir='../tests/test-multi-aperture'
	Case=Case(casedir)
	Case.plot_initial_layout()
	Case.plot_designed_layout()
	Case.plot_aiming()
	Case.plot_oelt(num_aperture=3, savefig=True)

	
