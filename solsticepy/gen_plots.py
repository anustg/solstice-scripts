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
		


if __name__=='__main__':
	casedir='../tests/test-multi-aperture'
	Case=Case(casedir)
	Case.plot_initial_layout()
	Case.plot_designed_layout()
	Case.plot_aiming()
	Case.plot_oelt(num_aperture=3, savefig=True)

	
