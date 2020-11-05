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

		if perf:
			# the initial large field
			plt.scatter(self.pos_and_aim[:,0], self.pos_and_aim[:,1], c=self.annual)
		else:
			plt.plot(self.pos_and_aim[:,0], self.pos_and_aim[:,1], '.')

		plt.show()
		plt.close()

	def plot_designed_layout(self, perf=True):
		'''
		casedir: str, the case directory
		'''
		x,y,annual=self.layout()

		if perf:
			# the design
			plt.scatter(x, y,c=annual)

		else:
			plt.plot(x, y, '.')

		plt.show()
		plt.close()

if __name__=='__main__':
	casedir='../tests/test-crs-design'
	Case=Case(casedir)
	Case.plot_initial_layout()
	Case.plot_designed_layout()


	
