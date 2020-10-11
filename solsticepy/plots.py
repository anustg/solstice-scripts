import numpy as np
import matplotlib.pyplot as plt



def get_layout(casedir):
	'''
	hst_fn: heliostat pos_and_aim.csv file
	pfm_fn: heliostat performance annual_hst.csv file
	idx_fn: the index of the selected heliostats, selected_hst.csv file
	'''
	hst_fn=casedir+'/des_point/pos_and_aiming.csv'
	pfm_fn=casedir+'/annual_hst.csv'
	idx_fn=casedir+'/selected_hst.csv'

	pos_and_aim=np.loadtxt(hst_fn, skiprows=2, delimiter=',')
	X=pos_and_aim[:,0]
	Y=pos_and_aim[:,1]
	
	annual=np.loadtxt(pfm_fn, delimiter=',')
	selected=np.loadtxt(idx_fn, delimiter=',')
	selected=selected.astype(int)

	return X[selected], Y[selected],annual[selected]


def plot_layout(casedir):
	'''
	hst_fn: heliostat pos_and_aim.csv file
	pfm_fn: heliostat performance annual_hst.csv file
	idx_fn: the index of the selected heliostats, selected_hst.csv file
	'''
	hst_fn=casedir+'/des_point/pos_and_aiming.csv'
	pfm_fn=casedir+'/annual_hst.csv'
	idx_fn=casedir+'/selected_hst.csv'

	pos_and_aim=np.loadtxt(hst_fn, skiprows=2, delimiter=',')
	X=pos_and_aim[:,0]
	Y=pos_and_aim[:,1]
	ROW=pos_and_aim[:,-1]
	
	annual=np.loadtxt(pfm_fn, delimiter=',')
	selected=np.loadtxt(idx_fn, delimiter=',')
	selected=selected.astype(int)

	plt.scatter(X[selected], Y[selected],c=annual[selected])
	#plt.plot(X[selected], Y[selected],'.')
	plt.show()
	plt.close()

if __name__=='__main__':
	casedir='../tests/test-crs-design-tiny-heliostats-philipe2'

	plot_layout(casedir)
	
