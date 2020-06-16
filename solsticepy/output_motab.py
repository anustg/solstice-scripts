import numpy as np
from datetime import datetime

def output_motab(table,savedir=None):
	'''
	output the .motab table fiel
	'''
	f=open(savedir, 'w')
	f.write('#1\n')


	# size of the lookup table
	m=np.shape(table)[0]-2
	n=np.shape(table)[1]-2
	f.write('double optics(%s, %s)\n'%(m,n))

	hour_angle=table[2, 3:]
	declination=table[3:,2]

	for i in range(m):
		if i ==0:
			row_i=np.append(0, hour_angle)
		else:
			row_i=np.append(declination[i-1], table[2+i, 3:])

		#content=np.array2string(row_i, formatter={'float_kind':lambda x: "%.2f" % row_i})
		#content=np.array2string(row_i, precision=2, separator=' ', suppress_small=True)
		#f.write(content+'\n')
		f.write(" ".join(map(str, row_i)))
		f.write("\n")

	f.close()


def output_matadata_motab(table, field_type, aiming, n_helios, A_helio, eff_design, H_rcv, W_rcv, H_tower, Q_in_rcv, A_land, savedir=None):
	'''
	output the .motab table fiel
	'''
	f=open(savedir, 'w')
	f.write('#1\n')
	f.write('#Comments: Field type: %s, Aiming Strategy: %s, Date:%s\n'%(field_type, aiming, datetime.now()))
	f.write('#METALABELS,n_helios,A_helio,Eff_design,H_rcv,W_rcv,H_tower, Q_in_rcv, A_land\n')
	f.write('##METAUNITS,real,m2,real,m,m,m,W,m2\n')
	f.write('#METADATA,%s,%s,%s,%s,%s,%s,%s,%s\n'%(n_helios,A_helio,eff_design,H_rcv,W_rcv,H_tower,Q_in_rcv,A_land))

	# size of the lookup table
	m=np.shape(table)[0]-2
	n=np.shape(table)[1]-2
	f.write('double optics(%s, %s)\n'%(m,n))

	hour_angle=table[2, 3:]
	declination=table[3:,2]

	for i in range(m):
		if i ==0:
			row_i=np.append(0, hour_angle)
		else:
			row_i=np.append(declination[i-1], table[2+i, 3:])

		#content=np.array2string(row_i, formatter={'float_kind':lambda x: "%.2f" % row_i})
		#content=np.array2string(row_i, precision=2, separator=' ', suppress_small=True)
		#f.write(content+'\n')
		f.write(" ".join(map(str, row_i)))
		f.write("\n")

	f.close()
