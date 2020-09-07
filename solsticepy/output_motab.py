import numpy as np
from datetime import datetime

def output_motab(table,savedir=None, title=None):
	'''
	output the .motab table fiel
	'''
	f=open(savedir, 'w')
	f.write('#1\n')

	if table.ndim==2:
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

	else:
		# 3D table, include the breakdown of the total energy
		a=len(table)
		m=np.shape(table[0])[0]-2
		n=np.shape(table[0])[1]-2

		hour_angle=table[0][2, 3:]
		declination=table[0][3:,2]

		for t in range(a):
			f.write('double %s(%s, %s)\n'%(title[t], m,n))

			for i in range(m):
				if i ==0:
					row_i=np.append(0, hour_angle)
				else:
					row_i=np.append(declination[i-1], table[t][2+i, 3:])
				f.write(" ".join(map(str, row_i)))
				f.write("\n")	
						
			f.write("\n")
			f.write("\n")

	f.close()


def output_matadata_motab(table, field_type, aiming, n_helios, A_helio, eff_design, H_rcv, W_rcv, H_tower, Q_in_rcv, A_land, savedir=None, details_en=None):
	"""Output the .motab file to work with the SolarTherm program

	``Arguments``
		* table (numpy array): the oelt that returned by design_crs.CRS.field_design_annual, listed by declination and hour angles
		* field_type (str): polar or surrounding
		* aiming (str): 'single' or 'isp' or others to specify the aiming strategy
		* n_helios (int): total number of heliostats
		* A_helio (float): area of each heliostat (m2)
		* eff_design (float): the optical efficiency at design point
		* H_rcv (float): height of the heliostat (m)
		* W_rcv (float): width of the heliostat (m)
		* H_tower (float), tower height (m)
		* Q_in_rcv (float): the incident power on the receiver (W)
		* A_land (float): the total land area of the field (m2)
		* savedir (str): the directory to save this .motab file
		* details_en (dict): the key of the dict is the name of the breakdown of energy (loss), the value of the dict is a numpy array that contains the value of the energy (loss) with the corresponding declination and hour angles
	
	``Return``
		write the table(s) to the .motab file
	"""
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

	if details_en!=None:
		for key in details_en:
			breakdown=key
			table=details_en[key]
			# size of the lookup table  
			m=np.shape(table)[0]-2
			n=np.shape(table)[1]-2
			f.write('double %s(%s, %s)\n'%(key,m,n))
			for i in range(m):
				if i ==0:
					row_i=np.append(0, hour_angle)
				else:
					row_i=np.append(declination[i-1], table[2+i, 3:])
			f.write(" ".join(map(str, row_i)))
			f.write("\n")				
	f.close()


def read_motab(filename):

	with open(filename) as f:
		content=f.read().splitlines()
	f.close()
	res=content[4].split(',')

	n_helios=float(res[1])
	A_helio=float(res[2])
	eff_des=float(res[3])
	Q_in_rcv=float(res[-2])
	A_land=float(res[-1])

	oelt=np.array([])
	solar_hour=np.array([])
	declination=np.array([])
	t=content[6].split(' ')

	for v in t[1:]:
		solar_hour=np.append(solar_hour, float(v))

	for t in content[7:]:
		v=t.split(' ')
		declination=np.append(declination, float(v[0]))
		oelt=np.append(oelt, np.array(v[1:],dtype=float))

	oelt=oelt.reshape(len(declination), len(solar_hour))

	return n_helios, A_helio, eff_des, Q_in_rcv, A_land, solar_hour, declination, oelt
	


