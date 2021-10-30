import numpy as np
from datetime import datetime
import re

def output_motab(table,savedir=None, title=None):
	'''
	output the .motab table file
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


def output_matadata_motab(table, field_type, aiming, n_helios, A_helio, eff_design, eff_annual, H_rcv, W_rcv, H_tower, Q_in_rcv, A_land, savedir=None, details_en=None):
	"""Output the .motab file to work with the SolarTherm program

	``Arguments``
		* table (numpy array): the oelt that returned by design_crs.CRS.field_design_annual, listed by declination and hour angles
		* field_type (str): polar or surrounding
		* aiming (str): 'single' or 'isp' or others to specify the aiming strategy
		* n_helios (int): total number of heliostats
		* A_helio (float): area of each heliostat (m2)
		* eff_design (float): the optical efficiency at design point
		* eff_annual (float): the annual optical efficiency (dni weighted)
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
	f.write('#METALABELS,n_helios,A_helio,Eff_design,Eff_annual,H_rcv,W_rcv,H_tower, Q_in_rcv, A_land\n')
	f.write('##METAUNITS,real,m2,real,real,m,m,m,W,m2\n')
	f.write('#METADATA,%s,%s,%s,%s,%s,%s,%s,%s,%s\n'%(n_helios,A_helio,eff_design,eff_annual,H_rcv,W_rcv,H_tower,Q_in_rcv,A_land))

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



def output_matadata_motab_multi_aperture(TABLE, eff_design, eff_annual, A_land, H_tower, A_helio, n_helios_total, Q_in_rcv_total, num_aperture, Q_in_rcv, n_helios, H_rcv, W_rcv, Z_rcv, mac, savedir=None):
	"""Output the .motab file to work with the SolarTherm program

	``Arguments``
		* TABLE   (dic): each key is the index of each aperture, and the value is the corresponding oelt (numpy array) that returned by design_crs.CRS.field_design_annual, listed by declination and hour angles
		* eff_design    (float): the optical efficiency at design point
		* eff_annual    (float): the annual optical efficiency (dni weighted)
		* A_land        (float): the total land area of the field (m2)
		* H_tower       (float): tower height (m)
		* A_helio       (float): area of each heliostat (m2)
		* n_helios_total  (int): total number of heliostats
		* Q_in_rcv_total(float): the total required incident power on the receiver at design point (W)
		* num_aperture    (int): number of apertures in the multi-aperture configuration
		* Q_in_rcv       (list): a list of the incident power on each of the aperture at design point (W)
		* n_helios       (list): a list of number of heliostats that associated with each aperture
		* H_rcv          (list): a list of height of each aperture (m)
		* W_rcv          (list): a list of width of each aperture (m)
		* Z_rcv          (list): a list of elevation height of each aperture (m)
		* mac                  : MultiApertureConfiguration class instance
		* savedir         (str): the directory to save this .motab file
		* details_en     (dict): the key of the dict is the name of the breakdown of energy (loss), the value of the dict is a numpy array that contains the value of the energy (loss) with the corresponding declination and hour angles
	
	``Return``
		write the table(s) to the .motab file
	"""
	labels='#METALABELS,eff_design,eff_annual,A_land,H_tower,A_helio,n_helios_total,Q_in_rcv_total,num_aperture, radial_distance'
	units='##METAUNITS,real,real,m2,m,m2,real,W,real, m'
	data='#METADATA,%s,%s,%s,%s,%s,%s,%s,%s,%s'%(eff_design, eff_annual, A_land, H_tower, A_helio, n_helios_total, Q_in_rcv_total, num_aperture, mac.r)

	# the sequence is correspondence to the level index, 
	for lv in range(1, num_aperture+1):
		i=mac.get_i_index(lv)

		labels+=',Q_in_rcv_%s,n_helios_%s,H_rcv_%s,W_rcv_%s, Z_rcv_%s'%(i, i, i, i, i)
		units+=',W,real,m,m'
		data+=',%s,%s,%s,%s,%s'%(Q_in_rcv[i], n_helios[i], H_rcv[i], W_rcv[i], Z_rcv[i])

	labels+='\n'
	units+='\n'
	data+='\n'

	f=open(savedir, 'w')
	f.write('#1\n')
	f.write('#Comments: Field type: multi-aperture, Aiming Strategy: to-each-aperture-centre, Date:%s\n'%(datetime.now()))
	f.write(labels)
	f.write(units)
	f.write(data)
	f.write('\n')

	for i in range(num_aperture+1):
		table=TABLE[i]

		# size of the lookup table  
		m=np.shape(table)[0]-2
		n=np.shape(table)[1]-2

		if i==num_aperture:
			f.write('double optical_efficiency_total(%s, %s)\n'%(m,n))
		else:
			lv=mac.get_lv_index(i)
			f.write('double optical_efficiency_level_%.0f(%s, %s)\n'%(lv, m,n))

		hour_angle=table[2, 3:]
		declination=table[3:,2]

		for i in range(m):
			if i ==0:
				row_i=np.append(0, hour_angle)
			else:
				row_i=np.append(declination[i-1], table[2+i, 3:])

			f.write(" ".join(map(str, row_i)))
			f.write("\n")
		f.write("\n")
		
	f.close()


def append_oelts(table, identifier, motabfile, mac=None):
	"""Append oelt(s) to an existing .motab file, 
	   for example, the oelt(s) of windy conditions

	``Arguments``
		* table (numpy array): the additional oelt(s) that is to be appended  
		* identifier (str): the identified name of the new oelts 
		* motabfile (str): the directory of the existing.motab file
		* mac : MultiApertureConfiguration class instance, if it is a multi-aperture configuration

	``Return``
		write the new table(s) to the .motab file
	"""
	f=open(motabfile, 'a+')
	f.write("\n")
	if mac!=None:
		# multi-aperture configuration
		TABLE=table
		num_aperture=mac.n
		for i in range(num_aperture+1):
			table=TABLE[i]

			# size of the lookup table  
			m=np.shape(table)[0]-2
			n=np.shape(table)[1]-2

			if i==num_aperture:
				f.write('double optical_efficiency_total_%s(%s, %s)\n'%(identifier, m,n))
			else:
				lv=mac.get_lv_index(i)
				f.write('double optical_efficiency_level_%.0f_%s(%s, %s)\n'%(lv, identifier, m,n))

			hour_angle=table[2, 3:]
			declination=table[3:,2]

			for i in range(m):
				if i ==0:
					row_i=np.append(0, hour_angle)
				else:
					row_i=np.append(declination[i-1], table[2+i, 3:])

				f.write(" ".join(map(str, row_i)))
				f.write("\n")
			f.write("\n")
		

	else:
		# single aperture configuration
		m=np.shape(table)[0]-2
		n=np.shape(table)[1]-2
		f.write('double optics_%s(%s, %s)\n'%(identifier, m,n))

		hour_angle=table[2, 3:]
		declination=table[3:,2]

		for i in range(m):
			if i ==0:
				row_i=np.append(0, hour_angle)
			else:
				row_i=np.append(declination[i-1], table[2+i, 3:])

			f.write(" ".join(map(str, row_i)))
			f.write("\n")


	f.close()



def read_motab(filename, multi_aperture=False):

	with open(filename) as f:
		content=f.read().splitlines()
	f.close()
	res=content[4].split(',')
	nr=0
	if multi_aperture:

		eff_des=float(res[1])
		eff_annual=float(res[2])
		A_land=float(res[3])
		H_tower=float(res[3])
		A_helio=float(res[5])
		n_helios=float(res[6])
		Q_in_rcv=float(res[7])
		num_aperture=int(res[8])
		radial_distance=float(res[9])

		OELT={}
		Q_in_rcv_i=[]
		n_helios_i=[]
		H_rcv_i=[]
		W_rcv_i=[]
		Z_rcv_i=[]
		for i in range(num_aperture+1):
			if i<num_aperture:
				Q_in_rcv_i.append(float(res[10+i*5]))
				n_helios_i.append(float(res[11+i*5]))
				H_rcv_i.append(float(res[12+i*5]))
				W_rcv_i.append(float(res[13+i*5]))
				Z_rcv_i.append(float(res[14+i*5]))

			oelt=np.array([])
			solar_hour=np.array([])
			declination=np.array([])

			t=re.findall("[-+]?\d*\.\d+|\d+", content[6+i*(nr+3)])

			nr=int(t[-2])-1
			nc=int(t[-1])-1

			t=content[7+i*(nr+3)].split(' ')

			for v in t[1:]:
				solar_hour=np.append(solar_hour, float(v))

			for t in content[8+i*(nr+3): 8+nr+i*(nr+3)]:
				v=t.split(' ')
				declination=np.append(declination, float(v[0]))
				oelt=np.append(oelt, np.array(v[1:],dtype=float))

			oelt=oelt.reshape(len(declination), len(solar_hour))
			OELT[i]=oelt

		return eff_des, eff_annual, A_land, H_tower, A_helio, n_helios, Q_in_rcv, num_aperture, radial_distance, Q_in_rcv_i, n_helios_i, H_rcv_i, W_rcv_i, Z_rcv_i, solar_hour, declination, OELT

	else:
		n_helios=float(res[1])
		A_helio=float(res[2])
		eff_des=float(res[3])
		eff_annual=float(res[4])
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
		return n_helios, A_helio, eff_des, eff_annual, Q_in_rcv, A_land, solar_hour, declination, oelt



