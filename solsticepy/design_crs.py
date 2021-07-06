import os
import sys
import time
import numpy as np
import matplotlib.pyplot as plt
from uncertainties import ufloat
from scipy.interpolate import interp1d,interp2d
import matplotlib.cm as cm
from scipy.optimize import curve_fit

from .process_raw import *
from .cal_layout import radial_stagger
from .cal_field import *
from .cal_sun import *
from .gen_yaml import gen_yaml, Sun
from .gen_vtk import *
from .input import Parameters
from .output_motab import output_matadata_motab, output_motab
from .master import *


class CRS:
	'''
	The Central Receiver System (CRS) model includes three parts:
	the sun, the field and the receiver.
	'''

	def __init__(self, latitude, casedir, nproc=None, verbose=False):
		'''
		Arguements:
			casedir : str, the directory of the case
			nproc (int): number of processors, e.g. nproc=1 will run in serial mode,
                                                    nproc=4 will run with 4 processors in parallel
											        nproc=None will run with any number of processors that are available
			verbose : bool, write results to files or not
		'''
		self.casedir=casedir
		self.verb=verbose

		if not os.path.exists(casedir):
			os.makedirs(casedir)
		self.latitude=latitude
		self.sun=SunPosition()
		self.master=Master(casedir, nproc)

	def receiversystem(self, receiver, rec_w=0., rec_h=0., rec_x=0., rec_y=0., rec_z=100., rec_tilt=0., rec_grid_w=10, rec_grid_h=10, rec_abs=1., num_aperture=1, gamma=0.):

		'''
		Arguements:
		    (1) receiver  :   str, type of the receiver, i.e. 'flat', 'cylinder', 'multi_aperture' or directory of the 'stl',
		    (2) rec_w     : float, width of a flat receiver or radius of a cylindrical receiver (m)
		    (3) rec_h     : float, height of the receiver (m)
		    (4) rec_x     : float, x location of the receiver (m)
		    (5) rec_y     : float, y location of the receiver (m)
		    (6) rec_z     : float, z location of the receiver (m), or a list (for multi-aperture) of elevation heights of apertures
		    (7) rec_tilt  : float, tilt angle of the receiver (deg), 0 is where the receiver face to the horizontal
		    (8) rec_grid_w  :   int, number of elements in the horizontal(x)/circumferential direction
		    (9) rec_grid_h  :   int, number of elements in the vertical(z) direction
		    (10) rec_abs    : float, receiver surface absorptivity, e.g. 0.9
		    (11) num_aperture :   int, number of apertures if it is a multi-aperture receiver
		    (12) gamma   : float, the anangular range of the multi-aperture configration (deg)
		'''
		self.receiver=receiver
		self.rec_abs=rec_abs
		self.rec_w=rec_w
		self.rec_z=rec_z
		self.num_aperture=num_aperture
		self.gamma=gamma

		if receiver[-3:]=='stl':
			self.rec_param=[rec_w, rec_h, receiver, rec_x, rec_y, rec_z, rec_tilt]
		elif receiver=='flat':
			self.rec_param=[rec_w, rec_h, rec_grid_w, rec_grid_h, rec_x, rec_y, rec_z, rec_tilt]
		elif receiver=='cylinder':
			self.rec_param=[rec_w, rec_h, rec_grid_w, rec_grid_h, rec_x, rec_y, rec_z, rec_tilt]
		elif receiver=='multi-aperture':
			# rec_w and rec_h is the size of one aperture
			# rec_grid_w and rec_gird_h is the number of elements of one aperture
			# rec_x, rec_y, rec_z is the location of the central point of the multi-aperture receiver
			# rec_tilt is the tilt angle of each aperture
			# num_aperture is the number of apertures
			# gamma is the anangular range of the multi-aperture configration (deg)
			self.rec_param=[rec_w, rec_h, rec_grid_w, rec_grid_h, rec_z, rec_tilt, num_aperture, gamma]
			self.rec_w=rec_w



	def heliostatfield(self, field, hst_rho, slope, hst_w, hst_h, tower_h, tower_r=0.01, hst_z=0., num_hst=0., R1=0., fb=0., dsep=0.):

		'''
		Arguements:
		    (1) field     : str,
		        -- 'polar', 'polar-half', 'surround' or 'surround-half' or 'multi-aperture' for desiging a new field
		        -- or the directory of the layout file
		            the layout file is a 'csv' file, (n+2, 7)
		           - n is the total number of heliostats
		           - the 1st row is the number of each column
		           - the 2nd row is the unit
		           - the first three columns are X, Y, Z coordinates of each heliostat
		           - the fourth column is the focal length
		           - the last three columns are the aiming points
		    (2) hst_w     : float, width of a heliostat (m)
		    (3) hst_h     : float, height of a heliostat (m)
		    (4) hst_z     : float, the installation height of the heliostat
		    (5) hst_rho   : float, reflectivity of heliostat
		    (6) slope     : float, slope error(radians)
		    (7) R1        : float, layout parameter, the distance of the first row of heliostat
		    (8) dsep      : float, layout parameter, the separation distance of heliostats (m)
		    (9) tower_h    : float, tower height (m)
		    (10)tower_r   : float, radius of tower (m)
		    (11)num_hst  :   int, number of heliostats used in the field design

		 '''

		if field[-3:]=='csv':
			print('KNOWN FIELD')
			layout=np.loadtxt(field, delimiter=',', skiprows=2)

		else:
			# design a new field
			savefolder=self.casedir+'/des_point'
			if not os.path.exists(savefolder):
				os.makedirs(savefolder)

			pos_and_aiming, self.Nzones, self.Nrows=radial_stagger(latitude=self.latitude, num_hst=num_hst, width=hst_w, height=hst_h, hst_z=hst_z, towerheight=tower_h, R1=R1, fb=fb, dsep=0., field=field, num_aperture=self.num_aperture, gamma=self.gamma, rec_w=self.rec_w, rec_z=self.rec_z, savedir=savefolder, verbose=self.verb )

			layout=pos_and_aiming[2:, :]


		self.hst_w=hst_w
		self.hst_h=hst_h
		self.hst_rho=hst_rho
		self.slope=slope
		self.tower_h=tower_h
		self.tower_r=tower_r

		self.hst_pos=layout[:,:3].astype(float)
		self.hst_foc=layout[:,3].astype(float)
		self.hst_aims=layout[:,4:7].astype(float)

		self.hst_aim_idx=layout[:,7].astype(float)

		self.hst_zone=layout[:,9].astype(float)     # zone number
		self.hst_row=layout[:,10].astype(float)      # row index in the zone


	def yaml(self, dni=1000,sunshape=None,csr=0.01,half_angle_deg=0.2664,std_dev=0.2):
		'''
		Generate YAML files for the Solstice simulation
		'''
		outfile_yaml = self.master.in_case(self.casedir, 'input.yaml')
		outfile_recv = self.master.in_case(self.casedir, 'input-rcv.yaml')

		#att_factor=self.get_attenuation_factor()
		att_factor=1e-6
		print('attenuation', att_factor)
		sun = Sun(sunshape=sunshape, csr=csr, half_angle_deg=half_angle_deg, std_dev=std_dev)

		if self.latitude>0:
			hemisphere='North'
		else:
			hemisphere='South'
		gen_yaml(sun, self.hst_pos, self.hst_foc, self.hst_aims, self.hst_w
		, self.hst_h, self.hst_rho, self.slope, self.receiver, self.rec_param
		, self.rec_abs, outfile_yaml=outfile_yaml, outfile_recv=outfile_recv
		, hemisphere='North', tower_h=self.tower_h, tower_r=self.tower_r
		, spectral=False , medium=att_factor, one_heliostat=False)


	def field_design_annual(self,  dni_des, num_rays, nd, nh, weafile, method, Q_in_des=None, n_helios=None, zipfiles=False, gen_vtk=False, plot=False):
		'''
		Design a field according to the ranked annual performance of heliostats
		(DNI weighted)

		'''
		print('')
		print('Start field design')
		system=self.receiver

		AZI, ZENITH,table,case_list=self.sun.annual_angles(self.latitude, casefolder=self.casedir,nd=nd, nh=nh, verbose=self.verb)
		case_list=case_list[1:]
		SOLSTICE_AZI, SOLSTICE_ELE=self.sun.convert_convention('solstice', AZI, ZENITH)

		run=np.r_[0]
		nhst=len(self.hst_pos)

		#oelt=table
		ANNUAL=np.zeros(nhst)
		annual_solar=0.
		hst_annual={}

		for i in range(len(case_list)):
			c=int(case_list[i,0].astype(float))
			if c not in run:
				# the morning positions
				azimuth=SOLSTICE_AZI[c-1]
				elevation= SOLSTICE_ELE[c-1]

				if np.sin(elevation*np.pi/180.)>=1.e-5:
				    dni=1618.*np.exp(-0.606/(np.sin(elevation*np.pi/180.)**0.491))
				else:
				    dni=0.

				sys.stderr.write("\n"+green('Sun position: %s \n'%c))
				print('azimuth: %.2f'% azimuth, ', elevation: %.2f'%elevation)

				onesunfolder=os.path.join(self.casedir,'sunpos_%s'%(c))

				# run solstice
				if elevation<1.:
					efficiency_total=0
					performance_hst=np.zeros((nhst, 9))
					efficiency_hst=np.zeros(nhst)
				else:
					efficiency_total, performance_hst=self.master.run(azimuth, elevation, num_rays, self.hst_rho, dni, folder=onesunfolder, gen_vtk=gen_vtk, printresult=False, verbose=self.verb, system=system)

					#res=np.loadtxt(onesunfolder+'/result-formatted.csv', dtype=str, delimiter=',')
					#res_hst=np.loadtxt(onesunfolder+'/heliostats-raw.csv', dtype=str, delimiter=',')
					#efficiency_total=res[-2,1].astype(float)/res[1,1].astype(float)
					#performance_hst=res_hst[1:,-9:].astype(float)

					efficiency_hst=performance_hst[:,-1]/performance_hst[:,0]

				hst_annual[c]=performance_hst
				sys.stderr.write(yellow("Total efficiency: {:f}\n".format(efficiency_total)))
				run=np.append(run,c)

			cc=0
			for a in range(len(table[3:])):
				for b in range(len(table[0,3:])):
					val=re.findall(r'\d+', table[a+3,b+3])
					if str(c) in val:
						if cc==0:
							# i.e. morning positions
							ANNUAL+=dni*efficiency_hst
							annual_solar+=dni
							cc+=1
						else:
							# the symetrical points (i.e. afternoon)
							eff_symetrical=np.array([])
							for e in range(self.Nzones):
								idx_z=(self.hst_zone==e)
								eff_zone=efficiency_hst[idx_z]
								row_zone=self.hst_row[idx_z]

								nr=int(self.Nrows[e])
								for r in range(nr):
									idx_r=(row_zone==r)
									eff_row=eff_zone[idx_r]
									if r%2==0:
										eff_row=eff_row[::-1]
									else:
										eff_row[1:]=eff_row[1:][::-1]

									eff_symetrical=np.append(eff_symetrical, eff_row)

							#print(np.shape(eff_symetrical))
							#check=np.append(self.hst_zone, (self.hst_row, self.hst_num_idx, efficiency_hst, eff_symetrical))
							#print(np.shape(check))
							#check=check.reshape(5,int(len(check)/5))
							#np.savetxt('./check.csv', check.T, fmt='%.5f', delimiter=',')
							ANNUAL+=dni*eff_symetrical
							annual_solar+=dni

		ANNUAL/=annual_solar
		if self.verb:
			np.savetxt(self.casedir+'/annual_hst.csv',ANNUAL, fmt='%.2f', delimiter=',')

		designfolder=self.casedir+'/des_point'
		day=self.sun.days(21, 'Mar')
		dec=self.sun.declination(day)
		hra=0. # solar noon
		zen=self.sun.zenith(self.latitude, dec, hra)
		azi=self.sun.azimuth(self.latitude, zen, dec, hra)
		azi_des, ele_des=self.sun.convert_convention('solstice', azi, zen)

		sys.stderr.write("\n"+green('Design Point: \n'))
		efficiency_total, performance_hst_des=self.master.run(azi_des, ele_des, num_rays, self.hst_rho, dni_des, folder=designfolder, gen_vtk=gen_vtk, printresult=False, verbose=self.verb, system=system)

		#res=np.loadtxt(designfolder+'/result-formatted.csv', dtype=str, delimiter=',')
		#res_hst=np.loadtxt(designfolder+'/heliostats-raw.csv', dtype=str, delimiter=',')
		#efficiency_total=res[-2,1].astype(float)/res[1,1].astype(float)
		#performance_hst_des=res_hst[1:,-9:].astype(float)


		Qin=performance_hst_des[:,-1]
		Qsolar=performance_hst_des[0,0]

		#ID=ANNUAL.argsort()
		#ID=ID[::-1]
		ann_rank=ANNUAL/np.max(ANNUAL)
		ann_rank=np.around(ann_rank, decimals=1)
		#ID=ann_rank.argsort()
		#ID=ID[::-1]

		ID=np.lexsort((self.hst_foc,-ann_rank))
		#ID=np.lexsort((-ann_rank, self.hst_foc))

		if method==1:
			hst_aim_idx=self.hst_aim_idx[ID]
			print('')
			print('Method 1')
			self.Q_in_rcv=Q_in_des
			if self.receiver=='multi-aperture':
				self.Q_in_rcv_i=[] # the incident power on each aperture
				for ap in range(self.num_aperture):
					self.Q_in_rcv_i.append(0.)
			power=0.
			select_hst=np.array([])
			if self.receiver=='multi-aperture-individual':
				# selecting heliostats based on the required heat from individual receiver
				# initial selection
				assert isinstance(Q_in_des, list), "Q_in_des should be a list that specify the reuquired incident power to each aperture"

				for ap in range(self.num_aperture):
					power_i=0.
					idx_apt_i=(hst_aim_idx==ap)
					id_i=ID[idx_apt_i]

					for i in range(len(id_i)):
						if power_i<Q_in_des[ap]:
							idx=id_i[i]
							select_hst=np.append(select_hst, idx)
							power_i+=Qin[idx]
					power+=power_i
				self.Q_in_rcv_i=Q_in_des

			else:
				# initial selection
				# for single-aperture receiver
				# or multi-aperture receiver configuration that selects heliostats based on the total required heat
				assert isinstance(Q_in_des, float), "Q_in_des should be float, which is the total required incident power to the receiver"

				for i in range(len(ID)):
					if power<Q_in_des:
						idx=ID[i]
						select_hst=np.append(select_hst, idx)
						power+=Qin[idx]
						ap_idx=int(hst_aim_idx[i])
						self.Q_in_rcv_i[ap_idx]+=Qin[idx]


		else:
			select_hst=np.array([])
			print('')
			print('Method 2')
			#TODO the Method 2 does not include multi-aperture option
			num_hst=0
			power=0.
			for i in range(len(ID)):
			    if num_hst<n_helios:
			        idx=ID[i]

			        select_hst=np.append(select_hst, idx)
			        num_hst+=1
			        power+=Qin[idx]

			self.Q_in_rcv=power

		select_hst=select_hst.astype(int)

		self.hst_pos= self.hst_pos[select_hst,:]
		self.hst_foc=self.hst_foc[select_hst]
		self.hst_aims=self.hst_aims[select_hst,:]
		self.hst_aim_idx=self.hst_aim_idx[select_hst]

		self.n_helios=len(select_hst) # total number of heliostats
		self.eff_des=power/float(self.n_helios)/Qsolar

		print('num_hst', self.n_helios)
		print('power   @design', power)
		print('opt_eff @design', self.eff_des)
		Xmax=max(self.hst_pos[:,0])
		Xmin=min(self.hst_pos[:,0])
		Ymax=max(self.hst_pos[:,1])
		Ymin=min(self.hst_pos[:,1])
		A_land=(Xmax-Xmin)*(Ymax-Ymin)
		print('land area', A_land)

		if self.verb:
			title=np.array([['x', 'y', 'z', 'foc', 'aim x', 'aim y', 'aim z', 'aim_rec_index'], ['m', 'm', 'm', 'm', 'm', 'm', 'm', '-']])
			design_pos_and_aim=np.hstack((self.hst_pos, self.hst_foc.reshape(self.n_helios, 1)))
			design_pos_and_aim=np.hstack((design_pos_and_aim, self.hst_aims))
			design_pos_and_aim=np.hstack((design_pos_and_aim, self.hst_aim_idx.reshape(self.n_helios, 1)))
			#symmetric=design_pos_and_aim
			#symmetric[:, 0]=-symmetric[:, 0]
			#design_pos_and_aim=np.vstack((design_pos_and_aim, symmetric))
			#designed_field=design_pos_and_aim
			design_pos_and_aim=np.vstack((title, design_pos_and_aim))
			np.savetxt(self.casedir+'/pos_and_aiming.csv', design_pos_and_aim, fmt='%s', delimiter=',')
			np.savetxt(self.casedir+'/selected_hst.csv', select_hst, fmt='%.0f', delimiter=',')

		annual_solar=0.
		annual_field=0.

		oelt={}
		QTOT=np.zeros(np.shape(table))
		QIN=np.zeros(np.shape(table))
		self.n_helios_i=[]
		for ap in range(self.num_aperture):
			# lookup table
			print(ap)
			oelt[ap]=np.zeros(np.shape(table))

			idx_apt_i=(self.hst_aim_idx==ap)
			self.n_helios_i.append(np.sum(idx_apt_i))
			run=np.r_[0]

			print('')
			print('Aperture %s'%ap)
			print('num helios', np.sum(idx_apt_i))

			for i in range(len(case_list)):
				c=int(case_list[i,0].astype(float))
				if c not in run:
					#sundir=designfolder+'/sunpos_%s'%c
					res_hst=hst_annual[c]
					Qtot=res_hst[select_hst,0]
					Qin=res_hst[select_hst,-1]

					eff=np.sum(Qin[idx_apt_i])/np.sum(Qtot[idx_apt_i])

					print('sun position:', (c), 'eff', eff)

					azimuth=SOLSTICE_AZI[c-1]
					elevation= SOLSTICE_ELE[c-1]

					if np.sin(elevation*np.pi/180.)>=1.e-5:
						dni=1618.*np.exp(-0.606/(np.sin(elevation*np.pi/180.)**0.491))
					else:
						dni=0.

				for a in range(len(table[3:])):
					for b in range(len(table[0,3:])):

						val=re.findall(r'\d+',table[a+3,b+3])
						if val==[]:
							oelt[ap][a+3,b+3]=0
						else:
							if c==float(val[0]):
								oelt[ap][a+3,b+3]=eff
								QTOT[a+3,b+3]+=np.sum(Qtot[idx_apt_i])
								QIN[a+3,b+3]+=np.sum(Qin[idx_apt_i])
								annual_solar+=dni
								annual_field+=dni*eff


			oelt[ap][2, 3:]=table[2, 3:].astype(float)
			oelt[ap][3:,2]=table[3:,2].astype(float)


		self.eff_annual=annual_field/annual_solar

		if self.num_aperture==1:
			return oelt[0], A_land
		else:
			oelt[self.num_aperture]= np.divide(QIN, QTOT, out=np.zeros(QIN.shape, dtype=float), where=QTOT!=0)
			oelt[self.num_aperture][2, 3:]=table[2, 3:].astype(float)
			oelt[self.num_aperture][3:,2]=table[3:,2].astype(float)

			if self.num_aperture==3:
				oelt[0][3:,3+int(nh/2):]=oelt[2][3:, 3:3+int(nh/2)][:,::-1]
				oelt[2][3:,3+int(nh/2):]=oelt[0][3:, 3:3+int(nh/2)][:,::-1]

			if self.verb:
				for ap in range(self.num_aperture+1):
					if ap==self.num_aperture:
						np.savetxt(self.casedir+'/lookup_table_total.csv', oelt[ap], fmt='%s', delimiter=',')
					else:
						np.savetxt(self.casedir+'/lookup_table_%s.csv'%ap, oelt[ap], fmt='%s', delimiter=',')
			return oelt, A_land


	def annual_oelt(self, dni_des, num_rays, nd, nh, zipfiles=False, gen_vtk=False, plot=False):
		'''
		Annual performance of a known field
		'''
		self.n_helios=len(self.hst_pos)
		oelt, ANNUAL=self.master.run_annual(nd=nd, nh=nh, latitude=self.latitude, num_rays=num_rays, num_hst=self.n_helios,rho_mirror=self.hst_rho, dni=dni_des, verbose=self.verb)

		Xmax=max(self.hst_pos[:,0])
		Xmin=min(self.hst_pos[:,0])
		Ymax=max(self.hst_pos[:,1])
		Ymin=min(self.hst_pos[:,1])
		A_land=(Xmax-Xmin)*(Ymax-Ymin)
		print('land area', A_land)

		if self.verb:
			np.savetxt(self.casedir+'/lookup_table.csv', oelt, fmt='%s', delimiter=',')


		designfolder=self.casedir+'/des_point'
		day=self.sun.days(21, 'Mar')
		dec=self.sun.declination(day)
		hra=0. # solar noon
		zen=self.sun.zenith(self.latitude, dec, hra)
		azi=self.sun.azimuth(self.latitude, zen, dec, hra)
		azi_des, ele_des=self.sun.convert_convention('solstice', azi, zen)

		sys.stderr.write("\n"+green('Design Point: \n'))
		efficiency_total, performance_hst_des=self.master.run(azi_des, ele_des, num_rays, self.hst_rho, dni_des, folder=designfolder, gen_vtk=False, printresult=False, verbose=self.verb)
		self.eff_des=efficiency_total.n

		return oelt, A_land


	def dni_TMY(self, weafile, nd, nh, plot=False):
		'''
		Argument:

		weafile : str, directory of the weather file .motab
		# col0 - time (s)
		# col1 - GHI
		# col2 -DNI (W/m2)
		# col3 -DHI
			...
		'''
		with open(weafile) as f:
			content=f.read().splitlines()
		f.close()

		lines=len(content)
		seconds=np.array([])
		dni=np.array([])
		for i in range(2, lines-1):
			l=content[i].split(",")
			seconds=np.append(seconds, l[0])
			dni=np.append(dni, l[2])
		seconds=seconds.astype(float)+1800.
		days=seconds/3600/24

		wea_dec=np.array([])
		for d in days:
			d=int(d)+1
			wea_dec=np.append(wea_dec, self.sun.declination(d)) #deg

		wea_hra=((seconds/3600.)%24-12.)*15. #deg
		wea_dni=dni.astype(float)

		dh=360./float(nh)
		dd=23.45*2./float(nd)

		hra_lim=180.+dh/2.
		dec_lim=23.45+dd/2.

		hra_bin=np.linspace(-hra_lim, hra_lim, nh+1)
		dec_bin=np.linspace(-dec_lim, dec_lim, nd+1)
		bins=np.array([hra_bin, dec_bin])

		dni_weight, xbins, ybins=np.histogram2d(wea_hra, wea_dec, bins, weights=wea_dni)
		dni_n, xbins, ybins=np.histogram2d(wea_hra, wea_dec, bins)

		dni_avg=np.divide(dni_weight, dni_n, out=np.zeros_like(dni_weight), where=dni_n!=0)

		if plot:
			X=np.linspace(-180.,  180. , nh)
			Y=np.linspace(-23.45, 23.45, nd)
			plt.pcolormesh(hra_bin, dec_bin, dni_avg.T)
			cb=plt.colorbar()
			plt.xlabel('Solar hour angle')
			plt.ylabel('Delination angle')
			plt.savefig(self.casedir+'/DNI_weather.png', bbox_inches='tight')
			plt.close()

			# check the symmetricity of the weather dni data
			data=dni_avg.T
			plt.pcolormesh(hra_bin[:int(nh/2)], dec_bin, data[:,:int(nh/2), ]-np.fliplr(data[:,-int(nh/2):]))
			cb=plt.colorbar()
			plt.xlabel('Solar hour angle')
			plt.ylabel('Delination angle')
			plt.savefig(self.casedir+'/DNI_weather_diff.png', bbox_inches='tight')
			plt.close()

			np.savetxt(self.casedir+'/DNI_weather.csv', data, fmt='%.4f', delimiter=',')

		return dni_weight.T, dni_avg.T

	def get_attenuation_factor(self):

		foc=self.hst_foc.astype(float)

		# to get the attenuation factor
		def func(x, b):
		    return np.exp(-b * x)
		def fun_two(x):
		    return 0.99321-0.0001176*x+1.97e-8*x**2
		xdata = np.linspace(0, np.max(foc), np.max(foc)*100)
		y = fun_two(xdata)
		ydata = y
		popt, pcov = curve_fit(func, xdata, ydata)
		y2 = [func(i, popt[0]) for i in xdata]
		att_factor =popt[0]
		return att_factor

	def generateVTK(self,eta_hst, savevtk):

		FPF=FieldPF(0., 0., np.r_[0., 1., 0.])
		norms=np.zeros(np.shape(self.hst_pos))
		norms[:,-1]=1.
		COORD, TRI, element, nc=FPF.view_heliostats(width=self.hst_w, height=self.hst_h, normals=norms, hstpos=self.hst_pos)
		NORMS=np.repeat(norms, element, axis=0)

		#field performance
		hst_tot  =eta_hst[:,0]
		rec_abs=eta_hst[:,1]

		hst_eff=rec_abs/hst_tot

		savevtk=savevtk+'/results-field.vtk'
		TOT=np.repeat(hst_tot, element)
		ABS=np.repeat(rec_abs, element)
		EFF=np.repeat(hst_eff, element)

		DATA={'tot': TOT, 'rec_abs': ABS, 'efficiency':EFF}
		gen_vtk(savedir=savevtk, points=COORD.T, indices=TRI, norms=NORMS, colormap=True, DATA=DATA)

if __name__=='__main__':
	start=time.time()
	casedir='./test-crs-design'
	tablefile=casedir+'/OELT_Solstice.motab'
	if os.path.exists(tablefile):
		print('')
		print('Load exsiting OELT')

	else:

		pm=Parameters()
		pm.Q_in_rcv=565e6
		pm.nd=5
		pm.nh=5
		pm.H_tower=250.
		pm.H_rcv=24.
		pm.W_rcv=24.
		pm.dependent_par()
		pm.saveparam(casedir)
		print(pm.fb)
		print(pm.H_tower)
		crs=CRS(latitude=pm.lat, casedir=casedir)
		weafile='/home/yewang/solartherm-master/SolarTherm/Data/Weather/gen3p3_Daggett_TMY3.motab'
		crs.heliostatfield(field=pm.field_type, hst_rho=pm.rho_helio, slope=pm.slope_error, hst_w=pm.W_helio, hst_h=pm.H_helio, tower_h=pm.H_tower, tower_r=pm.R_tower, hst_z=pm.Z_helio, num_hst=pm.n_helios, R1=pm.R1, fb=pm.fb, dsep=pm.dsep)

		crs.receiversystem(receiver=pm.rcv_type, rec_w=float(pm.W_rcv), rec_h=float(pm.H_rcv), rec_x=float(pm.X_rcv), rec_y=float(pm.Y_rcv), rec_z=float(pm.Z_rcv), rec_tilt=float(pm.tilt_rcv), rec_grid_w=int(pm.n_W_rcv), rec_grid_h=int(pm.n_H_rcv),rec_abs=float(pm.alpha_rcv))

		crs.yaml(sunshape=pm.sunshape,csr=pm.crs,half_angle_deg=pm.half_angle_deg,std_dev=pm.std_dev)

		oelt, A_land=crs.field_design_annual(dni_des=900., num_rays=int(1e6), nd=pm.nd, nh=pm.nh, weafile=weafile, method=1, Q_in_des=pm.Q_in_rcv, n_helios=None, zipfiles=False, gen_vtk=False, plot=False)


		if (A_land==0):
		    tablefile=None
		else:
		    A_helio=pm.H_helio*pm.W_helio
		    output_matadata_motab(table=oelt, field_type=pm.field_type, aiming='single', n_helios=self.n_helios, A_helio=A_helio, eff_design=self.eff_des, H_rcv=pm.H_rcv, W_rcv=pm.W_rcv, H_tower=pm.H_tower, Q_in_rcv=pm.Q_in_rcv, A_land=A_land, savedir=tablefile)


	end=time.time()
	print('total time %.2f'%((end-start)/60.), 'min')
