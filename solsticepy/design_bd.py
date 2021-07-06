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


class BD:
	'''
        The Beam-Down (BD) model includes five parts:
        the sun, the field, the secondary reflector, the CPC and the receiver.
	'''

	def __init__(self, latitude, casedir):
		'''
		Arguements:
			casedir : str, the directory of the case
		'''
		self.casedir=casedir

		if not os.path.exists(casedir):
			os.makedirs(casedir)
		self.latitude=latitude
		self.sun=SunPosition()
		self.master=Master(casedir)

	def cpcradius(self, rec_w=2., rec_l=2., cpc_nfaces=4):
		'''
		Ouput
		rec_radius_ref:	receiver radius in the direction of theta (CPC half acceptance angle), it is used for the design criteria
		rec_radius:			receiver radius in the opposite direction (equal to rec_radius_theta if nber of CPC faces is different to 4, because symmetry applies)
		'''
		if cpc_nfaces is 4:
			rec_radius_ref = np.minimum(rec_w, rec_l)/2.
			rec_radius = np.maximum(rec_w, rec_l)/2.
		else:
			rec_radius = np.sqrt(rec_w*rec_w+rec_l*rec_l)/2.
			rec_radius_ref = rec_radius

		return rec_radius_ref, rec_radius

	def hyperboloidparameters(self, rec_z, cpc_theta_deg, cpc_h, field_rim_angle, aim_z):
		'''
		'''
		#intersection point of line passing by aiming point and furthest heliostat & line passing by secondary real foci (aperture of cpc) forming an angle theta with the vertical (cpc acceptance angle)
		theta = cpc_theta_deg*np.pi/180.
		phy = field_rim_angle*np.pi/180.
		x_inter = ((rec_z+cpc_h) - aim_z) / (1/np.tan(phy) - 1/np.tan(theta))
		z_inter = x_inter/np.tan(phy)+aim_z
		# translate the intersection point to place the mid distance of teh 2 focci at the origin
		foci_dist = (aim_z - (rec_z+cpc_h))/2.
		z_inter -= foci_dist

		assert z_inter>0., 'there is no corresponding hyperbola to fit the design criteria'

		foci_dist_sqrt = foci_dist**2
		z_inter_sqrt = z_inter**2
		x_inter_sqrt = x_inter**2
		# c0, c1, c2 are the hyperbola equation parameters with the distance to the vertex as the unknown parameter
		c2 = 1.
		c1 = -(z_inter_sqrt + x_inter_sqrt + foci_dist_sqrt)
		c0 = z_inter_sqrt*foci_dist_sqrt
		delta = c1**2 - 4*c0*c2
		vertex_dist_sqrt = (-c1-np.sqrt(delta))/2.
		check = z_inter_sqrt/vertex_dist_sqrt - x_inter_sqrt/(foci_dist_sqrt-vertex_dist_sqrt)
		print('equation resolution should be 1, and it is: ', check)
		vertex_dist = np.sqrt(vertex_dist_sqrt)

		return vertex_dist, x_inter, foci_dist

	def intersectionlinehyperbol(self, a_hyper, b_hyper, m_line, z0_line):
		'''
		intersection point between hyperbola and line:
		z^2/a^2 - y^2/b^2 = 1
		z = m*y + z0
		'''
		b_hyper_sqrt = b_hyper**2
		a_hyper_sqrt = a_hyper**2
		m_line_sqrt = m_line**2
		aCoef_poly = b_hyper_sqrt - a_hyper_sqrt/m_line_sqrt
		bCoef_poly = 2*z0_line*a_hyper_sqrt/m_line_sqrt
		cCoef_poly = -a_hyper_sqrt*((z0_line**2)/m_line_sqrt+b_hyper_sqrt)
		delta = bCoef_poly**2 - 4*aCoef_poly*cCoef_poly
		root1 = (-bCoef_poly+np.sqrt(delta))/(2*aCoef_poly)
		root2 = (-bCoef_poly-np.sqrt(delta))/(2*aCoef_poly)
		if root1>0. and root2>0.:	# the line intersect the upper sheet of the hyperbola twice
			z_inter = min(root1,root2)
		else:	# the line intersect the upper and lower sheet of the hyperbola
			z_inter = max(root1,root2)
		y_inter = (z_inter-z0_line)/m_line
		check = (z_inter**2)/a_hyper_sqrt - (y_inter**2)/b_hyper_sqrt
		print('equation resolution should be 1, and it is: ', check)

		return y_inter

	def zhyperboloid(self, x_hyper, y_hyper, vertex_dist, foci_dist):
		'''
		calculate the z coordinate of the hyperboloid:
		(x^2+y^2)/a^2 - (z+z0-g/2)^2/b^2 + 1 = 0
		with the vetex of the hyperboloid located in the xOy plan.
		'''
		g_para = 2*foci_dist
		f_para = (vertex_dist+foci_dist)/g_para
		a_para_sqrt = (g_para**2)*(f_para-f_para**2)
		b_para = g_para*(f_para-1/2.)
		z0_hyper = abs(b_para) + g_para/2.

		z_hyper = b_para*np.sqrt((x_hyper**2+y_hyper**2)/a_para_sqrt+1) - z0_hyper+g_para/2.
		return z_hyper

	def receiversystem(self, receiver, rec_abs=1., rec_w=1.2, rec_l=10., rec_z=0., rec_grid=200, cpc_nfaces=4, cpc_theta_deg=20., cpc_h_ratio=1.,
	cpc_nZ=20., field_rim_angle=30., aim_z=62., secref_inv_eccen=None, secref_vert=np.array([[-15,25],[-15,-25],[15,-25],[15,25]]), refl_sec=0.95, slope_error=0.0):

		'''
		Variables:
		- cpc_theta_deg
		- field_rim_angle
		- aim_z (tower height)
		- rec_z (optional)
		Arguments:
		    (1) receiver  :   str, type of the receiver; 'beam-down',
		    # Arguments for flat receiver
		    (2) rec_abs    : float, receiver surface absorptivity, e.g. 0.9
		    (3) rec_w     : float,  width of the receiver (m)
		    (4) rec_l     : float, length of the receiver (m)
		    (5) rec_z     : float, z (vertical) location of the receiver (m)
		    (6) rec_grid  : number of slices for solstice flux map
		    # Arguments for Compound Parabolic Concentrator (CPC)
		    (7) cpc_nfaces : int, number of faces of the CPC
		    (8) cpc_theta_deg : float, acceptance angle of CPC (deg)
		    (9) cpc_h_ratio     : float, ratio of critical CPC height calculated with cpc_theta_deg
		    in local coordinate system of the parabola in the xOy plan
		    WARNING: cpc_h is different from the total height of the CPC
		    (10) cpc_nZ     : int, number of number of incrementation for the clipping polygon for the construction of each CPC face
		    # Arguments for Secondary Reflector
		    (11) field_rim_angle : float, rim angle of the field formed wih vertical axis in the zOx plan (deg)
		    (12) aim_z   : float, z (vertical) location of the heliostats' aiming point (m)
		    (13) secref_inv_eccen    : flaot, hyperboloid inverse eccentricity: ratio of the apex distance over the foci distance to the origin, must be between 0 and 1
		    (14) secref_vert	: array, clipping polygon of the secondary reflector
		    (15) r_f	: float, secondary mirror and CPC reflectivity, e.g. 0.9
		    (16) slope_error	: float, slope error of secondary mirror and CPC refelctivity (?)
		'''
		self.receiver=receiver
		self.rec_abs=rec_abs

		assert cpc_nfaces > 2, 'The number of faces for the CPC should be minimum 3, and it is {cpc_nfaces=}'
		if secref_inv_eccen != None:
			assert 0<=secref_inv_eccen<=1, 'The inverse eccentricity of the hyperbole must be between 0 and 1, and it is {secref_inv_eccen=}'
		rec_rad_ref, rec_rad = self.cpcradius(rec_w, rec_l, cpc_nfaces)

		# Calculate the maximum height of the CPC if not defined
		cpc_theta = cpc_theta_deg*np.pi/180.
		cpc_h = rec_rad_ref*(1+1/np.sin(cpc_theta))/np.tan(cpc_theta)*cpc_h_ratio

		assert aim_z > (cpc_h+rec_z), 'The imaginary foci of the hyperbol is lower than its real foci'

		# Calculate the characteristics of the secondary reflector
		if secref_inv_eccen is None:
			vertex_dist, x_inter, foci_dist = self.hyperboloidparameters(rec_z, cpc_theta_deg, cpc_h, field_rim_angle, aim_z)
		else:
			foci_dist = (aim_z-(cpc_h+rec_z))/2.
			vertex_dist = secref_inv_eccen*foci_dist
			x_inter = None
		secref_z = rec_z + cpc_h + vertex_dist + foci_dist
		print('hyperboloid distances ratio: ', vertex_dist/foci_dist)
		print('hyperboloid eccentricity: ', foci_dist/vertex_dist)

		# Calculate the clipping polygon of the secondary reflector
		if secref_vert is None:
			b_hyper = np.sqrt(foci_dist**2-vertex_dist**2)
			if  x_inter is None:
				m_line = -1/np.tan(field_rim_angle*np.pi/180.)
				x_inter = self.intersectionlinehyperbol(vertex_dist, b_hyper, m_line, foci_dist)

			if rec_rad_ref != rec_rad:
				m_line = cpc_h/rec_rad
				# if the coefficient of the line is outside the tangent of the hyperbola, there will be no intersection
				hyper_tangent = vertex_dist/b_hyper
				if abs(m_line) <= hyper_tangent:
					y_inter = rec_rad + (x_inter-rec_rad_ref)
					print('tangent')
				else:
					y_inter = self.intersectionlinehyperbol(vertex_dist, b_hyper, m_line, -foci_dist)
			else:
				y_inter = x_inter
			secref_vert = np.array([[-x_inter,y_inter],[-x_inter,-y_inter],[x_inter,-y_inter],[x_inter,y_inter]])

		# calculate the field maxium radius along x and y axis
		self.x_max = aim_z*np.tan(field_rim_angle*np.pi/180.)
		y_inter = max(secref_vert[:,1])
		z_hyper = secref_z + self.zhyperboloid( 0.0, y_inter, vertex_dist, foci_dist)
		if z_hyper<aim_z:
			phy = np.arctan(y_inter/(aim_z-z_hyper))
		else:
			phy = 80.*np.pi/180.
		self.y_max = aim_z*np.tan(phy)

		self.rec_param=np.array([rec_w, rec_l, rec_z, rec_grid, cpc_nfaces, cpc_theta_deg, cpc_h, cpc_nZ, aim_z, secref_z, refl_sec, slope_error, secref_vert])

	def heliostatfield(self, field, hst_rho, slope, hst_w, hst_h, tower_h, tower_r=0.01, hst_z=0., num_hst=0., R1=0., fb=0., dsep=0., x_max=-1.0, y_max=-1.0):

		'''
                Arguments:
		    (1) field     : str,
		        -- 'polar', 'polar-half', 'surround' or 'surround-half' for desiging a new field
		        -- the directory of the layout file
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
		    (12)Q_in_rcv :   int, required heat of the receiver in the field design
		 '''

		if field[-3:]=='csv':
		    print('KNOWN FIELD')
		    layout=np.loadtxt(field, delimiter=',', skiprows=2)

		else:
		    # design a new field
		    savefolder=self.casedir+'/des_point'
		    if not os.path.exists(savefolder):
		        os.makedirs(savefolder)

		    if (x_max>0) and (y_max>0):
		        x_max = self.x_max
		        y_max = self.y_max
		        r_max = np.sqrt(x_max**2 + y_max**2)*1.5
		        pos_and_aiming=radial_stagger(latitude=self.latitude, num_hst=num_hst, width=hst_w, height=hst_h, hst_z=hst_z, towerheight=tower_h, R1=R1, fb=fb, dsep=0., field=field, savedir=savefolder, plot=False, r_field_max= r_max)
		    else:
		        pos_and_aiming=radial_stagger(latitude=self.latitude, num_hst=num_hst, width=hst_w, height=hst_h, hst_z=hst_z, towerheight=tower_h, R1=R1, fb=fb, dsep=0., field=field, savedir=savefolder, plot=False)

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
		self.hst_row=layout[:,-1].astype(float)

		Xmax = max(self.hst_pos[:,0])
		if (x_max>R1) and (Xmax>x_max*1.5):
			select_hst=np.array([])
			for i in range(len(self.hst_pos)):
				Xhelio = self.hst_pos[i,0]
				if np.abs(Xhelio)<x_max*1.5:
					select_hst=np.append(select_hst, i)

			select_hst=select_hst.astype(int)

			self.hst_pos= self.hst_pos[select_hst,:]
			self.hst_foc=self.hst_foc[select_hst]
			self.hst_aims=self.hst_aims[select_hst,:]
			self.hst_row=self.hst_row[select_hst]

		Ymax = max(self.hst_pos[:,1])
		if (y_max>R1) and (Ymax>y_max*1.5):
			select_hst=np.array([])
			for i in range(len(self.hst_pos)):
				Yhelio = self.hst_pos[i,1]
				if np.abs(Yhelio)<y_max*1.5:
					select_hst=np.append(select_hst, i)

			select_hst=select_hst.astype(int)

			self.hst_pos= self.hst_pos[select_hst,:]
			self.hst_foc=self.hst_foc[select_hst]
			self.hst_aims=self.hst_aims[select_hst,:]
			self.hst_row=self.hst_row[select_hst]

		np.savetxt(self.casedir+'/trimmed_field.csv', self.hst_pos, fmt='%s', delimiter=',')


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
		dni_weight=self.dni_TMY(weafile, nd, nh)

		AZI, ZENITH,table,case_list=self.sun.annual_angles(self.latitude, casefolder=self.casedir, nd=nd, nh=nh)
		case_list=case_list[1:]
		SOLSTICE_AZI, SOLSTICE_ELE=self.sun.convert_convention('solstice', AZI, ZENITH)

		oelt=table
		run=np.r_[0]
		nhst=len(self.hst_pos)
		ANNUAL=np.zeros(nhst)
		hst_annual={}

		for i in range(len(case_list)):
			c=int(case_list[i,0].astype(float))
			if c not in run:
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
					efficiency_total, performance_hst=self.master.run(azimuth, elevation, num_rays, self.hst_rho, dni, folder=onesunfolder, gen_vtk=gen_vtk, printresult=False, system='beamdown')

					efficiency_hst=performance_hst[:,-1]/performance_hst[:,0]


				hst_annual[c]=performance_hst
				sys.stderr.write(yellow("Total efficiency: {:f}\n".format(efficiency_total)))

			cc=0
			for a in range(len(table[3:])):
				for b in range(len(table[0,3:])):
					val=re.findall(r'\d+', table[a+3,b+3])
					if val==[]:
						table[a+3,b+3]=0
					else:
						if c==float(val[0]):
							#table[a+3,b+3]=efficiency_total # this line cause wired problem
							if cc==0: # avoid adding the symetrical point
								ANNUAL+=dni_weight[a,b]*efficiency_hst
								cc+=1
			run=np.append(run,c)
		np.savetxt(self.casedir+'/annual_hst.csv',ANNUAL, fmt='%.2f', delimiter=',')

		designfolder=self.casedir+'/des_point'
		day=self.sun.days(21, 'Mar')
		dec=self.sun.declination(day)
		hra=0. # solar noon
		zen=self.sun.zenith(self.latitude, dec, hra)
		azi=self.sun.azimuth(self.latitude, zen, dec, hra)
		azi_des, ele_des=self.sun.convert_convention('solstice', azi, zen)

		sys.stderr.write("\n"+green('Design Point: \n'))
		efficiency_total, performance_hst_des=self.master.run(azi_des, ele_des, num_rays, self.hst_rho, dni_des, folder=designfolder, gen_vtk=gen_vtk, printresult=False, system='beamdown')

		Qin=performance_hst_des[:,-1]
		Qsolar=performance_hst_des[0,0]

		np.savetxt(self.casedir+'/annual_hst_design.csv',Qin, fmt='%.2f', delimiter=',') # New, to check the NRJ at design point

		ID=ANNUAL.argsort()
		ID=ID[::-1]

		select_hst=np.array([])
		if method==1:
			print('')
			print('Method 1')
			power=0.
			self.Q_in_rcv=Q_in_des
			for i in range(len(ID)):
				if power<Q_in_des:
					idx=ID[i]
					if Qin[idx]>0.:
						select_hst=np.append(select_hst, idx)
						power+=Qin[idx]

		else:
			print('')
			print('Method 2')
			num_hst=0
			power=0.
			for i in range(len(ID)):
			    if num_hst<n_helios:
			        idx=ID[i]

			        select_hst=np.append(select_hst, idx)
			        num_hst+=1
			        power+=Qin[idx]

		self.Q_in_rcv=power
		print('')
		print(power)
		print('')
		select_hst=select_hst.astype(int)

		self.hst_pos=self.hst_pos[select_hst,:]
		self.hst_foc=self.hst_foc[select_hst]
		self.hst_aims=self.hst_aims[select_hst,:]

		self.n_helios=len(select_hst)
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

		num_hst=len(self.hst_pos)
		title=np.array([['x', 'y', 'z', 'foc', 'aim x', 'aim y', 'aim z'], ['m', 'm', 'm', 'm', 'm', 'm', 'm']])
		design_pos_and_aim=np.hstack((self.hst_pos, self.hst_foc.reshape(num_hst, 1)))
		design_pos_and_aim=np.hstack((design_pos_and_aim, self.hst_aims))
		#symmetric=design_pos_and_aim
		#symmetric[:, 0]=-symmetric[:, 0]
		#design_pos_and_aim=np.vstack((design_pos_and_aim, symmetric))
		#designed_field=design_pos_and_aim
		design_pos_and_aim=np.vstack((title, design_pos_and_aim))
		np.savetxt(self.casedir+'/pos_and_aiming.csv', design_pos_and_aim, fmt='%s', delimiter=',')
		np.savetxt(self.casedir+'/selected_hst.csv', select_hst, fmt='%.0f', delimiter=',')


		# lookup table
		run=np.r_[0]

		for i in range(len(case_list)):
			c=int(case_list[i,0].astype(float))
			if c not in run:
			    #sundir=designfolder+'/sunpos_%s'%c
			    res_hst=hst_annual[c]
			    Qtot=res_hst[:,0]
			    Qin=res_hst[:,-1]
			    Qtot_sum=np.sum(Qtot[select_hst])
			    if Qtot_sum>0.:
			        eff=np.sum(Qin[select_hst])/Qtot_sum
			    else:
			        eff=0.
			    print('')
			    print('sun position:', (c), 'eff', eff)

			for a in range(len(oelt[3:])):
			    for b in range(len(oelt[0,3:])):
			        val=re.findall(r'\d+',oelt[a+3,b+3])
			        if val==[]:
			            oelt[a+3,b+3]=0
			        else:
			            if c==float(val[0]):
			                oelt[a+3,b+3]=eff

		np.savetxt(self.casedir+'/lookup_table.csv', oelt, fmt='%s', delimiter=',')

		return oelt, A_land


	def annual_oelt(self, dni_des, num_rays, nd, nh, zipfiles=False, gen_vtk=False, plot=False):
		'''
		Annual performance of a known field
		'''
		self.n_helios=len(self.hst_pos)
		oelt, ANNUAL=self.master.run_annual(nd=nd, nh=nh, latitude=self.latitude, num_rays=num_rays, num_hst=self.n_helios, rho_mirror=self.hst_rho, dni=dni_des, system='beamdown')

		Xmax=max(self.hst_pos[:,0])
		Xmin=min(self.hst_pos[:,0])
		Ymax=max(self.hst_pos[:,1])
		Ymin=min(self.hst_pos[:,1])
		A_land=(Xmax-Xmin)*(Ymax-Ymin)
		print('land area', A_land)

		np.savetxt(self.casedir+'/lookup_table.csv', oelt, fmt='%s', delimiter=',')

		designfolder=self.casedir+'/des_point'
		day=self.sun.days(21, 'Mar')
		dec=self.sun.declination(day)
		hra=0. # solar noon
		zen=self.sun.zenith(self.latitude, dec, hra)
		azi=self.sun.azimuth(self.latitude, zen, dec, hra)
		azi_des, ele_des=self.sun.convert_convention('solstice', azi, zen)

		sys.stderr.write("\n"+green('Design Point: \n'))
		efficiency_total, performance_hst_des=self.master.run(azi_des, ele_des, num_rays, self.hst_rho, dni_des, folder=designfolder, gen_vtk=True, printresult=False, system='beamdown')
		self.eff_des=efficiency_total.n

		return oelt, A_land


	def dni_TMY(self, weafile, nd, nh):
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
		seconds=seconds.astype(float)
		days=seconds/3600/24
		wea_dec=np.array([])
		for d in days:
			wea_dec=np.append(wea_dec, self.sun.declination(d)) #deg

		wea_hra=((seconds/3600.)%24-12.)*15. #deg
		wea_dni=dni.astype(float)

		hra_lim=180.*(float(nh)/float(nh-1))
		dec_lim=23.45*(float(nd)/float(nd-1))
		hra_bin=np.linspace(-hra_lim, hra_lim, nh+1)
		dec_bin=np.linspace(-dec_lim, dec_lim, nd+1)
		bins=np.array([hra_bin, dec_bin])


		dni_weight, xbins, ybins=np.histogram2d(wea_hra, wea_dec, bins, weights=wea_dni)

		return dni_weight.T

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
		hst_tot=eta_hst[:,0]
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
		bd=BD(latitude=pm.lat, casedir=casedir)
		weafile='./demo_TMY3_weather.motab'

		# Polygon receiver
		rec_w=10.
		rec_l=10.
		rec_z=-60.
		rec_grid=200
		# 2D-crossed cpc with n faces
		n_CPC_faces=4
		theta_deg=20.
		n_Z=30
		# Secondary refector 'hyperboloid'
		rim_angle = 45.
		refl_sec = 0.95
		slope_error = 0.0


		bd.receiversystem(receiver='beam_down', rec_abs=float(pm.alpha_rcv), rec_w=float(rec_w), rec_l=float(rec_l), rec_z=float(rec_z), rec_grid=int(rec_grid), cpc_nfaces=int(n_CPC_faces), cpc_theta_deg=float(theta_deg), cpc_h=None, cpc_nZ=float(n_Z), field_rim_angle=float(rim_angle), aim_z=float(pm.H_tower), secref_inv_eccen=None, secref_vert=None, refl_sec=float(refl_sec), slope_error=float(slope_error))

		bd.heliostatfield(field='surround', hst_rho=pm.rho_helio, slope=pm.slope_error, hst_w=pm.W_helio, hst_h=pm.H_helio, tower_h=pm.H_tower, tower_r=pm.R_tower, hst_z=pm.Z_helio, num_hst=num_hst, R1=pm.R1, fb=pm.fb, dsep=pm.dsep, x_max=150., y_max=150.)

		bd.yaml(sunshape=pm.sunshape,csr=pm.crs,half_angle_deg=pm.half_angle_deg,std_dev=pm.std_dev)

		oelt, A_land=bd.field_design_annual(dni_des=900., num_rays=int(1e6), nd=pm.nd, nh=pm.nh, weafile=weafile, method=1, Q_in_des=pm.Q_in_rcv, n_helios=None, zipfiles=False, gen_vtk=False, plot=False)


		if (A_land==0):
		    tablefile=None
		else:
		    A_helio=pm.H_helio*pm.W_helio
		    output_matadata_motab(table=oelt, field_type='surround', aiming='single', n_helios=self.n_helios, A_helio=A_helio, eff_design=self.eff_des, H_rcv=pm.H_rcv, W_rcv=pm.W_rcv, H_tower=pm.H_tower, Q_in_rcv=pm.Q_in_rcv, A_land=A_land, savedir=tablefile)


	end=time.time()
	print('total time %.2f'%((end-start)/60.), 'min')
