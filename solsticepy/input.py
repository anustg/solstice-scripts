import numpy as np
import os

class Parameters:

	def __init__(self):
		'''
		default parameters

		e.g. a 50 MWe power plant (Yogi Goswami, Principles of Solar Engineering, Third Edition, page 480, section 8.8)
			 design point DNI 950 W/m2, SM=1.8
			 the required field rating is 305 MW
			 the total reflector area is approximated as 483,000 m2
			 tower height 137 m for a surround field, 183 m for a polar field

			for PS10, the annual capacity of the field is around 120GWh
		'''

		self.simulation()
		self.Sun()
		self.Heliostat()
		self.Receiver()
		self.dependent_par()

	def Sun(self):
		'''
			(1) lat       : float, latitude of the plant lcation
			(2) dni_des       : float, dni at design point (W/m2)   
			(3) sunshape  :   str, 'pillbox' or 'Buie' 
			(4) sunsize   : float, 
				-- for 'pillbox': half angle of the pillbox sunshape (degree)
				-- for 'Buie': circumsolar ratio (CSR)  

			(5) extinction: float, around 1e-6, the extinction coefficient if the atmosphere is surrounded by participant medium
		'''
		self.lat=37.44 # latitude, default: location of PS10
		self.dni_des=1000.
		self.sunshape='pillbox' # pillbox or buie
		# csr for buie, half_angle for pillbox
		self.sunshape_param=4.65*1.e-3*180./np.pi # convert rad to degree --> solstice convention
		self.extinction=1e-6
		self.wea_file=None


	def Heliostat(self):
		'''
			(1) field_type : str,
				-- 'polar', 'surround', 'multi-aperture', 'polar-half' or 'surround-half for desiging a new field 
				    the 'half' option is for simulation of a symmetric field
				-- the directory of the layout file
				    the layout file is a 'csv' file, (n+2, 7)
				   - n is the total number of heliostats 
				   - the 1st row is the number of each column
				   - the 2nd row is the unit 
				   - the first three columns are X, Y, Z coordinates of each heliostat
				   - the fourth column is the focal length
				   - the last three columns are the aiming points
			(2) Q_in_rcv      : float, required heat of the receiver (W)
			(3) W_helio       : float, width of a heliostat (m) 
			(4) H_helio       : float, height of a heliostat (m)
			(5) slope_error   : float, slope error(radians)
			(6) helio_rho     : float, reflectivity of heliostat 
			(7) helio_soil    : float, percentage of the heliostat surface that is not soiled
			(8) helio_sf_ratio: float, percentage of avaiable heliostat reflective surface area 
			(9) H_tower       : float, tower height (m)
			(10) R_tower      : float, radius of tower (m)
			(11) concret_tower:  bool, True-a solid concrete tower, or False-a truss tower
			(12)single_field  :  bool, True-one tower one field, or False-multi-tower
			(13)R1            : float, layout parameter, the distance of the first row of heliostat 
			(14)dsep          : float, layout parameter, the separation distance of heliostats (m)
			(15)fb            : float, in (0-1), a factor to expand the field to reduce block 
			---(*) n_helios   :   int, number of heliostats for the designed field 
			---(*) Z_helio    : float, the installation height of the heliostat
		'''

		self.field_type='polar'
		if self.method==1:
			self.Q_in_rcv=10e6 # required heat of the receiver  
		else:
			self.n_helios=1000
		self.W_helio=10.
		self.H_helio=10.
		self.slope_error=1.53e-3 # radian
		self.slope_error_windy=2.e-3 # radian, a largher optical error in windy conditions
		self.windy_optics=0 # simulate the windy oelt or not? 1 is yes, 0 is no
		#self.helio_rho=0.95 # heliostat reflectivity
		#self.helio_soil=0.95 # soiling factor of heliostats
		#self.helio_sf_ratio=0.97 # heliostat reflective surface availability
		#self.helio_refl=self.helio_rho*self.helio_soil*self.helio_sf_ratio
		self.helio_refl=0.9
		self.H_tower=100.
		self.R_tower=0.001  # shading effect of tower is neglected at the moment
		self.concret_tower=False
		self.single_field=True 
		self.R1=90.
		self.dsep=0.
		self.fb=0.7
		self.aimingstrategy=0 # use some sophisticated aiming strategy, e.g. MDBA that used for sodium receivers, 1 is yes, 0 is no
		self.aim_pm1=0. # parameter 1 in aiming strategy
		self.aim_pm2=0. # parameter 2 in aiming strategy
		self.f_oversize=1. # oversising factor
		self.delta_r2=0. # field expanding for zone2
		self.delta_r3=0. # field expanding for zone3
		self.SM=0. # solar multiple
		

	def Receiver(self):
		'''
			(1) rcv_type  :   str, 'flat', 'cylinder', 'particle', 'multi-aperture' or 'stl', type of the receiver
			(2) num_aperture: int, number of apertures if it is a multi-aperture receiver
			(3) alpha  : float, the angular space between two adjacent apertures (except the aperture faces to the South) (deg)	 
			(4) H_rcv     : float, height of the receiver
			(5) W_rcv     : float, width of a flat receiver or diameter of a cylindrical receiver (m)
			(6) tilt_rcv  : float, tilt angle of the receiver (deg), 0 is where the receiver face to the horizontal
			(7) alpha_rcv : float, receiver surface absorptivity (0-1), set as 1 if rcv_type='particle'
			(8) n_H_rcv   :   int, number of discretisation of the receiver in the height direction
			(9) n_W_rcv   :   int, number of discretisation of the receiver in the width direction (for a flat receiver) 
				                    or in the circular direction (for a cylindrical receiver) 
			(10) X_rcv    : float, x location of the receiver (m)
			(11) Y_rcv    : float, y location of the receiver (m)
			(12) Z_rcv    : float, z location of the receiver (m)
			(13) num_aperture: int, number of apertures for a multi-aperture configuration
			(14) gamma    : float, the angular range of the multi-aperture configuration (deg)

		'''
		self.rcv_type='flat'
		self.H_rcv=10. #  height of the receiver (m)
		self.W_rcv=10. #  width of a flat receiver or diameter of a cylindrical receiver
		self.tilt_rcv=0. # receiver tilt angle
		self.alpha_rcv=1. 
		self.n_H_rcv=10
		self.n_W_rcv=10
		self.X_rcv=0. # receiver location
		self.Y_rcv=0.
		#self.Z_rcv=None
		self.num_aperture=1
		self.gamma=0.
		self.therm=0 # integration with receiver thermal performance, 1 is yes, 0 is no
		self.Nb=0. # number of banks
		self.Nfp=0. # number of flow paths
		self.Do=0. # tube outer diameter
		self.fluxlimitpath='' # directory of the files for flux limits
		self.T_in=290+273.15 # K, receiver inlet temperature
		self.T_out=565+273.15 # K, receiver outlete temperature
		self.HTF='salt' # 'salt' or 'sodium'
		self.rcv_material='Incoloy800H' # 'Haynes230' or 'Incoloy800H' or 'Inconel740H'

	def simulation(self):
		'''
		n_row_oelt: int, number of rows of the lookup table (i.e. simulated days per year)
		n_col_oelt: int, number of columns of the lookup table (i.e. simulated number of hours per day)
		n_rays    : int, number of rays for the simulation
		n_procs   : int, number of processors for the mcrt simulation 
		casedir   : str, the directory for saving the result files
		'''  
		self.n_row_oelt=5
		self.n_col_oelt=5
		self.n_rays=int(5e6)
		self.n_procs=0
		self.casedir='.'
		self.method=1 # 1 - design the field based on the Q_in_rcv
				      # 2 - design the field based on the n_helios
		self.verbose=0 # save all the simulation details or not? 1 is yes, 0 is no
		self.gen_vtk=0 # visualise the simulation scene or not? 1 is yes, 0 is no	
    

	def dependent_par(self):
		'''
		'''
		self.Z_helio=self.H_helio*0.7       
		if self.lat>=0:
			self.hemisphere='North'
		elif self.lat<0:
			self.hemisphere='South'

		if self.num_aperture==1:
			self.Z_rcv=self.H_tower


		# estimate a rough number of large field
		eta_field=0.4 # assumed field effieicy at design point
		if self.method==1:
			if isinstance(self.Q_in_rcv, float):
				self.n_helios=self.Q_in_rcv/self.W_helio/self.H_helio/self.dni_des/eta_field
			else:
				self.n_helios=sum(self.Q_in_rcv)/self.W_helio/self.H_helio/self.dni_des/eta_field
			  
			if self.field_type!='surround':
				self.n_helios*=1.5  


	def saveparam(self, savedir):
		if not os.path.exists(savedir):
			os.makedirs(savedir)
    
		param=np.array([
				['method', self.method, '-'],    
				['windy optics', bool(self.windy_optics), '-'],   
				['','',''], 
				['field', self.field_type, '-'],  
				['Q_in_rcv', self.Q_in_rcv, 'W'],   
				['SM', self.SM, 'W'],    
				['n_helios(pre_des if method ==1)', self.n_helios, '-'],    
				['W_helio', self.W_helio, 'm'],    
				['H_helio', self.H_helio, 'm'],
				['helio effective reflectivity', self.helio_refl, '-'],     
				['slope_error', self.slope_error, 'rad'],
				['slope_error_windy', self.slope_error_windy, 'rad'],
				['H_tower', self.H_tower, 'm'],    
				['R_tower', self.R_tower, 'm'], 
				['concret_tower', self.concret_tower, '-'],    
				['single_field', self.single_field, '-'],  
				['fb factor', self.fb, '-'],     
				['R1', self.R1, 'm'],    
				['dsep', self.dsep, 'm'],  
				['delta_r2', self.delta_r2, 'm'],  
				['delta_r3', self.delta_r3, 'm'],  
				['aiming strategy',self.aimingstrategy,''],
				['self.aim_pm1', self.aim_pm1, ''],
				['self.aim_pm2', self.aim_pm2, ''],
				['f_oversize', self.f_oversize, ''],
				['','',''],   
				['rcv_type', self.rcv_type, '-'],  
				['num_aperture', self.num_aperture, '-'],
				['aperture angular range',self.gamma , 'deg'],
				['H_rcv', self.H_rcv, 'm'],
				['W_rcv', self.W_rcv, 'm'],   
				['tilt_rcv', self.tilt_rcv, 'deg'],  
				['alpha_rcv', self.alpha_rcv, '-'],  
				['n_H_rcv', self.n_H_rcv, '-'],   
				['n_W_rcv', self.n_W_rcv, '-'],    
				['X_rcv', self.X_rcv, 'm'],   
				['Y_rcv', self.Y_rcv, 'm'],   
				['Z_rcv', self.Z_rcv, 'm'],   
				['','',''], 
				['receiver thermal pm', self.therm, ''],
				['Nb',self.Nb, ''],
				['Nfp', self.Nfp, ''],
				['Do', self.Do, ''],
				['fluxlimitpath', self.fluxlimitpath, ''],
				['receiver material', self.rcv_material, ''],
				['heat transfer fluid', self.HTF, ''],
				['rcv inlet temperature', self.T_in, 'K'],
				['rcv outlet temperature', self.T_out, 'K'],
				['','',''], 			
				['n_row_oelt', self.n_row_oelt, '-'] , 
				['n_col_oelt', self.n_col_oelt, '-'] ,
				['n_rays', self.n_rays, '-'] ,
				['n_procs', self.n_procs, '-'] 
				])
		np.savetxt(savedir+'/simulated_parameters.csv', param, delimiter=',', fmt='%s')
    

