import numpy as np

class SunPosition:

	def __init__(self):
		""" Calculate sun position according to a certain location, date and time. Reference: John A. Duffie and William A. Beckman. Solar Engineering of Thermal Processes, 4th edition.
		
		``Example``

		  * Location: latitude = 37.44
		  * Date    : 21 Jun
		  * Time    : solar noon or two hours after sunrise
									
			>>> from solsticepy.cal_sun import *
			>>> sun=SunPosition()
			>>> day=sun.days(21, 'Jun')
			>>> latitude=37.44
			>>> delta=sun.declination(day)

			Solar noon

			>>> omega=0.
			>>> theta=sun.zenith(latitude, delta, omega)
			>>> phi=sun.azimuth(latitude, theta, delta, omega)
			>>> print(theta)
			13.987953925483858
			>>> print(phi)
			0.0

			Two hours after sunrise

			>>> daytime, sunrise=sun.solarhour(delta, latitude)
			>>> omega=sunrise+15.*2. # solar noon
			>>> theta=sun.zenith(latitude, delta, omega)
			>>> phi=sun.azimuth(latitude, theta, delta, omega)
			>>> print(theta)
			67.91774797592434
			>>> print(phi)
			-103.31434583806747 

		"""		
		pass


	def days(self, dd, mm):
		"""Calculate the day-of-the-year for specified date and month, ref. J Duffie page 14, Table 1.6.1

		``Arguments``

			* dd (int): the day in the month
			* mm (str): the month, e.g. 'Jan', 'Feb', Mar' etc

		``Return`` 

			* days (int): the-day-of-the year for the specified dd-mm, e.g. 1 Jan is day 1

		"""

		if mm=='Jan':
		    days=dd

		elif mm=='Feb':
		    days=31+dd

		elif mm=='Mar':
		    days=59+dd

		if mm=='Apr':
		    days=90+dd

		elif mm=='May':
		    days=120+dd

		elif mm=='Jun':
		    days=151+dd

		if mm=='Jul':
		    days=181+dd

		elif mm=='Aug':
		    days=212+dd

		elif mm=='Sep':
		    days=243+dd

		if mm=='Oct':
		    days=273+dd

		elif mm=='Nov':
		    days=304+dd

		elif mm=='Dec':
		    days=334+dd

		return days


	def declination(self, days, form='detail'):
		"""Calculate the solar declination angle for a specified day-of-the-year, ref. J Duffie page 13, declination angle: delta=23.45*sin(360*(284+day)/365)

		``Arguments``

		  * day (int): day of the year, i.e from 1 to 365
		  * form (str): 'detail' or simple' model

		``Return``

		  * delta (float): declination angle (deg)
		"""

		if form=='detail':
		    B=float(days-1)*360./365.*np.pi/180.

		    delta=(180./np.pi)*(0.006918 - 0.399912*np.cos(B) +0.070257*np.sin(B)- 0.006758*np.cos(2.*B) + 0.000907*np.sin(2.*B)- 0.002697*np.cos(3.*B) + 0.00148*np.sin(3.*B))

		else:
		    delta=23.45*np.sin(360.*float(284+days)/365.*np.pi/180.) # deg

		return delta


	def solarhour(self, delta, latitude):
		"""Calculate day length and sunrise hour-angle, ref. J Duffie page 17

		``Arguments``

		  * delta (float): declination angle (deg)
		  * latitude (float): latitude angle (deg)

		``Returns``

		  * hour (float): length of the daylight hour (h)
		  * sunrise (float): the solar hour angle at sunrise (deg)

		"""

		sunset=np.arccos(-np.tan(latitude*np.pi/180.)*np.tan(delta*np.pi/180.))*180./np.pi # deg

		sunrise=-sunset

		hour=(sunset-sunrise)/15.

		return hour, sunrise


	def zenith(self, latitude, delta, omega):
		"""Calculate the zenith angle, ref. J Duffie eq.1.6.5

		``Arguments``

		  * latitude (float): latitude angle at the location (deg)
		  * delta (float):  declination angle (deg)
		  * omega (float): solar hour angle (deg)

		``Return``

		  * theta (float): the zenith angle (deg)

		"""            
		latitude*=np.pi/180.
		delta*=np.pi/180.
		omega*=np.pi/180.

		theta=np.arccos(np.cos(latitude)*np.cos(delta)*np.cos(omega)+np.sin(latitude)*np.sin(delta))*180./np.pi

		return theta
		
	def azimuth(self, latitude, theta, delta, omega):
		"""Calculate azimuth angle at specified location and sun position, ref. J Duffie eq. 1.6.6

		``Arguments``

		  * latitude (float): latitude angle (deg)
		  * delta (float): declination angle (deg)
		  * theta (float): zenith angle (deg)
		  * omega (float): solar hour angle (deg)

		``Return``
 
		  * phi (float): azimuth angle (deg), counted from South towards to West 
		"""
		latitude*=np.pi/180.
		delta*=np.pi/180.
		theta*=np.pi/180.

		a1=np.cos(theta)*np.sin(latitude)-np.sin(delta)
		a2=np.sin(theta)*np.cos(latitude)
		b=a1/a2

		if abs(b+1.)<1e-10:
		    phi=np.pi

		elif abs(b-1.)<1e-10:
		    phi=0.
		else:
		    phi=abs(np.arccos((np.cos(theta)*np.sin(latitude)-np.sin(delta))/(np.sin(theta)*np.cos(latitude)))) # unit radian

		if omega<0:
		    phi=-phi

		phi*=180./np.pi

		return phi

	def convert_AZEL_to_declination_hour(self, theta, phi, latitude):
		""" Convert azimuth-elevation angle to declination-hour angle

		``Arguments``

		  * theta (float): zenith angle (deg)
		  * phi (float): azimuth angle (deg), counted from South towards to West
		  * latitude (float): latitude latitude (deg)

		``Returns``

		  * delta: declination angle (deg)
		  * omega: solar hour angle (deg)
		"""
     
		phi*=np.pi/180.
		theta*=np.pi/180.
		latitude*=np.pi/180.

		delta=np.arcsin(np.cos(theta)*np.sin(latitude)-np.cos(abs(phi))*np.sin(theta)*np.cos(latitude))

		omega=np.arccos((np.cos(theta)-np.sin(latitude)*np.sin(delta))/(np.cos(latitude)*np.cos(delta)))
		if phi<0:
		    omega=-omega

		delta*=180./np.pi
		omega*=180./np.pi

		return delta, omega

	def convert_convention(self, tool, azimuth, zenith):
		"""Return azimuth-elevation angle using the angle convention of Solstice or SolarTherm

		``Arguments``

		  * tool (str): 'solstice' or 'solartherm'
		  * azimuth (float): azimuth angle (deg), counted from South towards to West
		  * zenith (float): zenith angle (deg)


		``Returns``

		  * sol_azi (float): azimuth angle 
		  * sol_ele (float): elevation angle
		"""
		if tool=='solstice':
		    # azimuth: from East to North 
		    sol_azi=-(90.+azimuth)
		    sol_ele=90.-zenith 

		elif tool=='solartherm':
		    # azimuth: from North to Ease (clockwise)
		    # ref Ali email (20/3/18)
		    sol_azi=180.+azimuth
		    sol_ele=90.-zenith  

		if isinstance(sol_azi, np.ndarray):
		    for i in range(len(sol_azi)):
		        azi=sol_azi[i]
		        ele=sol_ele[i]
		        if (azi>=360. or azi<0.):
		            sol_azi[i]=(azi+360.)%360.
		        if ele<=1e-20:
		            sol_ele[i]=0.  
		else:
		    if (sol_azi>=360. or sol_azi<0.):
		        sol_azi=(sol_azi+360.)%360.
		    if sol_ele<=1e-20:
		        sol_ele=0.  
		
		return sol_azi, sol_ele        
		                   

	def annual_angles(self, latitude, casefolder='NOTSAVE', nd=5, nh=5):
		"""Generate a range of sun positions (azimuth-zenith angles and declination-solarhour angles) for annual optical lookup table simulation. Automatically detect the time when the sun is below the horizon (elevation<0), where a ray-tracing simulation is not required.
		
		``Arguments``

		  * latitude (float): latitude of the location (deg)
		  * casefolder (str): directory to save the table and case_list in .csv files, or 'NOTSAVE' (by default) to not write the output to files
		  * nd (int): number of rows of the lookup table (points in the declination movement, suggest nd>=5)
		  * nh (int): number of columns of the lookup table (hours in a day, i.e. 24h)

		``Returns``

		  * AZI (1 D numpy array): a list of azimuth angles (deg), counted from South towards to West
		  * ZENITH (1D numpy array): a list of zenith angles (deg)
		  * table (numpy array): declination (row) - solarhour (column) lookup table to be simulated 
		  * case_list (numpy array): a flatten list of cases to be simulated, with the correspondence between declination-solarhour angles and azimuth-zenith angles for each case


		``Example``

			>>> from solsticepy.cal_sun import *
			>>> sun=SunPosition()
			>>> latitude=37.44
			>>> casefolder='.'
			>>> nd=5
			>>> nh=9
			>>> sun.annual_angles(latitude, casefolder, nd, nh)

			A 5x9 lookup table will be generated to the current directory:

			+--------------+----+-------+-------+--------+--------+--------+--------+--------+--------+--------+--------+
			| Lookup table |    |       |                      Solar hour angles (deg)                                  | 
			+==============+====+=======+=======+========+========+========+========+========+========+========+========+
			|              |    |       |   1   |    2   |   3    |    4   |    5   |    6   |   7    |    8   |    9   |
			+              +----+-------+-------+--------+--------+--------+--------+--------+--------+--------+--------+
			|              |    |       | -180  |  -135  |  -90   |  -45   |    0   |   45   |   90   |   135  |  180   |
			+              +----+-------+-------+--------+--------+--------+--------+--------+--------+--------+--------+
			|              | 1  |-23.45 |  `-`  |   `-`  |   `-`  | case 1 | case 2 | ***1   |   `-`  |   `-`  |   `-`  |
			+ Declination  +----+-------+-------+--------+--------+--------+--------+--------+--------+--------+--------+
			|              | 2  |-11.73 |  `-`  |   `-`  |   `-`  | case 3 | case 4 | ***3   |   `-`  |   `-`  |   `-`  |
			+ angles       +----+-------+-------+--------+--------+--------+--------+--------+--------+--------+--------+
			|              | 3  |   0   |  `-`  |   `-`  | case 5 | case 6 | case 7 | ***6   | ***5   |   `-`  |   `-`  |
			+ (deg)        +----+-------+-------+--------+--------+--------+--------+--------+--------+--------+--------+
			|              | 4  | 11.73 |  `-`  |   `-`  | case 8 | case 9 | case 10| ***8   | ***9   |   `-`  |   `-`  |
			+              +----+-------+-------+--------+--------+--------+--------+--------+--------+--------+--------+
			|              | 5  | 23.45 |  `-`  |   `-`  | case 11| case 12| case 13| ***11  | ***12  |   `-`  |   `-`  |
			+--------------+----+-------+-------+--------+--------+--------+--------+--------+--------+--------+--------+

			  * case `n`: this sun position is the `n`-th case to be simulated
			  * ***`n`  : it is the symmetric case with the case `n`, no re-simulation is required
			  * `-`     : the sun is below the horizon, no simulation is required
		"""

		# declination angle (deg)  
		# -23.45 ~ 23.45
		DELTA=np.linspace(-23.45, 23.45, nd)
 
		# solar time
		solartime=np.linspace(-180., 180., nh) # deg

		table=np.zeros(((nh+3)*(nd+3)))
		table=table.astype(str)
		for i in range(len(table)):
			table[i]=''

		table=table.reshape(nd+3, nh+3)
		table[0,0]='Lookup table'
		table[3,0]='Declination (deg)'
		table[0,3]='Hour angle (deg)'

		table[3:,1 ]=np.arange(1,nd+1)
		table[1 ,3:]=np.arange(1,nh+1)

		table[3:,2 ]=DELTA   
		table[2 ,3:]=solartime

		c=1
		AZI=np.array([])
		ZENITH=np.array([])

		case_list=np.array(['Case','declination (deg)','solar hour angle (deg)', 'azimuth (deg) S-to-W ', 'zenith (deg)'])

		for i in range(nd):
			delta=DELTA[i]
			hour, sunrise=self.solarhour(delta, latitude)
			sunset=-sunrise
			for j in range(nh):
				omega=solartime[j]
				if (omega>sunset or omega<sunrise):
					table[3+i,3+j]='-' 

				else:
					if omega<0:

						table[3+i, 3+j]=' case %s'%(c)
						table[3+i, -(1+j)]='***%s'%(c)

						#zenith angle
						theta=self.zenith(latitude, delta, omega)
						# azimuth        
						phi=self.azimuth(latitude, theta, delta, omega)
						AZI=np.append(AZI, phi)
						ZENITH=np.append(ZENITH, theta)

						case_list=np.append(case_list, (c, delta, omega, phi, theta)) 

						#zenith angle (afternoon)
						theta=self.zenith(latitude, delta, -omega)
						# azimuth        
						phi=self.azimuth(latitude, theta, delta, -omega)    
						case_list=np.append(case_list, (c, delta, -omega, phi, theta))                 

						c+=1

					elif omega==0:
						table[3+i, 3+j]=' case %s'%(c)

						#zenith angle
						theta=self.zenith(latitude, delta, omega)
						# azimuth        
						phi=self.azimuth(latitude, theta, delta, omega)
						AZI=np.append(AZI, phi)
						ZENITH=np.append(ZENITH, theta)

						case_list=np.append(case_list, (c, delta, omega, phi, theta)) 
						c+=1

		              	                                     
		case_list=case_list.reshape(int(len(case_list)/5),5)
		#azimuth=case_list[1:,-2].astype(float)
		#zenith=case_list[1:,-1].astype(float)

		if casefolder!='NOTSAVE':    
		    np.savetxt(casefolder+'/table_view.csv', table, fmt='%s', delimiter=',')  
		    np.savetxt(casefolder+'/annual_simulation_list.csv', case_list, fmt='%s', delimiter=',')          

		return AZI, ZENITH, table,case_list



if __name__=='__main__':
    # example: PS10, summer solstice, solar noon
    latitude=37.44

    sun=SunPosition()
    dd=sun.days(21, 'Jun')
    delta=sun.declination(dd)
    print('Declination angle', delta)


    daytime,sunrise=sun.solarhour(delta, latitude)
    print('Day timeS', daytime)
    print('sun rise', sunrise)
    print('')

    omega=0. # solar noon
    theta=sun.zenith(latitude, delta, omega)
    phi=sun.azimuth(latitude, theta, delta, omega)
    print('elevation', 90.-theta)
    print('azimuth', phi)

    azi, zen, table,caselist=sun.annual_angles(latitude, hemisphere='North', casefolder='.',nd=5, nh=7)


    
    



        
        
        
