from srcPy.gen_YAML import gen_YAML
from srcPy.cal_sun import *
from srcPy.get_raw import proces_raw_results
import numpy as N
import os
import matplotlib.pyplot as plt
from uncertainties import ufloat
import re

# please SET:
# set the directory of the Solstice software in your local system
solstice_dir='/home/yewang/Solstice-0.8.1-GNU-Linux64' 
# set the folder for saving the current ray-tracing case
casefolder='./2-PS10/annual'
#
#
# the sun
# =========
# S1. DNI
DNI=1000 # W/m2
# S2. sunshape
sunshape='pillbox' # or 'buie'
sunsize=0.2664 # the half angle of the pillbox sunshape distribution, in degree 
               # or CSR value of Buie sunshape

# S3. sun position - annual performance 

# set the mesh of the sky
nd=7 # number of days in a year (minimum is 5)
nh=7 # number of time in a day

latitude=37.44 # latitude at the location of the plant
#

# S4. number of rays for the ray-tracing simulation
num_rays=1000
#
# the field
# ==========     
# F1.Layout
layout=N.loadtxt('./PS10_layout.csv', delimiter=',', skiprows=2)
hst_pos=layout[:,:3]
hst_foc=layout[:,3] # F2.5
hst_aims=layout[:,4:] # F4.
# F2. Heliostat
hst_w=12.84 # m
hst_h=9.45 # m
rho_refl=0.88 # mirror reflectivity
slope_error=2.9e-3 # radians
#
# the receiver
# ============
# R1. shape
receiver='flat'
# R2. Size
rec_w=16. # width, m
rec_h=15. # height, m
# R3. tilt angle
tilt=0. # deg
# R4. position
loc_x=0. # m
loc_y=0. # m
loc_z=107.7# m
# R5. Abosrptivity
rec_abs=0.9
# receiver mesh, for binning the flux distribution
rec_mesh=100


#
# ===============================================================
# the lines below will automatically create the ray-tracing scene
# based on the settings above

if not os.path.exists(casefolder):
    os.makedirs(casefolder) 
# generate the YAML file (just once)
rec_param=N.r_[rec_w, rec_h, rec_mesh, loc_x, loc_y, loc_z, tilt]
gen_YAML(DNI, sunshape, sunsize, hst_pos, hst_foc, hst_aims,hst_w, hst_h, rho_refl, slope_error, receiver, rec_param, rec_abs, casefolder, spectral=False, medium=0., OneHeliostat=False )

# annual sun position calculation
#===========================================
sun=SunPosition()
azimuth, zenith,table,caselist=sun.annual_angles(latitude, hemisphere='North', casefolder=casefolder ,nd=nd, nh=nh)
# convert to solstice convention
sol_azi=270.+azimuth
sol_ele=90.-zenith                             
for i in xrange(len(sol_azi)):
    azi=sol_azi[i]
    ele=sol_ele[i]
    if (azi>=360. or azi<0.):
        sol_azi[i]=(azi+360.)%360.
    if ele<=1e-20:
        sol_ele[i]=0.
#==========================================

ANNUALRES=N.array([ 'Solstice_Azi (deg)','Solstice_Elevation (deg)','Optical efficiecy +/- errors'])
EFF=N.array([])

cases=caselist[1:]
run=N.r_[0]
for i in xrange(len(cases)):    
    c=int(cases[i,0].astype(float))

    if c not in run:
        azi=sol_azi[c-1]
        ele=sol_ele[c-1]
        print ''
        print ''
        print 'sun position:', (c)
        print 'azimuth', azi
        print 'elevation', ele
        onesunfolder=casefolder+'/sunpos_%s'%(c)
        if not os.path.exists(onesunfolder):
            os.makedirs(onesunfolder) 
        # run solstice
        os.system('solstice -D%s,%s -v -n %s -R %s/input-rcv.yaml -fo %s/simul %s/input.yaml'%(azi, ele, num_rays,casefolder, onesunfolder, casefolder))
        rawfile=onesunfolder+'/simul'
        savedir=onesunfolder
        rho_mirror=rho_refl
        eff=proces_raw_results(rawfile, savedir,rho_mirror)
        ANNUALRES=N.append(ANNUALRES, (azi, ele, eff.nominal_value))

    else:
        ANNUALRES=N.append(ANNUALRES, ('-', '-', eff.nominal_value))

    for a in xrange(len(table[3:])):
        for b in xrange(len(table[0,4:])):
            val=re.findall(r'\d+',    table[a+3,b+4])
            if len(val)!=0:
                if c==float(val[0]):
                    table[a+3,b+4]=eff.nominal_value
            else:
                    table[a+3,b+4]=0

    run=N.append(run,c)
ANNUALRES=ANNUALRES.reshape(len(ANNUALRES)/3, 3)
results=N.hstack((caselist, ANNUALRES))

N.savetxt(casefolder+'/annual_simulation_results.csv', results, fmt='%s', delimiter=',')
N.savetxt(casefolder+'/table_view_results.csv', table, fmt='%s', delimiter=',')

# plot

eta=table[3:,4:].astype(float)
declination=table[2,4:].astype(float)
solarhour=table[3:,3].astype(float)

dd=declination[1]-declination[0]
ds=solarhour[1]-solarhour[0]
dec=N.linspace(declination[0]-dd/2., declination[-1]+dd/2., nd+1)
sh=N.linspace(solarhour[0]-ds/2., solarhour[-1]+ds/2.,  2*nh)

plt.figure(1)
ax =plt.subplot()
im=ax.pcolormesh(dec, sh, eta)
plt.colorbar(im)
plt.tight_layout()
plt.xlabel('Declination angle (deg)')
plt.ylabel('Solar hour angle (deg)')
ax.set_title("Annual optical efficiency \n Lookup table")
plt.savefig(casefolder+"/eta_annual.png", bbox_inches='tight')
plt.close()    

