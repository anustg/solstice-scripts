#
# This is an example script to simulate a central receiver system (CRS)
# using Solstice via solsticepy
#
import solsticepy
from solsticepy.design_bd import *
from solsticepy.gen_vtk import read_vtk_annual
import numpy as np
import os
from glob import glob

#==================================================
# INPUT PARAMETERS

# define a unique case folder for the user
# =========

snum = 0
suffix = ""
while 1:
    import datetime,time
    dt = datetime.datetime.now()
    ds = dt.strftime("%a-%H-%M")
    casefolder = os.path.join(os.getcwd(),'case-%s-%s%s'%(os.path.basename(__file__),ds,suffix))
    if os.path.exists(casefolder):
        snum+=1
        suffix = "-%d"%(snum,)
        if snum > 200:
            raise RuntimeError("Some problem with creating casefolder")
    else:
        # good, we have a new case dir
        break

pm=Parameters()

# Variables
pm.H_tower=150.453433 # 80. # tower height or vertical distance to aiming point (located at the center of xOy plan)
pm.cpc_theta_deg=10.516461   # acceptance half angle of the CPC in degree
pm.cpc_h_ratio=1.0 # must be inferior or equal to 1
pm.rim_angle_x = 70.038185 # rim angle of the heliostat field in the xOz plan in degree
pm.rim_angle_y = 67.097593 # rim angle of the heliostat field in the xOz plan in degree
pm.secref_inv_eccen = 0.622431 # hyperboloid inverse eccentricity, [0,1]
pm.Z_rcv=0.
pm.fb=0.16153
pm.tilt_secref = 0. # angle of the tilted axis of the hyperboloid, from the vertical to the North (+) or South (-), [-180,180]


# fixed parameters
# =========
pm.sunshape='buie'
pm.crs=0.02
pm.Q_in_rcv=40e6
pm.n_row_oelt=5
pm.n_col_oelt=5
pm.H_rcv=10.
pm.W_rcv=1.2
pm.W_helio=6.1 #1.84 # ASTRI helio size
pm.H_helio=6.1 # 2.44
pm.field_type='surround'
pm.Z_helio=0.
pm.R1=15.
pm.lat=-27.85 #degree
print(pm.fb)
print(pm.dsep)
print(pm.R1)

# enter the parameters for the beam-down components
pm.rcv_type='beam_down'
# 2D-crossed cpc with n faces
pm.cpc_nfaces=4
pm.n_H_rcv=20
pm.cpc_nZ=30
# Properties for Mirrors surfaces
pm.rho_helio=0.9
pm.rho_secref=0.95
pm.rho_cpc=0.95
pm.slope_error=1.e-3 #heliostats # radian
pm.slope_error_bd=1.e-3 #beam-down components # radian

pm.saveparam(casefolder)
# create the environment and scene
# =========
tablefile=casefolder+'/OELT_Solstice.motab'
#weafile='./demo_TMY3_weather.motab'
weafile='./AUS_WA_Leinster_Airport_954480_TMY.motab'
bd=BD(latitude=pm.lat, casedir=casefolder)

# Design the field: Create an oversized field and trim the oversized field to the requested power at design point
bd.receiversystem(receiver=pm.rcv_type, rec_abs=float(pm.alpha_rcv), rec_w=float(pm.W_rcv), rec_l=float(pm.H_rcv), rec_z=float(pm.Z_rcv), rec_grid=int(pm.n_H_rcv),
cpc_nfaces=int(pm.cpc_nfaces), cpc_theta_deg=float(pm.cpc_theta_deg), cpc_h_ratio=float(pm.cpc_h_ratio), cpc_nZ=float(pm.cpc_nZ), rim_angle_x=float(pm.rim_angle_x),
rim_angle_y=pm.rim_angle_y, aim_z=float(pm.H_tower), secref_inv_eccen=pm.secref_inv_eccen, tilt_secref=float(pm.tilt_secref), rho_secref=float(pm.rho_secref),
rho_cpc=float(pm.rho_cpc), slope_error=float(pm.slope_error_bd))

bd.heliostatfield(field=pm.field_type, hst_rho=pm.rho_helio, slope=pm.slope_error, hst_w=pm.W_helio, hst_h=pm.H_helio, tower_h=pm.H_tower, tower_r=pm.R_tower,
hst_z=pm.Z_helio, num_hst=pm.n_helios, R1=pm.R1, fb=pm.fb, dsep=pm.dsep, x_max=300., y_max=300.)

bd.yaml(sunshape=pm.sunshape,csr=pm.crs,half_angle_deg=pm.half_angle_deg,std_dev=pm.std_dev)

oelt, A_land=bd.field_design_annual(dni_des=900., num_rays=int(1e5), nd=pm.n_row_oelt, nh=pm.n_col_oelt, weafile=weafile, method=1, Q_in_des=pm.Q_in_rcv, n_helios=None, zipfiles=False, gen_vtk=False, plot=False)

# Recalculate the annual efficiency of the designed field and create vtk
bd.yaml(sunshape=pm.sunshape,csr=pm.crs,half_angle_deg=pm.half_angle_deg,std_dev=pm.std_dev)

oelt, A_land=bd.annual_oelt(dni_des=900., num_rays=int(1e5), nd=pm.n_row_oelt, nh=pm.n_col_oelt, zipfiles=False, gen_vtk=True, plot=False, verbose=True)

if (A_land==0):
    tablefile=None
else:
    A_helio=pm.H_helio*pm.W_helio
    output_matadata_motab(table=oelt, field_type=pm.field_type, aiming='single', n_helios=bd.n_helios, A_helio=A_helio, eff_design=bd.eff_des, H_rcv=pm.H_rcv, W_rcv=pm.W_rcv, H_tower=pm.H_tower, Q_in_rcv=bd.Q_in_rcv, A_land=A_land, savedir=tablefile)

filename = glob(os.path.join(casefolder,'sunpos_*'))

# Read vtk and produce 1D flux map
dataname='Front_faces_Absorbed_flux'
read_vtk_annual(vtkfile=casefolder, vtkname='receiver', savedir=casefolder,  dataname=dataname, ncases=len(filename), gencsv=True)
