#! /bin/env python3

from __future__ import division
import unittest

from solsticepy.cal_sun import SunPosition
from solsticepy.gen_yaml import Sun, gen_yaml, gen_yaml_target_aligned
from solsticepy.master import Master
from solsticepy.cal_layout import radial_stagger
from solsticepy.gen_vtk import flux_reader, read_vtk
import os, sys
import numpy as np
import matplotlib.pyplot as plt
import colorama
from solsticepy.master import yellow, green, red


def test_alignment():
    hst_w=12.2
    hst_h=12.2
    
    sunshape='pillbox'
    sunshape_param=1e-4
    slope_error=1e-6

    target_aligned=True
    if target_aligned:
        casename='test-target-aligned-tracking-quadricxy'
    else:
        casename='test-azi-ele-tracking'


    hst_z=hst_w*0.7
    hst_pos=np.r_[0, 500, hst_z]
    hst_pos=hst_pos.reshape(1,3)

    #==================================================
    # INPUT PARAMETERS

    # the sun
    # =========
    DNI = 1000 # W/m2
    sun = Sun(dni=DNI, sunshape=sunshape, half_angle_deg=sunshape_param, std_dev=sunshape_param)

    latitude=34. +53./60.   # latitude of the crs plant
    # S4. number of rays for the ray-tracing simulation
    num_rays=int(20e6)
    #
    # the field
    # ==========
    rho_refl=1. # mirror reflectivity
    # F3. Tower
    tower_h=0.01 # tower height
    tower_r=0.01 # tower radius
    #
    # the receiver
    # ============
    # R1. shape
    receiver='flat' # 'flat' or 'stl'
    rec_w=20 # width, m
    rec_h=20 # height, m
    tilt=0.  # deg
    loc_x=0. # m
    loc_y=0. # m
    loc_z=62.# m
    rec_abs=1.
    rec_mesh=100
                                          
    casefolder=casename
    shape='flat'
    rec_param=np.r_[rec_w, rec_h, rec_mesh, rec_mesh, loc_x, loc_y, loc_z, tilt]
    hst_aims=np.r_[loc_x, loc_y, loc_z]
    hst_aims=hst_aims.reshape(1, 3)
    hst_foc=np.linalg.norm(hst_aims-hst_pos)

    master=Master(casedir=casefolder)
    outfile_yaml = master.in_case(folder=casefolder, fn='input.yaml')
    outfile_recv = master.in_case(folder=casefolder, fn='input-rcv.yaml')


    data=np.array(['pillbox', 'slope (mrad)', 'hst size (m)', 'position', 'azimuth', 'elevation', 'DNI (W/m2)', 'peak flux (W/m2)', 'beam area (m2)', 'aligned', 'H-incident zenith (deg)','H-incident azimuth (deg)'])
    poses=np.arange(45., 360., 45.)
    nd=5
    nh=22
    for p in poses:

        casedir=casefolder+'/P_%.0f'%(p)

        sunpos=SunPosition()
        AZI, ZENITH,table,case_list=sunpos.annual_angles(latitude, casefolder=casedir, nd=nd, nh=nh)
        case_list=case_list[1:]
        SOLSTICE_AZI, SOLSTICE_ELE=sunpos.convert_convention('solstice', AZI, ZENITH)
        run=np.r_[0]
        for i in range(len(case_list)):     
            c=int(case_list[i,0].astype(float))
            if c not in run:
                azimuth=SOLSTICE_AZI[c-1]
                elevation=SOLSTICE_ELE[c-1]
                if np.sin(elevation*np.pi/180.)>=1.e-5:
    	            dni=1618.*np.exp(-0.606/(np.sin(elevation*np.pi/180.)**0.491))
                else:
    	            dni=0.

                onesunfolder=casedir+'/sunpos_%.0f'%(c)
                # run Solstice using the generate inputs, and run all required post-processing
                azimuth=azimuth+p
                if azimuth>=360:
                    azimuth-=360.

                azimuth=int(azimuth)
                elevation=int(elevation)
                hst_vtk_fn=onesunfolder+'/%.0f-%.0f-primaries.vtk'%(azimuth, elevation)

                sun_azi=azimuth/180.*np.pi
                sun_ele=elevation/180.*np.pi
                # generate the YAML file from the input parameters specified above
                gen_yaml_target_aligned(sun, hst_pos, hst_aims, hst_w, hst_h
                    ,  rho_refl, slope_error, receiver, rec_param, rec_abs
                    , outfile_yaml=outfile_yaml, outfile_recv=outfile_recv, sun_azi=sun_azi, sun_ele=sun_ele
                    , hemisphere='North', tower_h=tower_h, tower_r=tower_r, spectral=False
                    , medium=0, one_heliostat=True)


                if not os.path.exists(hst_vtk_fn):
                    print(onesunfolder, 'dni', dni)
                    print(hst_vtk_fn)
                    master.run(azimuth, elevation, num_rays, rho_refl, dni, folder=onesunfolder, gen_vtk=True,verbose=True)
        

                area, peak=plot_flux(onesunfolder, dni, name='%.0f-%.0f-target_e'%(azimuth, elevation), casename='sunpos %s'%c,loc_z=loc_z)

                testx, testy, zenith_angle_deg, azimuth_angle_deg=check_alignment(onesunfolder, azimuth, elevation, plot=False)
                data=np.append(data, (sunshape_param, slope_error*1000., hst_w, p,  azimuth, elevation, dni, peak, area, testx, zenith_angle_deg, azimuth_angle_deg))

                sys.stderr.write(yellow("Case %s, %s\n"%(p, onesunfolder))) 
                if testx:
                    sys.stderr.write(green("Targed aligned\n")) 
                else:
                    sys.stderr.write(red("Not aligned\n"))    
                print()             


        data=data.reshape(int(len(data)/12), 12)
        np.savetxt(casedir+'/data-summary.csv', data, fmt='%s', delimiter=',')


def plot_flux(casefolder, dni, name, casename, loc_z):
	vtkfile=casefolder+'/'+name+'.vtk'

	points, mesh, FLUX_IN, FLUX_ABS, FLUX_IN_back, FLUX_ABS_back=flux_reader(vtkfile, casedir=casefolder, check=False) #,  dataname='Front_faces_Absorbed_flux')

	#points=np.loadtxt(casefolder+'/%s_points.csv'%(name), skiprows=1, delimiter=',')
	#mesh=np.loadtxt(casefolder+'/%s_mesh_data.csv'%(name), skiprows=1, delimiter=',')

	flux=FLUX_ABS #mesh[:,3]*dni/1000.
	peak=np.max(flux)

	if 1:
		fts=14
		plt2,ax=plt.subplots()
		cm = plt.cm.get_cmap('jet')

		#plt.triplot(values[0], values[1], indices) 
		plt.tripcolor(points[:,0], points[:,2]-loc_z, mesh[:,:3] , facecolors=flux, cmap=cm)  	

		clb=plt.colorbar()
		clb.ax.set_title(' Flux \n W/$m^2$',y=1.,fontsize=fts)
		plt.xlabel('X',fontsize=fts)
		plt.ylabel('Y',fontsize=fts)
		plt.xticks(fontsize=fts)
		plt.yticks(fontsize=fts)
		plt.title(casename+'\nDNI=%.0f W/m$^2$'%dni)
		clb.ax.tick_params(labelsize=fts)
		#plt.legend(bbox_to_anchor=(1.05, 1), loc=2, borderaxespad=0.)
		#plt.legend(loc='lower left',fontsize=fts)
		plt.savefig(open(casefolder+'_target_e_absorbed_flux.png','wb'), bbox_inches='tight', dpi=300)
		#plt.savefig(open('pillbox_%.0f_slope_%.0f_%s_%s_%sm_target_e_absorbed_flux.png'%(sunshape_param, slope_error*1000., nst_name, time, hst_w),'wb'), bbox_inches='tight', dpi=300)

		#plt.savefig(open('pillbox_%.0f_slope_%.0f_%s_%s_%sm_target_e_absorbed_flux.png'%(sunshape_param, slope_error*1000., nst_name, time, hst_w),'wb'), bbox_inches='tight', dpi=300)
		#plt.show()
		plt.close()

	flux=np.sort(flux)
	flux=flux[::-1]
	Qtot=np.sum(flux)
	i=0
	ne=0
	Q=0.
	frac=0.
	
	while frac<0.9:
		q=flux[i]			
		Q+=q
		frac=Q/Qtot
		i+=1
		ne+=1

	ntot=len(flux)
	area=float(ne)/float(ntot)*20.*20. # the targe size is 20 by 20 m
	return area, peak

def calculate_normal(triangle):
    # Get the vertices of the triangle
    p1, p2, p3 = triangle
    # Calculate two edges of the triangle
    v1 = p2 - p1
    v2 = p3 - p1
    # Compute the cross product to get the normal vector
    normal = np.cross(v1, v2)
    # Normalize the normal vector
    normal = normal / np.linalg.norm(normal)
    return normal

def check_alignment(casefolder, sun_azi, sun_ele, plot=False):

    hst_vtk_fn=casefolder+'/%.0f-%.0f-primaries.vtk'%(sun_azi, sun_ele)
    path_vtk_fn=casefolder+'/%.0f-%.0f-paths.vtk'%(sun_azi, sun_ele)

    # get the heliostat surface mesh
    points, poly, lines=read_vtk(hst_vtk_fn)
    check=True
    i=0
    while check and i<len(poly):
        triangle=points[poly[i,1:].astype(int)]
        normal=calculate_normal(triangle)
        normal = normal / np.linalg.norm(normal) 

        x=points[:,0]
        y=points[:,1]
        z=points[:,2]

        p0=triangle[0]
        p1=triangle[1]
        p2=triangle[2]
        # Calculate two edges of the triangle
        edge1 = p2 - p0  # First edge of the triangle
        edge2 = p2 - p1  # Second edge of the triangle
        edge3 = p1 - p0

        e1_T_e2=np.dot(edge1, edge2)/np.linalg.norm(edge1)/np.linalg.norm(edge2)
        e2_T_e3=np.dot(edge3, edge2)/np.linalg.norm(edge3)/np.linalg.norm(edge2)

        if abs(e1_T_e2)<1e-4:
            check=False
            perpend=e1_T_e2

        elif abs(e2_T_e3)<1e-4:
            check=False
            edge1=edge3
            perpend=e2_T_e3
        i+=1

    print('edge 1 perpendicular to edge 2: ', abs(perpend)<1e-4, perpend)


    # Define the local x-axis as the first edge (edge1)
    local_y = edge1 / np.linalg.norm(edge1)  
    local_x = edge2 / np.linalg.norm(edge2)  

    local_z = np.cross(local_x, local_y)
    local_z = local_z / np.linalg.norm(local_z)  # Normalize the vector
    z_T_n=np.dot(local_z, normal)/np.linalg.norm(normal)/np.linalg.norm(local_z)
    print('local z parallel to normal: ', (abs(z_T_n-1.) <1e-4 or abs(z_T_n+1.) <1e-4), z_T_n)

    # Calculate the centroid of the plane (for plotting the axes)
    centroid = np.array([np.mean(x), np.mean(y), np.mean(z)])



    if plot:
        # Create the figure and 3D axis
        fig = plt.figure()
        ax = fig.add_subplot(111, projection='3d')
        ax.plot(triangle[0,0], triangle[0,1], triangle[0,2], 'bo') #p0
        ax.plot(triangle[1,0], triangle[1,1], triangle[1,2], 'ro') #p1
        ax.plot(triangle[2,0], triangle[2,1], triangle[2,2], 'go') #p2
        # Plot the plane (triangular surface)
        ax.plot_trisurf(x, y, z, triangles=poly, color='lightblue', alpha=0.5, edgecolor='k')


    #print(path_vtk_fn)
   # get the incident and reflection rays
    points, poly, lines=read_vtk(path_vtk_fn)
    x=points[:,0]
    y=points[:,1]
    z=points[:,2]
    for l in lines[:1]:

        l=l.astype(int)
        aa=points[l[0]]
        bb=points[l[1]]
        cc=points[l[2]]
        #print(aa, bb, cc)
 
        OR=bb-cc
        OS=bb-aa
        n_SR=np.cross(OR, OS)
        cosx=np.dot(n_SR, local_x)/np.linalg.norm(n_SR)/np.linalg.norm(local_x)
        cosy=np.dot(n_SR, local_y)/np.linalg.norm(n_SR)/np.linalg.norm(local_y)

        #print()
        #print('angle n_SR and local_y', cosy)
        #print('angle n_SR and local_x', cosx)   

    zenith_angle_deg, azimuth_angle_deg	= calculate_angles(OS, local_x, local_y, local_z)


    testx=abs(cosx-1)<1e-5 or abs(cosx+1)<1e-5 
    testy=abs(cosy)<1e-5 

    print('Target aligned', testx, cosx, testy, cosy)
    print('zenith = %.0f, azimuth =%.0f'%(zenith_angle_deg, azimuth_angle_deg))

    if plot:
  
        # Plot the normal vector
        ax.quiver(
            centroid[0], centroid[1], centroid[2],  # Starting point
            local_z[0], local_z[1], local_z[2],       # Direction
            color='red', length=3, normalize=True, label='Normal Vector'
        )

        # Plot the local x-axis
        ax.quiver(
            centroid[0], centroid[1], centroid[2],  # Starting point
            local_x[0], local_x[1], local_x[2],     # Direction
            color='green', length=3, normalize=True, label='Local X-axis'
        )

        # Plot the local y-axis
        ax.quiver(
            centroid[0], centroid[1], centroid[2],  # Starting point
            local_y[0], local_y[1], local_y[2],     # Direction
            color='blue', length=3, normalize=True, label='Local Y-axis'
        )

      

        # draw a plane come through local y and z
        # Define the range for the local y and z axes
        y_range = np.linspace(-5, 5, 10)  # Range for the local y-axis
        z_range = np.linspace(-5, 5, 10)  # Range for the local z-axis
        # Create a grid of points in the local y-z plane
        y_grid, z_grid = np.meshgrid(y_range, z_range)




        # Compute the points in the plane
        plane_points = centroid + y_grid[:, :, None] * local_y + z_grid[:, :, None] *local_z

        # Extract the x, y, z coordinates of the plane
        x_plane = plane_points[:, :, 0]
        y_plane = plane_points[:, :, 1]
        z_plane = plane_points[:, :, 2]

        ax.plot_surface(x_plane, y_plane, z_plane, alpha=0.5, color='cyan')
     

        # plot the rays
        #points, poly, lines=read_vtk(path_vtk_fn)
        #x=points[:,0]
        #y=points[:,1]
        #z=points[:,2]
        for l in lines:
            l=l.astype(int)
            aa=points[l[0]]
            bb=points[l[1]]
            cc=points[l[2]]
            line_points=np.array([aa, bb, cc])
            line_centroid = bb#np.mean(line_points, axis=0)
            translation_vector = centroid - line_centroid
            translated_line_points = line_points + translation_vector
            ax.plot(translated_line_points[:, 0], translated_line_points[:, 1], translated_line_points[:, 2], color='b')


        # Customize the view
        ax.set_xlabel('X')
        ax.set_ylabel('Y')
        ax.set_zlabel('Z')
        ax.set_title('Local Axes of a Plane')

        # Set axis limits
        ax.set_xlim([-50, 50])  # Set x-axis limits
        ax.set_ylim([480, 500])  # Set y-axis limits
        ax.set_zlim([0, 20])  # Set z-axis limits

        set_axes_equal(ax)
        # Add a legend
        ax.legend()

        # Show the plot
        plt.show()
    return testx, testy, zenith_angle_deg, azimuth_angle_deg

def set_axes_equal(ax):
    """Set equal scaling for 3D plot axes."""
    x_limits = ax.get_xlim()
    y_limits = ax.get_ylim()
    z_limits = ax.get_zlim()

    x_range = abs(x_limits[1] - x_limits[0])
    y_range = abs(y_limits[1] - y_limits[0])
    z_range = abs(z_limits[1] - z_limits[0])

    max_range = max(x_range, y_range, z_range)

    x_mid = np.mean(x_limits)
    y_mid = np.mean(y_limits)
    z_mid = np.mean(z_limits)

    ax.set_xlim([x_mid - max_range / 2, x_mid + max_range / 2])
    ax.set_ylim([y_mid - max_range / 2, y_mid + max_range / 2])
    ax.set_zlim([z_mid - max_range / 2, z_mid + max_range / 2])

def calculate_angles(S, x_axis, y_axis, z_axis):

    # Calculate magnitudes
    magnitude_S = np.linalg.norm(S)
    magnitude_z = np.linalg.norm(z_axis)
    magnitude_y = np.linalg.norm(y_axis)

    # Calculate the dot product for zenith angle
    dot_product_S_z = np.dot(S, z_axis)

    # Calculate the zenith angle
    zenith_angle = np.arccos(dot_product_S_z / (magnitude_S * magnitude_z))

    # Project S onto the xy-plane
    projection_S_z = (dot_product_S_z / magnitude_z**2) * z_axis
    S_xy = S - projection_S_z

    # Calculate the magnitude of the projection
    magnitude_S_xy = np.linalg.norm(S_xy)

    # Calculate the dot product for azimuth angle
    dot_product_S_xy_y = np.dot(S_xy, y_axis)
    dot_product_S_xy_x = np.dot(S_xy, x_axis)

    # Calculate the azimuth angle using atan2 for correct quadrant
    azimuth_angle = np.arctan2(dot_product_S_xy_x, dot_product_S_xy_y)

    # Convert angles from radians to degrees
    zenith_angle_deg = np.degrees(zenith_angle)
    azimuth_angle_deg = np.degrees(azimuth_angle)

    return zenith_angle_deg, azimuth_angle_deg	

def plot_hst_azi_zenith():
    '''
    plot the incident angle relate to the heliostat azi and zenith

    '''
    cases=['test-target-aligned-tracking-quadricxy1'] #'test-target-aligned-tracking', 'test-azi-ele-tracking']
    names=['target-aligned tracking']#, 'azi-ele tracking']
    poses=np.r_[45] #np.arange(0, 360., 45.)
    '''
    for p in poses:
        # Create the polar plot
        plt.rcParams.update({'font.size': 14})
        fig, ax = plt.subplots(subplot_kw={'projection': 'polar'})
        for i,c in enumerate(cases):
            casefolder=c
            casedir=casefolder+'/P_%.0f'%(p)
            fn=casedir+'/data-summary.csv'
            data=np.loadtxt(fn, delimiter=',', skiprows=1)
            zenith=data[:,-2]
            azimuth=data[:,-1]
        
            #r=np.sin(zenith*np.pi/180.)
            azimuth=azimuth*np.pi/180.    
            ax.plot(azimuth, zenith, '.', label=names[i])

        # Add labels and legend
        ax.set_title('P=%.0f°'%(p), va='bottom')
        ax.legend(loc='upper left', bbox_to_anchor=(1.05, 1.05))

        
        # Show the plot
        #plt.show()    
        plt.savefig(open('./azi_zen_P_%.0f.png'%p, 'wb'), bbox_inches='tight', dpi=300)
        plt.close()
    '''

    # Create a 4x2 grid of subplots (8 subplots in total)
    fig, axes = plt.subplots(2, 4, subplot_kw={'projection': 'polar'}, figsize=(18, 10))
 
    for j,p in enumerate(poses):
        # Create the polar plot
        plt.rcParams.update({'font.size': 14})
        ax=axes.flat[j]

        for i,c in enumerate(cases):
            casefolder=c
            casedir=casefolder+'/P_%.0f'%(p)
            fn=casedir+'/data-summary.csv'
            data=np.loadtxt(fn, delimiter=',', skiprows=1)
            zenith=data[:,-2]
            azimuth=data[:,-1]
        
            #r=np.sin(zenith*np.pi/180.)
            azimuth=azimuth*np.pi/180.    
            ax.plot(azimuth, zenith, '.', label=names[i])

        # Add labels and legend
        ax.set_title('P=%.0f°'%(p), va='bottom')
        #ax.legend(loc='upper left', bbox_to_anchor=(1.05, 1.05))

        
    # Show the plot
    #plt.show()    
    plt.savefig(open('./azi_zen_all-quadric.png', 'wb'), bbox_inches='tight', dpi=300)
    plt.close()







if __name__=='__main__':
    test_alignment()
    #plot_hst_azi_zenith()



