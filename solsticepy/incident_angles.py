import numpy as np
from solsticepy.cal_sun import SunPosition
import matplotlib.pyplot as plt
import os

def get_incident_angle(sun_ele, sun_azi, hst_pos, hst_aims, slant, one_heliostat=False):

    sun_x=np.cos(sun_ele*np.pi/180.)*np.cos(sun_azi*np.pi/180.)
    sun_y=np.cos(sun_ele*np.pi/180.)*np.sin(sun_azi*np.pi/180.)
    sun_z=np.sin(sun_ele*np.pi/180.)
    sun_vec=np.r_[sun_x, sun_y, sun_z]
    sun_vec=sun_vec / np.linalg.norm(sun_vec)
    sun_vec=sun_vec.reshape(1, 3)

    OH=hst_aims-hst_pos # heliostat and target vectors

    if one_heliostat:
        norms = np.linalg.norm(OH) #, axis=0, keepdims=True) 
    else:
        norms = np.linalg.norm(OH, axis=1, keepdims=True) 
    unit_OH =OH/norms

    dot_products = np.sum(unit_OH * sun_vec, axis=1)
    cos_theta = dot_products
    angles_rad = np.arccos(np.clip(cos_theta, -1.0, 1.0))
    angles=angles_rad/2.

    Fx=slant*np.cos(angles)
    Fy=slant/np.cos(angles)
    return Fx, Fy, angles   


def get_clear_sky_DNI(elevation): 

    if np.sin(elevation*np.pi/180.)>=1.e-5:
        dni=1618.*np.exp(-0.606/(np.sin(elevation*np.pi/180.)**0.491))
    else:
        dni=0.
    return dni

def get_annual_incident(hst_pos, hst_aims, latitude, casename='', verbose=False):

    summary=np.array([])
    title=np.array(['sun ele','sun azi','Fx','Fy','incident angle','dni'])
    title=title.reshape(1,6)
    sun=SunPosition()
    inc_avg=np.zeros((len(hst_pos))) # DNI weighted average incident angles
    i=0
    for day in range(365):
        for h in range(24):
            delta=sun.declination(day)
            omega=-180+h*15
            theta=sun.zenith(latitude, delta, omega)
            phi=sun.azimuth(latitude, theta, delta, omega)
            sun_azi, sun_ele=sun.convert_convention(tool='solstice', azimuth=phi, zenith=theta)
            sun_azi=sun_azi
            slant=np.linalg.norm(hst_aims-hst_pos)
            Fx, Fy, angles=get_incident_angle(sun_ele, sun_azi, hst_pos, hst_aims, slant, one_heliostat=False)
            dni=get_clear_sky_DNI(sun_ele)
            inc_avg=(inc_avg+angles*dni)

            summary=np.append(summary, (sun_ele, sun_azi, dni))    
            i+=1
    summary=summary.reshape(int(len(summary)/3), 3)
    #incident_angles=summary[:,4]
    DNI=summary[:,2]
    inc_avg=inc_avg/np.sum(DNI)*180./np.pi

    return inc_avg

def master_angles():
    poses=np.arange(0, 360., 45.) # heliostat angular position
    rec_down=np.arange(15, 90, 15) # receiver lookdown angle
    H_tower=60 
    distance=H_tower*np.tan(rec_down*np.pi/180.)
    hst_aims=np.r_[0., 0., 60]
    hst_aims=hst_aims.reshape(1, 3)
    latitude=34. +53./60.
    incident_annual=np.array([])
    for d in distance:
        for p in poses:
            hst_pos=np.r_[d*np.sin(p*np.pi/180.), d*np.cos(p*np.pi/180.), 0]
            hst_pos=hst_pos.reshape(1,3)
            inc_avg=get_annual_incident(hst_pos, hst_aims,  latitude, casename='d=%.0f_p=%.0f'%(d, p), verbose=False)
            incident_annual=np.append(incident_annual, inc_avg)#np.average(incident_angles)) #
            print(d, p, hst_pos,inc_avg)

    incident_annual=incident_annual.reshape(len(distance), len(poses))
    r_edges=rec_down-(rec_down[1]-rec_down[0])/2.
    r_edges=np.append(r_edges, rec_down[-1]+(rec_down[1]-rec_down[0])/2.)
    a_edges=np.append(poses-45/2., poses[-1]+45./2.) #np.arange(-22.5, 382.5, 45)
   
    if 1: 
        fig, ax = plt.subplots(subplot_kw={'projection': 'polar'})
        ax.set_theta_zero_location('N')
        ax.set_theta_direction(-1)
        c = ax.pcolormesh(a_edges*np.pi/180., r_edges, incident_annual, cmap='viridis', shading='auto')
        cb=plt.colorbar(c, ax=ax)
        # Add colorbar with label and fontsize
   
        cb.set_label('Annualy weighted incident angle (°)', fontsize=14)
        cb.ax.tick_params(labelsize=14)
        # Add colorbar with label and fontsize
        # Set title font size
        #ax.set_title('2D Polar Colormap Example', fontsize=14)

        # Set tick label font size for both axes
        ax.tick_params(axis='both', labelsize=14)

        # Set radial label font size
        ax.set_rlabel_position(0)  # Optional: set position of radial labels
        for label in ax.get_yticklabels():
            label.set_fontsize(14)

        # Set theta label font size (angle labels)
        for label in ax.get_xticklabels():
            label.set_fontsize(14)
        #plt.show()
        plt.savefig(open('./annual_incident.png', 'wb'), bbox_inches='tight', dpi=300)
        plt.close()

    if 0:
        # Create polar plot
        fig, ax = plt.subplots(subplot_kw={'projection': 'polar'})
        xx, yy=np.meshgrid(distance, poses*np.pi/180.)
        # Scatter plot in polar coordinates with colormap
        sc = ax.scatter(yy,xx, c=incident_annual, cmap='viridis')
        # Add colorbar
        plt.colorbar(sc, ax=ax, label='Value')
        ax.set_theta_zero_location('N')
        ax.set_theta_direction(-1)
        plt.title('Polar Colormap Example')
        plt.show()

        



if __name__=='__main__':
    master_angles() 
