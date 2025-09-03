import numpy as np
from cal_sun import SunPosition

def get_incident_angle(sun_ele, sun_azi, hst_pos, hst_aims, slant, one_heliostat=False):

    sun_x=np.cos(sun_ele)*np.cos(sun_azi)
    sun_y=np.cos(sun_ele)*np.sin(sun_azi)
    sun_z=np.sin(sun_ele)
    sun_vec=np.r_[sun_x, sun_y, sun_z]
    sun_vec=sun_vec / np.linalg.norm(sun_vec)
    sun_vec=sun_vec.reshape(1, 3)

    OH=hst_aims-hst_pos # heliostat and target vectors
    if one_heliostat:
        norms = np.linalg.norm(OH) #, axis=0, keepdims=True) 
    else:
        norms = np.linalg.norm(OH, axis=0, keepdims=True) 
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


def master_angles():
    poses=np.arange(45., 360., 45.) # heliostat angular position
    rec_down=np.arange(15, 90, 15) # receiver lookdown angle
    H_tower=60 
    distance=H_tower*np.tan(rec_down*np.pi/180.)
    hst_aims=np.r_[0., 0., 60]
    hst_aims=hst_aims.reshape(1, 3)
    sun=SunPosition()
    latitude=34. +53./60.
    for d in distance:
        hst_pos=np.r_[0., d, 0]
        hst_pos=hst_pos.reshape(1,3)
        slant=np.linalg.norm(hst_aims-hst_pos)
        for p in poses:
            summary=np.array([])
            for day in range(365):
                for h in range(24):
                    delta=sun.declination(day)
                    omega=-180+h*15
                    theta=sun.zenith(latitude, delta, omega)
                    phi=sun.azimuth(latitude, theta, delta, omega)
                    sun_azi, sun_ele=sun.convert_convention(tool='solstice', azimuth=phi, zenith=theta)
                    sun_azi=sun_azi+p
                    if sun_azi>=360:
                        sun_azi-=360.

                    Fx, Fy, angles=get_incident_angle(sun_ele, sun_azi, hst_pos, hst_aims, slant, one_heliostat=True)
                    dni=get_clear_sky_DNI(sun_ele)
                    summary=np.append(summary, (d, p, Fx[0], Fy[0], angles[0], dni))
            summary=summary.reshape(int(len(summary)/6), 6)
            np.savetxt('./p_%.0f-d_%.0f-angles.csv'%(p, d), summary, delimiter=',')
        



if __name__=='__main__':
    master_angles() 
