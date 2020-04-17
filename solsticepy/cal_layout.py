import numpy as np
import sys
import matplotlib.pyplot as plt
from .cal_field import *
from .cal_sun import *
from .gen_vtk import gen_vtk

def radial_stagger(num_hst, width, height, hst_z, towerheight, R1, dsep=0., field='polar-half', savedir='.', plot=False):
    '''
    Ref. (Collado and Guallar, 2012), Campo: Generation of regular heliostat field.

    Generate a rather large field

    Arguements:
    num_hst    : int, number of heliostats
    width      : float, mirror width (m)
    height     : float, mirror height (m)
    hst_z      : float, the vertical location of each heliostat (m)
    towerheight: float, tower height (m)
    R1         : float, distance from the first row to the bottom of the tower (0, 0, 0)
    dsep       : float, separation distance (m)
    field      : str, 'polar-half' or 'surround-half' or 'polar' or 'surround' field, 
                       the 'half' option is for simulation a symmetric field
    savedir    : str, directory of saving the pos_and_aiming.csv

    Return:
    pos_and_aiming: (n, 7) array, position, focal length and aiming point of each generated heliostat
    '''
 
    # heliostat diagonal distantce
    DH=np.sqrt(height**2+width**2) 

    # distance between contiguous helistat center on the X and Y plane
    DM=DH+dsep

    # minimum radial increment
    delta_Rmin=0.866*DM

    # number of heliostats in the first row
    Nhel1 =int(2.*np.pi*R1/DM)


    # the total number of zones (estimated)
    #Nzones=int(np.log(5.44*3*(num_hst/az_rim*np.pi)/Nhel1**2+1)/np.log(4))+1

    X=np.array([])
    Y=np.array([])

    num=0
    i=0
    sys.stderr.write('DM '+repr(DM)+'\n')
    sys.stderr.write('dRm'+repr(delta_Rmin)+'\n')

    if field[-4:]=='half':

        while num<num_hst/2:
            Nrows= int((2.**(i))*Nhel1/5.44)
            Nhel=(2**(i))*Nhel1
            R=Nhel/2./np.pi*DM
            delta_az=2.*np.pi/Nhel
            sys.stderr.write('\nzone %d\n' % (i,))
            sys.stderr.write("R %f\n"%(R,))
            sys.stderr.write("rows %d\n"%(nrows,))
            sys.stderr.write("hst %d\n"%(Nhel,))
            sys.stderr.write("daz %f\n"%(delta_az,))

            for row in xrange(Nrows):

                r=R+float(row)*delta_Rmin

                for nh in xrange(Nhel):

                    if num<num_hst/2:   
               
                        if row%2==0:
                            # the odd row
                            azimuth=delta_az/2.+float(nh)*delta_az

                        else:
                            # the even row
                            azimuth=float(nh)*delta_az  
                 
                        xx=r*np.sin(azimuth)
                        yy=r*np.cos(azimuth)
    
                        if field=='polar-half':
                            if (xx>0 and yy>0):  
                                X=np.append(X, xx)
                                Y=np.append(Y, yy)
                                num+=1
                        elif field=='surround-half':
                            if (yy>0):  
                                X=np.append(X, xx)
                                Y=np.append(Y, yy)
                                num+=1
            i+=1

    else:
        while num<num_hst:
            Nrows= int((2.**(i))*Nhel1/5.44)
            Nhel=(2**(i))*Nhel1
            R=Nhel/2./np.pi*DM
            delta_az=2.*np.pi/Nhel
            
            sys.stderr.write("\nzone %d"%(i,))
            sys.stderr.write("R %f\n"%(R,))
            sys.stderr.write("rows %d\n"%(nrows,))
            sys.stderr.write("hst %d\n"%(Nhel,))
            sys.stderr.write("daz %f\n"%(delta_az,))		

            for row in xrange(Nrows):

                r=R+float(row)*delta_Rmin

                for nh in xrange(Nhel):
              
                    if num<num_hst:

                        if row%2==0:
                            # the odd row
                            azimuth=delta_az/2.+float(nh)*delta_az

                        else:
                            # the even row
                            azimuth=float(nh)*delta_az  
                 
                        xx=r*np.sin(azimuth)
                        yy=r*np.cos(azimuth)
             
                        if field=='polar':
                            if (yy>0):  
                                X=np.append(X, xx)
                                Y=np.append(Y, yy)
                                num+=1
                        elif field=='surround':
                                X=np.append(X, xx)
                                Y=np.append(Y, yy)
                                num+=1
            i+=1

    sys.stderr.write("\ntotal hst %d\n"%(num,))
    hstpos=np.zeros(num*3).reshape(num, 3)
    hstpos[:, 0]=X
    hstpos[:, 1]=Y
    hstpos[:,2]=hst_z

    aim_x=np.zeros(num)
    aim_y=np.zeros(num)
    aim_z=np.ones(num)*towerheight

    foc=np.sqrt((X-aim_x)**2+(Y-aim_y)**2+(hstpos[:,2]-aim_z)**2)

    pos_and_aiming=np.append(X, (Y, hstpos[:,2], foc, aim_x, aim_y, aim_z))
    title=np.array(['x', 'y', 'z', 'foc', 'aim x', 'aim y', 'aim z', 'm', 'm', 'm', 'm', 'm', 'm', 'm'])
    pos_and_aiming=pos_and_aiming.reshape(7,len(pos_and_aiming)/7)
    pos_and_aiming=np.append(title, pos_and_aiming.T)
    pos_and_aiming=pos_and_aiming.reshape(len(pos_and_aiming)/7, 7)
    
    np.savetxt('%s/pos_and_aiming.csv'%savedir, pos_and_aiming, fmt='%s', delimiter=',')

    if plot:
        fts=24
        plt.figure(dpi=100.,figsize=(12,9))
        plt.plot(X, Y, '.')
        plt.xlim(-1000, 1000)
        plt.ylim(-1000, 1000)
        plt.xticks(fontsize=fts)
        plt.yticks(fontsize=fts)
        plt.savefig(open('%s.png'%field,'w'), bbox_inches='tight')
        plt.close()

    return pos_and_aiming


if __name__=='__main__':
    
    pos_and_aim=radial_stagger(num_hst=4000, width=12.3, height=9.75, hst_z=1., towerheight=100, R1=87.5, dsep=0., field='polar', savedir='.', plot=True)

        
