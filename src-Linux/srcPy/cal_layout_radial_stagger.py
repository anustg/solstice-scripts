import numpy as N
from cal_field import *
from cal_sun import *
from gen_vtk import gen_vtk

def radial_stagger(Target_A, R1, width, height, towerheight, hst_z, receiver_norm, latitude, az_rim=2*N.pi, dsep=0., savedir='.'):
    '''
    Ref. (Collado and Guallar, 2012), Campo: Generation of regular heliostat field.

    Generate a rather large field

    Arguements:

    R1: distance from the first row to the bottom of the tower (0, 0, 0)
    width, height - float, heliostat dimension (m)
    az_rim - float, (rad), estimated shape
    dsep - float, separation distance
    savedir - directory of saving the pos_and_aiming.csv

    Return:
    pos_and_aiming: position, focul length and aiming points.
    
    '''

        
    # heliostat diagonal distantce
    DH=N.sqrt(height**2+width**2) 

    # distance between contiguous helistat center on the X and Y plane
    DM=DH+dsep

    # minimum radial increment
    delta_Rmin=0.866*DM

    # number of heliostats in the first row
    Nhel1 =int(2.*N.pi*R1/DM)

    target_num=Target_A/width/height

    # the total number of zones (estimated)
    az_rim_lim=N.pi/6.
    Nzones=int(N.log(5.44*3*(Target_A/az_rim_lim*N.pi/width/height)/Nhel1**2+1)/N.log(4))+1

    X=N.array([])
    Y=N.array([])

    for i in xrange(Nzones):
        
        R=(2.**(i))*(DM)*Nhel1/(2.*N.pi)

        Nrows= int((2.**(i))*Nhel1/5.44)

        Nhel=(2**(i))*Nhel1

        print ''
        print 'Zone :', i+1
        print 'R    :', R
        print 'N row:', Nrows
        print 'N hel:', Nhel

        delta_az=DM/R
        half_Nhel=Nhel

        for row in xrange(Nrows):

            r=R+float(row)*delta_Rmin

            for nh in xrange(half_Nhel):
          
                if row%2==0:
                    # the odd row
                    azimuth=delta_az/2.+float(nh)*delta_az
                else:
                    # the even row
                    azimuth=float(nh)*delta_az  
         
                xx=r*N.sin(azimuth)
                yy=r*N.cos(azimuth)   
                X=N.append(X, xx)
                Y=N.append(Y, yy)

    num=len(X)
    hstpos=N.zeros(num*3).reshape(num, 3)
    hstpos[:, 0]=X
    hstpos[:, 1]=Y
    hstpos[:,2]=hst_z

 

    sun=SunPosition()
    sun_azi, sun_zenith=sun.annual_angles(latitude, hemisphere='North', nd=5, nh=5)
    

    pf=FieldPF(sun_azi, sun_zenith, receiver_norm)            

    vis_idx=pf.get_rec_view(towerheight, hstpos)
    vis_hst=hstpos[vis_idx, :]

    cosf=pf.get_cosine(towerheight, vis_hst)
    sort_idx=N.argsort(cosf)
    sort_targ=sort_idx[len(sort_idx)-(int(target_num)+1):]

    hst=vis_hst[sort_targ]
    total_real=len(hst)

    total_area=width*height*float(total_real)
    print ''
    print 'Total number:', total_real
    print 'Total area  :', total_area/1.e6
    print 'Total design area:', Target_A/1.e6
    print 'diff', (total_area-Target_A)/Target_A

    Xv=hst[:,0]
    Yv=hst[:,1]
    Zv=N.ones(len(Xv))*hst_z
    
    #plt.figure(1)
    #plt.plot(X, Y, '.')
    #plt.scatter(Xv, Yv, c='r')
    #plt.show()

    aim_x=N.zeros(len(Xv))
    aim_y=N.zeros(len(Xv))
    aim_z=N.ones(len(Xv))*towerheight

    foc=N.sqrt((Xv-aim_x)**2+(Yv-aim_y)**2+(Zv-aim_z)**2)

    pos_and_aiming=N.append(Xv, (Yv,Zv, foc, aim_x, aim_y, aim_z))
    title=N.array(['x', 'y', 'z', 'foc', 'aim x', 'aim y', 'aim z', 'm', 'm', 'm', 'm', 'm', 'm', 'm'])
    pos_and_aiming=pos_and_aiming.reshape(7,len(pos_and_aiming)/7)
    pos_and_aiming=N.append(title, pos_and_aiming.T)
    pos_and_aiming=pos_and_aiming.reshape(len(pos_and_aiming)/7, 7)
    
    N.savetxt('%s/pos_and_aiming.csv'%savedir, pos_and_aiming, fmt='%s', delimiter=',')

    # get the view
    # TODO will be organised
    field=FieldPF(azimuth=N.r_[0.], zenith=N.r_[12.], receiver_norm=receiver_norm)
    sun_vec=field.get_solar_vector(azimuth=N.r_[0.], zenith=N.r_[12.])
    norms=field.get_normals(towerheight=towerheight, hstpos=hstpos, sun_vec=sun_vec)
    COORD, TRI, ele, nc=field.view_heliostats(width, height, norms, hstpos)
    cos=field.get_cosine(towerheight=towerheight, hstpos=hstpos)
    COS=N.repeat(cos, ele)
    DATA={'cos':COS}
    NORMS=N.repeat(norms, ele, axis=0)
    gen_vtk(savedir+'/pre_field_large.vtk', COORD.T, TRI, NORMS, True, DATA)


    norms=field.get_normals(towerheight=towerheight, hstpos=vis_hst, sun_vec=sun_vec)
    COORD, TRI, ele, nc=field.view_heliostats(width, height, norms, vis_hst)
    cos=field.get_cosine(towerheight=towerheight, hstpos=vis_hst)
    COS=N.repeat(cos, ele)
    DATA={'cos':COS}
    NORMS=N.repeat(norms, ele, axis=0)
    gen_vtk(savedir+'/pre_field_rec_vis.vtk', COORD.T, TRI, NORMS, True, DATA)

    pos=pos_and_aiming[2:,:3].astype(float)
    norms=field.get_normals(towerheight=towerheight, hstpos=pos, sun_vec=sun_vec)
    COORD, TRI, ele, nc=field.view_heliostats(width, height, norms, pos)
    cos=field.get_cosine(towerheight=towerheight, hstpos=pos)
    COS=N.repeat(cos, ele)
    DATA={'cos':COS}
    NORMS=N.repeat(norms, ele, axis=0)
    gen_vtk(savedir+'/pre_field_design.vtk', COORD.T, TRI, NORMS, True, DATA)


    return pos_and_aiming

if __name__=='__main__':
    ## a PS10 like field
    design_area=624.*1.*12.926*9.6 # m2
    R1=40.
    hst_w=13.
    hst_h=10.
    tower_h=110.

    # receiver orientation
    # receiver is mounted at the tower height, (0, 0, tower_h)
    # face to the North
    receiver_normal=N.r_[0., 1., 0.] 

    print 'Target field area:', '%.2f'%design_area, 'm2'
    print 
    pos_and_aim=radial_stagger(Target_A=design_area, R1=R1, width=hst_w, height=hst_h, towerheight=tower_h, hst_z=5., receiver_norm=receiver_normal, latitude=37.44, az_rim=2*N.pi, dsep=0.)

    


    








        
