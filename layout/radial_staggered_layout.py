import numpy as N

def radial_stagger(Nhel1, Nzones, width, height, towerheight, hst_z, az_rim=2*N.pi, dsep=0.):

    '''
    Ref. (Collado and Guallar, 2012), Campo: Generation of regular heliostat field.

    Arguements:
    Nhel1 - int, number of heliostat in the first zone
    Nzones - int, number of zones
    width, height - float, heliostat dimension (m)
    az_rim - float, (rad), start from y (North) clockwise, for trimming the surrounding field to a polar field 
    dsep - float, separation distance

    Return:
    pos_and_aiming: position, focul length and aiming points.
    
    '''

    # heliostat diagonal distantce
    DH=N.sqrt(height**2+width**2) 

    # distance between contiguous helistat center on the X and Y plane
    DM=DH+dsep

    # minimum radial increment
    delta_Rmin=0.866*DM

    total=0.
    total_real=0.
    X=N.array([])
    Y=N.array([])

    for i in xrange(Nzones):

        R=(2.**(i))*(DM)*Nhel1/(2.*N.pi)

        Nrows= int((2.**(i))*Nhel1/5.44)

        Nhel=(2**(i))*Nhel1
        total+=Nhel*Nrows
        print ''
        print 'Zone :', i+1
        print 'R    :', R
        print 'N row:', Nrows
        print 'N hel:', Nhel

        delta_az=DM/R
        half_Nhel=Nhel
        X_zone=N.array([])
        Y_zone=N.array([])

        for row in xrange(Nrows):

            r=R+float(row)*delta_Rmin

            for nh in xrange(half_Nhel):
          
                if row%2==0:
                    # the odd row
                    azimuth=delta_az/2.+float(nh)*delta_az
                else:
                    # the even row
                    azimuth=float(nh)*delta_az  

                if (azimuth <= az_rim) or (azimuth >=2.*N.pi-az_rim):                  

                    xx=r*N.sin(azimuth)
                    yy=r*N.cos(azimuth)   
                    X=N.append(X, xx)
                    Y=N.append(Y, yy)

                    X_zone=N.append(X_zone, xx)
                    Y_zone=N.append(Y_zone, yy)
                    total_real+=1

    total_area=width*height*float(total_real)
    print ''
    print 'Total number:', total_real
    print 'Total area  :', total_area/1.e6

    Z=N.ones(len(X))*hst_z

    aim_x=N.zeros(len(X))
    aim_y=N.zeros(len(X))
    aim_z=N.ones(len(X))*towerheight
    foc=N.sqrt((X-aim_x)**2+(Y-aim_y)**2+(Z-aim_z)**2)

    pos_and_aiming=N.append(X, (Y,Z, foc, aim_x, aim_y, aim_z))
    title=N.array(['x', 'y', 'z', 'foc', 'aim x', 'aim y', 'aim z', 'm', 'm', 'm', 'm', 'm', 'm', 'm'])
    pos_and_aiming=N.append(title, pos_and_aiming)
    pos_and_aiming=pos_and_aiming.reshape(7, len(pos_and_aiming)/7)
    N.savetxt('./pos_and_aiming.csv', pos_and_aiming.T, fmt='%s', delimiter=',')
    return pos_and_aiming

if __name__=='__main__':
    radial_stagger(Nhel1=52, Nzones=4, width=10., height=10., towerheight=200., hst_z=5., az_rim=N.pi*2., dsep=0.)
    
    #radial_stagger(Nhel1=35, Nzones=3, width=12.3, height=9.75, towerheight=200., hst_z=5., az_rim=N.pi*2., dsep=0.)
    








        
