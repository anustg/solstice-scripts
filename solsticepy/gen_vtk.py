
def gen_vtk(savedir, points, indices, norms, colormap=True, DATA=None):
    '''Generate 3D views of a heliostat field with triangular mesh in the VTK format that can be visualised in ParaView software
    
    ``Arguments``
    
      * savedir (str): directory to save the VTK file
      * points (nx3 numpy array): the vertices of the objects, each column of the array is X, Y, Z coordinates respectively
      * indices (nx3 numpy array): the indices of the triangular mesh
      * norms (nx3 numpy array): the normal vectors of the triangular mesh
      * colormap (bool): True - show the data of each heliostat (e.g. cosine factor), False - not show the data results, but only show the geometry of the heliostats
      * DATA (dic): key is 'cosine' or 'atm' or others performance parameter to be visualised

    ``Return``

      * No return value (a VTK file is created and written in the `savedir`)    
    '''
    num_points=len(points.T)
    num_tri=len(indices)
    f=open(savedir, 'w')
    f.write('# vtk DataFile Version 2.0\n')
    f.write('test\n')
    f.write('ASCII\n')
    f.write('DATASET POLYDATA\n')
    f.write('POINTS %s double\n'%num_points)
    for i in range(num_points):
        x=points[0,i]
        y=points[1,i]
        z=points[2,i]
        f.write('%.8f %.8f %.8f\n'%(x, y, z))
    f.write('POLYGONS %s %s\n'%(num_tri, num_tri*4))
    for i in range(num_tri):
        id1=indices[i,0]
        id2=indices[i,1]
        id3=indices[i,2]
        f.write('3 %.0f %.0f %.0f\n'%(id1, id2, id3))  
    f.write('CELL_DATA %s\n'%(num_tri))
    f.write('NORMALS cell_normals float\n')

    for i in range(num_tri):
        nx=norms[i,0]
        ny=norms[i,1]
        nz=norms[i,2]
        
        f.write('%.8f %.8f %.8f\n'%(nx, ny, nz))  
        

    if colormap:    

        totalnum=len(DATA)         
        f.write('FIELD PrimaryData %s\n'%(totalnum))
        interest=DATA.keys()
        for m in interest:
            f.write('%s 1 %s double\n'%(m, num_tri))
            for i in range(num_tri):
                f.write('%.8f\n'%(DATA[m][i]))            

       
    f.close()


if __name__=='__main__':
    gen_vtk()

