import os, re
import numpy as np


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


def read_vtk(vtkfile, dataname=None, savecsv=False, savedir='.'):
    '''This function reads a vtkfile (.vtk), obtains the points and indices of the mesh, and the data in each mesh according to the dataname that is given,
    and saves the information in a .csv format with a name that is the same as the prefix of the .vtk file

    ``Arguments``

    * vtkfile (str): directory to the .vtk file
    * dataname (str): the name of the data that is of interest
    * savecsv (bool): Wether to save the csv files or not in the `savedir`
    * savedir (str): directory to save the .csv file

    ``Returns``
    * points (3xfloat): coordinates of each points of the mesh (vertex of each polygon)
    * indices (n*integer): n index of n vertex of each polygon
    * data (float): flux of individual cells in W/m2

    '''

    f=open(vtkfile, 'r')
    content=f.readlines()
    f.close()

    l=len(content)
    i=0
    num_polygon=0
    normals=np.array([])
    while i<l:
    	line=content[i]
    	if 'POINTS' in line:
    		v= [int(s) for s in line.split() if s.isdigit()]
    		num_points=v[0]
    		start=i+1
    		j=start
    		end=i+num_points+1
    		points=np.array([])
    		while j<end:
    			line=content[j]
    			v =[float(s) for s in re.findall(r'-?\d+\.?\d*', line)]
    			points=np.append(points, v)
    			j+=1

    		points=points.reshape(num_points, 3)
    		title=np.array(['X', 'Y', 'Z'])
    		points=np.vstack((title, points))
    		i+=num_points

    	elif 'POLYGONS' in line:
    		v= [int(s) for s in line.split() if s.isdigit()]
    		num_polygon=v[0]
    		num_verts=int(v[1]/v[0])-1 # number of vertices
    		start=i+1
    		j=start
    		end=i+num_polygon+1
    		indices=np.array([])
    		while j<end:
    			line=content[j]
    			v =[float(s) for s in re.findall(r'-?\d+\.?\d*', line)]
    			indices=np.append(indices, v[1:])
    			j+=1

    		title=np.array([])
    		for t in range(num_verts):
    			title=np.append(title, 'polygon index %s'%(t+1))
    		indices=indices.reshape(num_polygon, num_verts)
    		indices=np.vstack((title, indices))
    		i+=num_polygon

    	elif dataname in line:
    		#print('\nData: %s\n'%dataname)
    		v= [int(s) for s in re.findall(r'-?\d+\.?\d*', line)]
    		name=v[0]
    		num_v=num_polygon
    		start=i+2
    		j=start
    		end=start+num_v
    		data=np.array([dataname])
    		while j<end:
    			line=content[j]
    			v =[float(s) for s in re.findall(r'-?\S+\.?\d*', line)]
    			data=np.append(data, v[0]) #save the nominal value, v[1] is the error value
    			j+=1

    		i+=num_v
    		data=data.reshape(num_polygon+1, 1)

    	elif 'NORMALS cell_normals float' in line:
    		num_verts=3 # number of vertices
    		start=i+1
    		j=start
    		end=i+num_polygon+1
    		while j<end:
    			line=content[j]
    			v =[float(s) for s in re.findall(r'-?\d+\.?\d*', line)]
    			normals=np.append(normals, v)
    			j+=1

    		title=np.array(['X', 'Y', 'Z'])
    		normals=normals.reshape(num_polygon, 3)
    		normals=np.vstack((title, normals))
    		i+=num_polygon

    	else:
    		i+=1

    output=np.hstack((indices, data))

    if savecsv:
        vtkname = os.path.splitext(os.path.split(vtkfile)[1])[0]
        np.savetxt(savedir+'/%s_mesh_data_'%vtkname+'.csv', output, fmt='%s', delimiter=',')
        np.savetxt(savedir+'/%s_mesh_vertices_'%vtkname+'.csv', points, fmt='%s', delimiter=',')
        if len(normals)>0:
            np.savetxt(savedir+'/%s_mesh_normals_'%vtkname+'.csv', normals, fmt='%s', delimiter=',')

    return points[1:].astype(float), indices[1:].astype(float), data[1:].astype(float)


if __name__=='__main__':
    #gen_vtk()
    read_vtk(vtkfile='../tests/data/example_rec.vtk', dataname='Front_faces_Absorbed_flux', savecsv=True, savedir='../tests/data')
