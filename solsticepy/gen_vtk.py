import numpy as np
import re

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


def read_vtk(vtkfile):

    '''

    '''
    f=open(vtkfile, 'r')
    content=f.readlines()
    f.close()

    l=len(content)
    i=0

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
            i+=num_points

        elif 'POLYGONS' in line:
            v= [int(s) for s in line.split() if s.isdigit()]   
            num_polygon=v[0]
            start=i+1
            j=start
            end=i+num_polygon+1
            poly=np.array([])
            while j<end:
                line=content[j]
                v =[float(s) for s in re.findall(r'-?\d+\.?\d*', line)]
                poly=np.append(poly, v)
                j+=1 
            i+=num_polygon 
        else:
            i+=1
    points=points.reshape(num_points, 3)
    poly=poly.reshape(num_polygon, 4)

    return points, poly


def flux_reader(vtkfile, casedir, check=False):

	'''
	vtkfile: str, the directory of the target vtk file
	'''
	scinot = re.compile('[+\-]?(?:0|[1-9]\d*)(?:\.\d*)?(?:[eE][+\-]?\d+)')
	f=open(vtkfile, 'r')
	content=f.readlines()
	f.close()

	l=len(content)
	i=0
	print("total lines", l)
	POINTS=np.array([])
	POLYGONS=np.array([])
	FLUX_IN=np.array([])
	FLUX_ABS=np.array([])
	FLUX_IN_back=np.array([])
	FLUX_ABS_back=np.array([])
	while i<l:
		line=content[i]
		if 'POINTS' in line:
			v= [int(s) for s in line.split() if s.isdigit()]
			num_points=v[0]
			start=i+1
			j=start
			end=i+num_points+1
			while j<end:
				line=content[j]
				v =[float(s) for s in re.findall(r'-?\d+\.?\d*', line)]
				POINTS=np.append(POINTS, v)
			
				j+=1
			print('num pos', num_points)
			print(len(POINTS)/3)
			#print(POINTS[0], POINTS[1], POINTS[2])
			POINTS=POINTS.reshape(num_points, 3)
			i+=num_points

		elif 'POLYGONS' in line:
			v= [int(s) for s in line.split() if s.isdigit()]   
			num_polygon=v[0]
			start=i+1
			j=start
			end=i+num_polygon+1
			while j<end:
				line=content[j]
				v =[float(s) for s in re.findall(r'-?\d+\.?\d*', line)]
	
				POLYGONS=np.append(POLYGONS, v[1:])
				j+=1			
			POLYGONS=POLYGONS.reshape(num_polygon, 3)
			POLYGONGS=POLYGONS.astype(int)
			i+=num_polygon 

		elif 'CELL_DATA' in line:
			v= [int(s) for s in line.split() if s.isdigit()] 
			num_data=v[0]
			print('num cell data', num_data)  
			i+=1

		elif 'Front_faces_Incoming_flux' in line:
			v = [int(s) for s in line.split() if s.isdigit()]   
			start=i+2
			j=start
			end=i+num_data+2
			while j<end:
				line=content[j]
				v =line.split(' ')# [float(s) for s in re.findall(scinot, line)]
				FLUX_IN=np.append(FLUX_IN, float(v[0])/1000.) #kW/m2
				j+=1			
			i+=num_data+1 

		elif 'Front_faces_Absorbed_flux' in line:
			v = [int(s) for s in line.split() if s.isdigit()]   
			start=i+2
			j=start
			end=i+num_data+2
			while j<end:
				line=content[j]
				v =line.split(' ')#[float(s) for s in re.findall(r'-?\d+\.?\d*', line)]
				FLUX_ABS=np.append(FLUX_ABS, float(v[0])/1000.) #kW/m2
				j+=1			
			i+=num_data+1 

		elif 'Back_faces_Incoming_flux' in line:
			v = [int(s) for s in line.split() if s.isdigit()]   
			start=i+2
			j=start
			end=i+num_data+2
			while j<end:
				line=content[j]
				v =line.split(' ')#[float(s) for s in re.findall(r'-?\d+\.?\d*', line)]
				FLUX_IN_back=np.append(FLUX_IN_back, float(v[0])/1000.) #kW/m2
				j+=1			
			i+=num_data+1 

		elif 'Back_faces_Absorbed_flux' in line:
			v = [int(s) for s in line.split() if s.isdigit()]   
			start=i+2
			j=start
			end=i+num_data+2
			while j<end:
				line=content[j]
				v =line.split(' ')#[float(s) for s in re.findall(r'-?\d+\.?\d*', line)]
				FLUX_ABS_back=np.append(FLUX_ABS_back, float(v[0])/1000.) #kW/m2
				j+=1			
			i+=num_data+1 

		else:
			i+=1

	if check:
		import matplotlib.pyplot as plt
		from mpl_toolkits.mplot3d import Axes3D
		from mpl_toolkits.mplot3d.art3d import Poly3DCollection, Line3DCollection
		from matplotlib import cm	
		import matplotlib as mpl	
		import pylab as pl	

		X=POINTS[:,0]
		Y=POINTS[:,1]
		Z=POINTS[:,2]

		flux=np.arange(len(POLYGONS))
		tri=POLYGONS.astype(int)
		fig = pl.figure()  
		ax = fig.add_subplot(111, projection = '3d')

		norm=mpl.colors.Normalize(vmin = np.min(flux), vmax =  np.max(flux))
		m = plt.cm.ScalarMappable(norm=norm, cmap='jet')
		m.set_array([])
		fcolors=m.to_rgba(flux)

		verts=POINTS[:, :3][tri]

		# the whole geomentry
		ax.plot_trisurf(X, Y, Z, triangles=tri, color=(0,0,0,0),linewidths=0.1, edgecolors='white')       
		# the visible part     
		ax.add_collection3d(Poly3DCollection(verts, 
		facecolors=fcolors, linewidths=0., edgecolors='r'))
		cbar=plt.colorbar(m)          
		plt.show()  


	return POINTS, POLYGONS, FLUX_IN, FLUX_ABS, FLUX_IN_back, FLUX_ABS_back



    

if __name__=='__main__':
    gen_vtk()

