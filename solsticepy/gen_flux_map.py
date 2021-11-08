import numpy as np

def basevectors(point1, point2, point3):
    '''
    Generate the base vectors and normal vector of the plane of the triangle T[point1,point2,point3]
    The base vectors in the plane are formed along two edges of the triangle if it has a right angle, else only one edge.
    '''
    vertices = np.vstack((point1,point2,point3))

    vectors = np.zeros((3,len(point1)))
    for i in range(3):
        vectors[i] = (vertices[i]-vertices[i-1]) / np.linalg.norm(vertices[i]-vertices[i-1])

    dotproducts = np.zeros((3))
    for i in range(3):
        dotproducts[i] = abs(np.vdot(vectors[i],-vectors[i-1]))

    closest_right_angle = min(dotproducts)
    i = 0
    while dotproducts[i] != closest_right_angle:
        i+=1

    origin = vertices[i-1]
    vbase_x = vectors[i-1]
    vbase_y = vectors[i]
    vbase_z = np.cross(vbase_x, vbase_y)
    vbase_z = vbase_z / np.linalg.norm(vbase_z)

    if dotproducts[i] != 0:
        vbase_y = np.cross(vbase_x, vbase_z)

    print('Cos vectors base: ', abs(np.vdot(vbase_x,vbase_y)))
    return origin, vbase_x, vbase_y, vbase_z

def projection(point3D, origin, vbase_x, vbase_y, vbase_z):
    '''
    Projection of 3D point onto 2D plane given by base and normal vectors
    '''
    x = point3D[0]
    y = point3D[1]
    z = point3D[2]
    projection_x = x * vbase_x[0] + y * vbase_y[0] + z * vbase_z[0] + origin[0]
    projection_y = x * vbase_x[1] + y * vbase_y[1] + z * vbase_z[1] + origin[1]

    return projection_x, projection_y


def trianglecentralpoint(vertices):
	'''
	Return the middle point of the 2 vertices on the longest edge of the triangle polygon (cell)
	'''
	max = 0
	index = 0
	n_vertices = len(vertices)
	for j in range(n_vertices):
		sqrt_length = np.linalg.norm(vertices[j]-vertices[j-1])
		if  sqrt_length >= max:
			max = sqrt_length
			index = j
	middle_point = (vertices[index] + vertices[index-1]) / 2

	return middle_point


def sortarray(y_coordinate, boundaries, coordinates, values, areas):
	'''
	coordinates: 1D array
	values: array of Power (W) of individual cells
	areas: array of area (m2) of individual cells
	Sum values of duplicated coordinates and sort values by ascending associated coordinate
	'''
	output_coordinates = np.array([])
	output_values = np.array([])
	i = 1
	while len(coordinates)>0:
		select_pt = (coordinates <= boundaries[i])
		not_select_pt = np.invert(select_pt)
		x_coordinate = (boundaries[i-1] + boundaries[i])/2.
		output_coordinates = np.append(output_coordinates, x_coordinate)
		flux = np.sum(values[select_pt]) / np.sum(areas[select_pt]) # Convert Power (W) in Flux (W/m2)
		output_values = np.append(output_values, flux)
		values = values[not_select_pt]
		coordinates = coordinates[not_select_pt]
		areas = areas[not_select_pt]
		i+=1

	n_array = len(output_coordinates)
	output_values = np.vstack((output_coordinates, np.full((n_array), y_coordinate), output_values))

	return output_values.T


def genfluxmap(points, indices, data):
	'''
	Gives 2D and 1D flux map
	points: 2D or 3D coordinates of the mesh (polygons)
	indices: index of polygon vertex
	data: flux of each polygon (mesh cell) in W/m2
	2D flux map (W/m2): nber of columns of pixels + list of flux distribution of each rectangle pixel (n column, m rows),
	flux distribution starts at the bottom left corner of the discretized surface
	1D flux map (W/m2): list of flux distribution of the surface discretized along its local second axes (1 column only, n rows),
	flux distribution starts at the bottom of the discretized surface
	'''
	n_dimension = len(points[0])
	n_vertices = len(indices[0])
	n_polygons = len(indices)

	# Convert 3D in 2D when surface is perpendicular to axis X, Y, or Z.
	if n_dimension > 2:
		vertices = np.zeros((n_vertices,n_dimension))
		for j in range(n_vertices):
			vertices[j] = points[int(indices[0,j].astype(float))]
		points = localcoordinates(points, vertices)
	n_dimension = len(points[0])

	# Calculate Center of each triangle mesh longest edge
	centers = np.zeros((n_polygons, n_dimension))
	cell_areas = np.zeros((n_polygons,1))
	for i in range(n_polygons):

		vertices = np.zeros((n_vertices,n_dimension))
		for j in range(n_vertices):
			vertices[j] = points[int(indices[i,j].astype(float))]
		center = trianglecentralpoint(vertices)
		for j in range(n_dimension):
			centers[i,j] = center[j]

		cell_areas[i] = trianglearea(vertices)
		# Convert the Flux (W/m2) of each cell in Power (W)
		data[i] = data[i].astype(float)*cell_areas[i]

	flux_map_2D = np.array([])
	flux_map_1D = np.array([])
	y_list = sorted(set(points[:,1].astype(float)))
	x_list = sorted(set(points[:,0].astype(float)))
	i = 1
	while len(centers)>0:
		yy = centers[:,1].astype(float)
		select_pt = (yy <= y_list[i])
		not_select_pt = np.invert(select_pt)

		## flux map 1D
		flux_all_row = sum(data[select_pt])/sum(cell_areas[select_pt])
		flux_map_1D = np.append(flux_map_1D, flux_all_row)
		## flux map 2D
		# xx = centers[select_pt,0].astype(float)
		# y = (y_list[i-1] + y_list[i])/2.
		# data_row = sortarray(y, x_list, xx, data[select_pt], cell_areas[select_pt])
		# flux_map_2D = np.append(flux_map_2D, data_row)

		centers = centers[not_select_pt,:]
		data = data[not_select_pt,:]
		cell_areas = cell_areas[not_select_pt]
		i+=1

	# flux_map_2D = flux_map_2D.reshape(len(flux_map_2D)/3, 3)
	# title = np.array(['X (m)', 'Y (m)', 'flux (W/m2)'])
	# flux_map_2D = np.vstack((title, flux_map_2D))

	flux_map_1D = np.hstack((len(flux_map_1D), flux_map_1D))

	return flux_map_2D, flux_map_1D

def trianglearea(vertices):
	'''
	Return the area of the triangle polygon following shoelace formula or algorithm
	NB:
	'''
	v1 = vertices[0]
	v2 = vertices[1]
	v3 = vertices[2]
	return 0.5*np.abs(v1[0]*(v2[1] - v3[1]) + v2[0]*(v3[1] - v1[1]) + v3[0]*(v1[1] - v2[1]))

def localcoordinates(points_3D, vertices):
	'''
	Convert 3D coordinate points in 2D local coordinates of the corresponding flat surface
	'''
	n_dimension = len(points_3D[0])
	i = 0
	while i < n_dimension:
		min_coord = min(points_3D[:,i].astype(float))
		max_coord = max(points_3D[:,i].astype(float))
		if ((min_coord == 0.0) and (max_coord == 0.0)):
			points_2D = np.delete(points_3D, i, 1)
			i += n_dimension
			n_dimension -= 1
		i += 1

	if n_dimension > 2:
		n_points = len(points_3D)
		p1 = vertices[0]
		p2 = vertices[1]
		p3 = vertices[2]
		origin, vbase_x, vbase_y, vbase_z = basevectors(p1, p2, p3)

		points_2D = np.zeros((n_points, 2))
		for i in range(n_points):
			p_x, p_y = projection(points_3D[i], origin, vbase_x, vbase_y, vbase_z)
			points_2D[i] = [p_x, p_y]

	return points_2D


if __name__=='__main__':

    p1=np.array([0,0,0])
    p2=np.array([0,1,0])
    p3=np.array([0,0,1])
    basevectors(p1,p2,p3)
