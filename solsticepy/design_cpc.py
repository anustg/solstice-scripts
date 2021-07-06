import os
import sys
import numpy as np
import time

class CPC:
        '''
        The Central Receiver System (CRS) model includes three parts:
        the sun, the field and the receiver.
        '''

        def __init__(self):
                '''
                Arguements:
                        casedir : str, the directory of the case
                '''
                self.y_min_para=0.0
                self.z_min_para=0.0
                self.z_max_para=0.0


        def regularpolygon(self, poly_half_w=2., poly_half_l=2., cpc_nfaces=6, bool=0):
            '''
            Arguments:
                (1) receiver  :   str, 'flat', 'cylinder' or directory of the 'stl', type of the receiver
                (2) poly_w   : float, width of the polygon (m)
                (3) poly_h   : float, length of the polygon (m)
                (4) cpc_nfaces     : int, number of faces of the 2D-crossed CPC
                (5) bool   : 0: gives the vertex of the polygon, 1: gives the center of the polygon edges.
            '''
            cpc_nfaces = int(cpc_nfaces)
            gamma = np.pi/cpc_nfaces
            poly_vertex=np.zeros((cpc_nfaces,2))
            for i in range(cpc_nfaces):
                temp = (gamma*(bool + i*2.))
                cosGamma = np.cos(bool*gamma)
                poly_vertex[i] = [-poly_half_w /cosGamma*np.sin(temp), poly_half_l /cosGamma*np.cos(temp)]

            return poly_vertex

        def intersectionpoint(self, p_A=1., p_B=1., p_C=1., focal_length=1.):
            '''
            Gives the intersection point between a parabola (z=y^2/(4f)) and a line (z=m*y+z0), with
            - f the focal length of the parabola
            - p_A = 1/(4f)
            - p_B = -m
            - p_C = -z0
            From equation: p_A y^2 + p_B y + p_C = 0
            '''
            delta = p_B*p_B - 4.*p_A*p_C
            y_inter_para = (-p_B + np.sqrt(delta))/(2.*p_A)
            z_inter_para  = self.parabola(y_inter_para, focal_length)

            return (y_inter_para, z_inter_para)

        def parabola(self, y_para, focal_length):
            '''
            Argument:
                (1) y_para  :   y coordinate of the parabola defined in yOz plan (m)
                (2) focal_length   : focal length of the parabola (m)
            '''
            z_para=y_para*y_para/(4.*focal_length)

            return z_para

        def hyperboloid(self, x_hyper, y_hyper, focal_image, focal_real):
            '''
            calculate the z coordinate of the hyperboloid:
            (x^2+y^2)/a^2 - (z+z0-g/2)^2/b^2 + 1 = 0
            with the vetex of the hyperboloid located in the xOy plan.
            '''
            g_para = focal_image+focal_real
            f_para = focal_real/g_para
            a_para_sqrt = (g_para**2)*(f_para-f_para**2)
            b_para = g_para*(f_para-1/2.)
            z0_hyper = abs(b_para) + g_para/2.

            z_hyper = b_para*np.sqrt((x_hyper**2+y_hyper**2)/a_para_sqrt+1) - z0_hyper+g_para/2.
            return z_hyper

        def thetagammarotations(self, y_para, z_para, theta, gamma):
            '''
            Gives transformation coordinates of rotation around X-axis of theta
            followed by rotation around Z-axis of gamma angle.

            Arguments:
              (1) y_para : y coordinate of the parabola defined in yOz plan (m)
              (2) z_para : z coordinate of the parabola defined in yOz plan (m)
              (3) theta : acceptance angle of the CPC (radian) (rotation angle around x-axis)
              (4) gamma : angle of rotation of the CPC face around the z-axiz (radian)
            '''
            x_rot_theta=0.
            y_rot_theta=y_para*np.cos(theta) - z_para*np.sin(theta)
            z_rot_theta=y_para*np.sin(theta) + z_para*np.cos(theta)

            x_rot_gamma=x_rot_theta*np.cos(gamma) - y_rot_theta*np.sin(gamma)
            y_rot_gamma=x_rot_theta*np.sin(gamma) + y_rot_theta*np.cos(gamma)
            z_rot_gamma=z_rot_theta
            rot_vector=np.array([x_rot_gamma, y_rot_gamma, z_rot_gamma])

            return rot_vector

        def xyztranslations(self, theta, gamma, x_position, y_position, z_position):
            '''
            Gives the translation of a CPC face for yaml file creation.

            Arguments:
              (1) theta : acceptance angle of the CPC (radian) (rotation angle around x-axis)
              (2) gamma : angle of rotation of the CPC face around the z-axis (radian)
              (3) x_position : x coordinate of the final position in the global coordinate system
              (4) y_position : y coordinate of the final position in the global coordinate system
              (4) z_position : z coordinate of the final position in the global coordinate system
            '''
            rot_vector = self.thetagammarotations(self.y_min_para, self.z_min_para, theta, gamma)
            x_trans=-rot_vector[0]+x_position
            y_trans=-rot_vector[1]+y_position
            z_trans=-rot_vector[2]+z_position

            return (x_trans, y_trans, z_trans)

        def cpcfaceclipping(self, cpc_nZ, theta, focal_length, cpc_nfaces, receiver_radius):
            '''
            Return the clipping polygon (vertex) of the CPC face

            Arguments:
                (1) cpc_nZ : number of iterations along the clipping polygon
                (2) theta : acceptance angle of the CPC (radian) (rotation angle around x-axis)
                (3) focal_length   : focal length of the parabola (m)
            '''
            delta_Z = (self.z_max_para - self.z_min_para)/cpc_nZ
            _, y_trans, z_trans= self.xyztranslations(theta, 0.0, 0.0, receiver_radius, 0.0)
            poly_vertex=np.zeros((2*(cpc_nZ+1),2))
            margin = 0.1
            increment = np.tan(np.pi/cpc_nfaces)
            z_CPC = 0.
            i = 0
            z_para = self.z_min_para - delta_Z
            #for i in range(cpc_nZ+1):
            #while z_CPC <= self.h_CPC:
            while abs(self.h_CPC - z_CPC) >= 0.0001:
                # Calculate the local y coordinate
                z_para += delta_Z
                y_para = np.sqrt(4*focal_length*z_para)
                # Calculate the radius of the CPC at the height z_para
                _, y_CPC, z_CPC = self.thetagammarotations(y_para, z_para, theta, 0.0)
                z_CPC += z_trans
                if z_CPC <= self.h_CPC + 0.000001:
                    # Calculate the local/global x coordinate
                    x_para = (y_CPC+y_trans)*increment+margin
                    poly_vertex[i]=[x_para, y_para]
                    poly_vertex[2*cpc_nZ+1-i]=[-x_para, y_para]
                    i += 1
                else:
                    z_para -= delta_Z
                    delta_Z/=2

            return poly_vertex

        def dependantparameters(self, rec_param):
            '''
            Return the clipping polygon (vertex) of the CPC face
            '''
            rec_w=rec_param[0]
            rec_l=rec_param[1]
            rec_z=rec_param[2]
            cpc_nfaces=int(rec_param[4])
            cpc_theta_deg=rec_param[5]
            cpc_h=rec_param[6]

            cpc_nfaces = int(cpc_nfaces)
            if cpc_nfaces==4:
                self.rec_x = rec_w/2.
                self.rec_y = rec_l/2.
                self.rec_radius = np.minimum(self.rec_x,self.rec_y)
            else:
                self.rec_x = np.sqrt(rec_w*rec_w+rec_l*rec_l)/2.
                self.rec_y = self.rec_x
                self.rec_radius = self.rec_x


            # Parameters of the parabola
            self.theta = cpc_theta_deg*np.pi/180.
            self.focal_length = self.rec_radius*(np.sin(self.theta) + 1)

            # Lowest and highest point of the parabola curve of the CPC face given in local coordinates of the parabola
            p_A = 1.
            p_B = 4*self.focal_length*np.tan(self.theta)
            p_C = -4*self.focal_length*self.focal_length
            self.y_min_para, self.z_min_para = self.intersectionpoint(p_A, p_B, p_C, self.focal_length)
            ''' OLD
            if cpc_h<=0.:
                p_B = -4.0*self.focal_length/np.tan(2.*self.theta)
                self.y_max_para, self.z_max_para = self.intersectionpoint(p_A, p_B, p_C, self.focal_length)
            else:
                self.z_max_para = self.z_min_para+cpc_h
                self.y_max_para = np.sqrt(4*self.focal_length*self.z_max_para)
            # height of the CPC
            _, _, z_trans = self.xyztranslations(self.theta, 0.0, 0.0, self.rec_radius, 0.0)
            _, _, z_CPC = self.thetagammarotations(self.y_max_para, self.z_max_para, self.theta, 0.0)
            self.h_CPC = z_trans + z_CPC
            print('calculated height: ', self.h_CPC)
            h_CPC = self.rec_radius*(1+1/np.sin(self.theta))/np.tan(self.theta)
            print('theoretical critical height: ', h_CPC)
            '''
            h_CPC = self.rec_radius*(1+1/np.sin(self.theta))/np.tan(self.theta)
            p_B = -4.0*self.focal_length/np.tan(2.*self.theta)
            _, self.z_max_para = self.intersectionpoint(p_A, p_B, p_C, self.focal_length)
            if cpc_h<=0.:
                self.h_CPC = h_CPC
            else:
                self.h_CPC = cpc_h

            print('theoretical critical height: ', h_CPC)
            print('selected height: ', self.h_CPC)

        def beamdowncomponents(self, rec_param):
            '''
            Generate YAML files for the Solstice simulation

            Arguments:
                # Receiver (flat)
                (1) rec_w     : float,  width of the receiver (m)
                (2) rec_l     : float, length of the receiver (m)
                (3) rec_z     : float, z (vertical) location of the receiver (m)
                (4) rec_grid  : number of slices for solstice flux map
                # Compound Parabolic Concentrator (CPC)
                (5) cpc_nfaces : int, number of faces of the CPC
                (6) cpc_theta_deg : float, acceptance angle of CPC (deg)
                (7) cpc_h     : float, vertical distance of the maximum and minimum point of the parabola
                in local coordinate system of the parabola in the xOy plan
                WARNING: cpc_h is different from the total height of the CPC
                (8) cpc_nZ     : int, number of number of incrementation for the clipping polygon for the construction of each CPC face
                # Secondary Reflector (hyperboloid)
                (9) aim_z   : float, z (vertical) location of the heliostats' aiming point (m)
                (10) secref_z    : flaot, z (vertical) location of the apex of the hyperboloid secondary reflector (m)
                (11) r_f	: float, secondary mirror and CPC reflectivity, e.g. 0.9
                (12) slope_error	: float, slope error of secondary mirror and CPC refelctivity (?)
                (13) secref_vert	: array, (x,y) clipping polygon of the secondary reflector
            '''

            rec_z=rec_param[2]
            rec_grid=int(rec_param[3])  # it assumes equal number of slices in x and y directions
            cpc_nfaces=int(rec_param[4])
            cpc_theta_deg=rec_param[5]
            cpc_nZ=int(rec_param[7])
            aim_z=rec_param[8]
            secref_z=rec_param[9]
            r_f=rec_param[10] # front
            slope_error=rec_param[11]
            secref_vert=rec_param[12]
            self.dependantparameters(rec_param[:7])

            iyaml='\n'

            #
            #    Materials
            #
            # CREATE a specular material for the secondary reflector
            iyaml+='- material: &%s\n' % 'material_sec_mirror'
            iyaml+='   front:\n'
            iyaml+='     mirror: {reflectivity: %6.4f, slope_error: %15.8e }\n' % (r_f, slope_error)
            iyaml+='   back:\n'
            iyaml+='     matte: {reflectivity: 0.0 }\n'
            iyaml+='\n'

            #
            #    Entities
            #
            # Receiver Polygon
            rec_vert=self.regularpolygon(self.rec_x, self.rec_y, cpc_nfaces, 1)
            iyaml+='- entity:\n'
            iyaml+='    name: %s\n' % 'receiver'
            iyaml+='    primary: 0\n'
            iyaml+='    transform: { translation: %s, rotation: %s }\n' % ([0, 0, rec_z], [0, 0, 0])
            iyaml+='    geometry:\n'
            iyaml+='    - material: *%s\n' % 'material_target'
            iyaml+='      plane:\n'
            iyaml+='        slices:  %d\n' % rec_grid
            iyaml+='        clip: \n'
            iyaml+='        - operation: AND \n'
            iyaml+='          vertices: \n'
            for i in range(cpc_nfaces):
                iyaml+='            - %s\n' % ([rec_vert[i][0],rec_vert[i][1]])
            iyaml+='\n'
            #
            # CREATE a virtual target entity at the receiver level
            virt_vert = np.array([ [-1.0, 1.0], [-1.0, -1.0], [1.0, -1.0], [1.0, 1.0] ])
            virt_vert = virt_vert*(self.rec_radius * 50.0)
            slices = 4
            iyaml+='\n- entity:\n'
            iyaml+='    name: virtual_target\n'
            iyaml+='    primary: 0\n'
            iyaml+='    transform: { translation: %s, rotation: %s }\n' % ([0, 0, -5], [0., 0, 0])
            iyaml+='    geometry: \n'
            iyaml+='      - material: *%s\n' % 'material_virtual'
            iyaml+='        plane: \n'
            iyaml+='          slices: %d\n' % slices
            iyaml+='          clip: \n'
            iyaml+='          - operation: AND \n'
            iyaml+='            vertices: \n'
            for i in range(len(virt_vert)):
                iyaml+='            - %s\n' % ([virt_vert[i][0],virt_vert[i][1]])
            iyaml+='\n'
            #
            # Secondary Reflector
            # focal image: distance between hyperbola vertex and aiming point of heliostats
            # focal real: distance between hyperbola vertex and aiming point of secondary reflector (aperture of the CPC)
            focal_image = aim_z - secref_z
            focal_real = abs(rec_z + self.h_CPC - secref_z) # aim at the aperture of the CPC
            print('HYPERBOLA focal IMAGE: ', focal_image)
            print('HYPERBOLA focal REAL: ', focal_real)
            iyaml+='- entity:\n'
            iyaml+='    name: %s\n' % 'secondary_reflector'
            iyaml+='    primary: 0\n'
            iyaml+='    transform: { translation: %s, rotation: %s }\n' % ([0, 0, secref_z], [0., 0, 0])
            iyaml+='    geometry:\n'
            iyaml+='    - material: *%s\n' % 'material_sec_mirror'
            iyaml+='      hyperbol:\n'
            iyaml+='        focals: &hyperbol_focals { real: %s, image: %s }\n' % (focal_real, focal_image)
            iyaml+='        clip: \n'
            iyaml+='        - operation: AND \n'
            iyaml+='          vertices: \n'
            for i in range(len(secref_vert)):
                iyaml+='            - %s\n' % ([secref_vert[i][0],secref_vert[i][1]])
            iyaml+='\n'
            #
            # CREATE a virtual target entity above the secondary reflector
            xmax=max(secref_vert[:,0])
            ymax=max(secref_vert[:,1])
            virt_vert = np.array([ [-xmax, ymax], [-xmax, -ymax], [xmax, -ymax], [xmax, ymax] ])*50
            zmax=self.hyperboloid(xmax, ymax, focal_image, focal_real)
            print('X HYPERBOLA: ', xmax)
            print('Y HYPERBOLA: ', ymax)
            print('Z HYPERBOLA: ', zmax)
            z_trans = secref_z+zmax+5
            if z_trans < aim_z:
                z_trans = aim_z
            iyaml+='\n- entity:\n'
            iyaml+='    name: virtual_sec_ref\n'
            iyaml+='    primary: 0\n'
            iyaml+='    transform: { translation: %s, rotation: %s }\n' % ([0, 0, z_trans], [0., 180., 0])
            iyaml+='    geometry: \n'
            iyaml+='      - material: *%s\n' % 'material_virtual'
            iyaml+='        plane: \n'
            iyaml+='          slices: %d\n' % slices
            iyaml+='          clip: \n'
            iyaml+='          - operation: AND \n'
            iyaml+='            vertices: \n'
            for i in range(len(virt_vert)):
                iyaml+='            - %s\n' % ([virt_vert[i][0],virt_vert[i][1]])
            iyaml+='\n'
            #
            # CREATE the receiver objects in receiver yaml
            # receiver at the CPC exit and 2 virtuals surfaces
            rcv = '\n'
            rcv+='- name: receiver \n'
            rcv+='  side: %s \n' % 'FRONT_AND_BACK'
            rcv+='  per_primitive: %s \n' % 'INCOMING_AND_ABSORBED'
            rcv+='- name: virtual_sec_ref \n'
            rcv+='  side: %s \n' % 'FRONT'
            rcv+='  per_primitive: %s \n' % 'INCOMING'
            rcv+='- name: virtual_target \n'
            rcv+='  side: %s \n' % 'FRONT'
            rcv+='  per_primitive: %s \n' % 'INCOMING'
            #
            # CREATE the facets of the CPC
            # For each face, the x,y,z translations and the clipping polygon are calculated
            # The clipping polygon is individual to a face if nCPC = 4 and the receiver is not a regular polygon (xRec!=yRec)
            cpc_face_pos=self.regularpolygon(self.rec_x, self.rec_y, cpc_nfaces, 0)
            cpc_face_polygon=self.cpcfaceclipping(cpc_nZ, self.theta, self.focal_length, cpc_nfaces, self.rec_radius)
            rec_vert = np.insert(rec_vert,0,rec_vert[cpc_nfaces-1],axis=0)
            increment = np.zeros(cpc_nfaces)
            if self.rec_x != self.rec_y:
                for i in range(cpc_nfaces):
                    increment[i] = np.linalg.norm(rec_vert[i]-rec_vert[i+1])/2 - self.rec_radius
            # implement individual facets of the CPC
            for i in range(cpc_nfaces):
                gamma = i*2.*np.pi/cpc_nfaces
                iyaml+='- entity:\n'
                iyaml+='    name: %s\n' % ('CPC_face_'+str(i))
                x_trans, y_trans, z_trans = self.xyztranslations(self.theta, gamma, cpc_face_pos[i][0], cpc_face_pos[i][1], rec_z)
                iyaml+='    transform: { translation: %s }\n' % ([x_trans, y_trans, z_trans])
                iyaml+='    children: \n'
                iyaml+='    - name: RotationGamma\n'
                iyaml+='      transform: { rotation: %s } \n' % ([0, 0, gamma*180.0/np.pi])
                iyaml+='      children: \n'
                iyaml+='      - name: RotationTheta\n'
                iyaml+='        primary: 0\n'
                iyaml+='        transform: {rotation: %s } \n' % ([cpc_theta_deg, 0, 0.])
                iyaml+='        geometry:\n'
                iyaml+='        - material: *%s\n' % 'material_sec_mirror'
                iyaml+='          parabolic-cylinder:\n'
                iyaml+='            focal:  %s\n' % self.focal_length
                iyaml+='            clip: \n'
                iyaml+='            - operation: AND \n'
                iyaml+='              vertices: \n'
                midlength = int(len(cpc_face_polygon)/2)
                for j in range(midlength):
                    iyaml+='               - %s\n' % ([cpc_face_polygon[j][0]+increment[i],cpc_face_polygon[j][1]])
                for j in range(midlength,midlength*2):
                    iyaml+='               - %s\n' % ([cpc_face_polygon[j][0]-increment[i],cpc_face_polygon[j][1]])
                iyaml+='\n'


            return iyaml, rcv

if __name__=='__main__':
        start=time.time()
#        casedir='./test-CPC-design'

# enter the parameters
        rec_w=11.
        rec_l=4.
        rec_z=0.
        rec_grid=20.
        n_CPC_faces=4
        theta_deg=20.
        cpc_h=0.
        n_Z=20
        aim_z = 30.
        sec_z = 22.
        r_f= 1. # front
        slope_error = 0.
        secref_vert = np.array([[-15,25],[-15,-10],[15,-10],[15,25]])

        cpc=CPC()
        cpc.beamdowncomponents(rec_w, rec_l, rec_z, rec_grid, n_CPC_faces, theta_deg, cpc_h, n_Z, aim_z, sec_z, r_f, slope_error, secref_vert)

        eta=1.
        print('total efficiency:', eta)

        end=time.time()
        print('total time %.2f'%((end-start)/60.), 'min')
