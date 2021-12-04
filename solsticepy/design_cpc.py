import os
import sys
import numpy as np
import time
from scipy.integrate import quad

class CPC:
        '''
        The Central Receiver System (CRS) model includes three parts:
        the sun, the field and the receiver.
        '''

        def __init__(self):
            '''
            Arguments:
                    casedir : str, the directory of the case
            '''


        def regularpolygon(self, poly_half_w=2., poly_half_l=2., bool=0):
            '''
            Gives the vertices of a polygon with n faces or the center of the polygon edges.
            Arguments:
                (1) poly_w   : float, width of the polygon (m)
                (2) poly_h   : float, length of the polygon (m)
                (3) bool   : 0: gives the vertices of the polygon, 1: gives the center of the polygon edges.
            '''
            gamma = np.pi/self.cpc_nfaces
            poly_vertex=np.zeros((self.cpc_nfaces,2))
            for i in range(self.cpc_nfaces):
                temp = (gamma*(bool + i*2.))
                cosGamma = np.cos(bool*gamma)
                poly_vertex[i] = [-poly_half_w /cosGamma*np.sin(temp), poly_half_l /cosGamma*np.cos(temp)]

            return poly_vertex

        def intersectionpoint(self, p_A=1., p_B=1., p_C=1.):
            '''
            Gives the intersection point between a parabola (z=y^2/(4f)) and a line (z=m*y+z0), with
            - f the focal length of the parabola
            - p_A = 1/(4f)
            - p_B = -m
            - p_C = -z0
            From equation: p_A y^2 + p_B y + p_C = 0
            '''
            delta = (p_B**2) - 4.*p_A*p_C
            y_inter_para = (-p_B + np.sqrt(delta))/(2.*p_A)
            z_inter_para  = (y_inter_para**2) / (4.*self.focal_length)

            return (y_inter_para, z_inter_para)


        def hyperboloid(self, x_hyper, y_hyper):
            '''
            Calculate the z coordinate of the hyperboloid:
            (x^2+y^2)/a^2 - (z+z0-g/2)^2/b^2 + 1 = 0
            with the apex of the hyperboloid located in the xOy plan.
            '''
            g_para = self.focal_image + self.focal_real
            f_para = self.focal_real / g_para
            a_para_sqrt = (g_para**2) * (f_para-f_para**2)
            b_para = g_para * (f_para-1/2.)
            z0_hyper = abs(b_para) + g_para/2.

            z_hyper = b_para*np.sqrt((x_hyper**2+y_hyper**2)/a_para_sqrt + 1) - z0_hyper+g_para/2.
            return z_hyper

        def xzrotations(self, y_para, z_para, theta, gamma):
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

        def cpcfacetranslation(self, y_min_para, z_min_para, gamma, x_position, y_position, z_position):
            '''
            Gives the translation of a CPC face for yaml file creation.
            Arguments:
              (1) gamma : angle of rotation of the CPC face around the z-axis (radian)
              (2) x_position : x coordinate of the final position in the global coordinate system
              (3) y_position : y coordinate of the final position in the global coordinate system
              (4) z_position : z coordinate of the final position in the global coordinate system
            '''
            rot_vector = self.xzrotations(y_min_para, z_min_para, self.theta, gamma)
            x_trans=-rot_vector[0]+x_position
            y_trans=-rot_vector[1]+y_position
            z_trans=-rot_vector[2]+z_position

            return (x_trans, y_trans, z_trans)

        def cpcfaceclipping(self, y_min_para, z_min_para, z_max_para):
            '''
            Gives the clipping polygon (vertices) of a CPC face
            Arguments:
                (1) y_min_para : y-coordinate of the lowest point of the parabola curve of the CPC face given in local coordinates of the parabola
                (2) z_min_para : z-coordinate of the lowest point of the parabola curve of the CPC face given in local coordinates of the parabola
                (3) z_max_para : z-coordinate of the highest point of the parabola curve of the CPC face given in local coordinates of the parabola
            '''
            delta_Z = (z_max_para - z_min_para)/self.cpc_nZ
            _, y_trans, z_trans= self.cpcfacetranslation(y_min_para, z_min_para, 0.0, 0.0, self.rec_radius, 0.0)

            margin = 0.1
            increment = np.tan(np.pi/self.cpc_nfaces)

            z_para = z_min_para
            y_para = np.sqrt(4*self.focal_length*z_para)
            _, y_CPC, z_CPC = self.xzrotations(y_para, z_para, self.theta, 0.0)
            x_para = (y_CPC+y_trans)*increment+margin
            pos_vertex = np.array([[x_para,y_para]])
            neg_vertex = np.array([[-x_para,y_para]])
            while abs(self.h_CPC - z_CPC) >= 0.0001:
                # Calculate the local y coordinate
                z_para += delta_Z
                y_para = np.sqrt(4*self.focal_length*z_para)
                # Calculate the radius of the CPC at the height z_para
                _, y_CPC, z_CPC = self.xzrotations(y_para, z_para, self.theta, 0.0)
                z_CPC += z_trans
                if z_CPC <= self.h_CPC + 0.000001:
                    # Calculate the local/global x coordinate
                    x_para = (y_CPC+y_trans)*increment+margin
                    pos_vertex = np.vstack((pos_vertex,[[x_para, y_para]]))
                    neg_vertex = np.vstack((neg_vertex,[[-x_para, y_para]]))
                else:
                    z_para -= delta_Z
                    delta_Z/=2

            pos_vertex = np.vstack((pos_vertex,neg_vertex[::-1]))
            return pos_vertex

        def intersectionlinehyperbol(self, a_hyper, b_hyper, m_angle, z0_line):
            '''
            intersection point between hyperbola and line:
            z^2/a^2 - y^2/b^2 = 1
            z = m*y + z0
            '''
            if m_angle == 0:
                y_inter = 0.

            else:
                m_line = -1 / np.tan(m_angle * np.pi/180.)
                a_over_m_sqrt = a_hyper**2 / m_line**2

                if abs(m_line) < 0.000000001: # abs(m_angle) == 90:
                    y_inter = np.sqrt(a_hyper**2 * ((z0_line**2)/(b_hyper**2)-1))

                else:
                    aCoef_poly = b_hyper**2 - a_over_m_sqrt
                    bCoef_poly = 2*z0_line * a_over_m_sqrt
                    cCoef_poly = - a_over_m_sqrt * (z0_line**2) - (a_hyper**2 * b_hyper**2)
                    delta = bCoef_poly**2 - 4*aCoef_poly*cCoef_poly
                    root1 = (-bCoef_poly+np.sqrt(delta)) / (2*aCoef_poly)
                    root2 = (-bCoef_poly-np.sqrt(delta)) / (2*aCoef_poly)
                    # If both roots are positives, the line intersects the upper sheet of the hyperbola twice,
                    # Else, the line intersects the upper sheet of the hyperbola once and potentially the lower sheet once
                    if root1>0. and root2>0. and abs(m_angle) < 90.:
                        z_inter = min(root1,root2)
                    else:
                        z_inter = max(root1,root2)
                    y_inter = (z_inter-z0_line)/m_line
                    #check = (z_inter**2)/a_hyper**2 - (y_inter**2)/b_hyper**2
                    #print('equation resolution should be 1, and it is: ', check)

            return y_inter


        def dependantparameters(self, rec_param):
            '''
            Calculate the characteristics parameters of CPC, receiver, secondary mirror
            '''
            rec_w=rec_param[0]
            rec_l=rec_param[1]
            self.rec_z=rec_param[2]
            self.rec_grid=int(rec_param[3])  # it assumes equal number of slices in x and y directions
            self.cpc_nfaces=int(rec_param[4])
            cpc_theta_deg=rec_param[5]
            cpc_h_ratio=rec_param[6]
            self.cpc_nZ=int(rec_param[7])
            self.aim_z=rec_param[8]
            secref_inv_eccen=rec_param[9]
            self.tilt_secref=rec_param[10]
            self.rho_secref=rec_param[11]
            self.rho_cpc=rec_param[12]
            self.slope_error=rec_param[13]

            assert self.cpc_nfaces > 2, 'The number of faces for the CPC should be minimum 3, and it is {cpc_nfaces=}'
            assert 0<=secref_inv_eccen<=1, 'The inverse eccentricity of the hyperbole must be between 0 and 1, and it is {secref_inv_eccen=}'

            self.receiverparameters(rec_w,rec_l)
            self.cpcparameters(cpc_theta_deg, cpc_h_ratio)
            self.secrefparameters(secref_inv_eccen)


        def receiverparameters(self, rec_w=2., rec_l=2.):
            '''
            Calculate Receiver critical dimensions
            rec_radius_ref:	receiver radius in the direction of theta (CPC half acceptance angle), it is used for the design criteria
            '''
            if self.cpc_nfaces is 4:
                self.rec_x = rec_w/2.
                self.rec_y = rec_l/2.
                self.rec_radius = np.minimum(self.rec_x,self.rec_y)
            else:
                self.rec_x = np.sqrt(rec_w**2 + rec_l**2) / 2.
                self.rec_y = self.rec_x
                self.rec_radius = self.rec_x


        def cpcparameters(self, cpc_theta_deg, cpc_h_ratio=1.):
            '''
            Calculate Parameters of the parabola (CPC height, acceotance angle in rad, focal length)
            '''
            self.theta = cpc_theta_deg * np.pi/180.
            self.h_CPC = self.rec_radius * (1 + 1/np.sin(self.theta)) / np.tan(self.theta)
            self.h_CPC *= cpc_h_ratio

            print('Heigth of the CPC: ', self.h_CPC)


        def secrefparameters(self, secref_inv_eccen):
            '''
            Calculate characteristics parameters of the hyperboloid with tilted axis: 2 foci and apex position with receiver positioned at the center
            focal image: distance between hyperbola apex and aiming point of heliostats
            focal real: distance between hyperbola apex and aiming point of secondary reflector (aperture of the CPC)
            '''
            secref_angle = self.tilt_secref * np.pi/180.
            real_foci_z = (self.h_CPC+self.rec_z)

            #assert self.aim_z > real_foci_z, 'The imaginary foci of the hyperbol is lower than its real foci'

            aim_y =  (self.aim_z-real_foci_z) * np.tan(secref_angle)
            foci_dist = np.sqrt((self.aim_z-real_foci_z)**2 + aim_y**2)/2.
            self.focal_image = (1-secref_inv_eccen) * foci_dist
            self.focal_real = 2*foci_dist - self.focal_image

            self.secref_z = real_foci_z + self.focal_real * np.cos(secref_angle)
            self.secref_y = self.focal_real * np.sin(secref_angle)

            print('HYPERBOLA focal IMAGE: ', self.focal_image)
            print('HYPERBOLA focal REAL: ', self.focal_real)


        def secrefpolygonsclipping(self, secref_vert):
            '''
            Calculate the clipping polygon of the secondary reflector (dimensions)
            Calculate the clipping polygon of the virtual surface associated to secondary reflector (dimensions)
            if length of secref_vert is 2 (aperture_angle_x, aperture_angle_y), the clipping polygon is calculated with the 2 rim angles
            if aperture_angle_y = None, the clipping polygon is a circle
            '''
            center=None
            ## Secondary reflector
            if len(secref_vert) is 3:
                aperture_angle_x=abs(secref_vert[0])
                aperture_angle_y=secref_vert[1]
                secref_offset=secref_vert[2]

                foci_dist = (self.focal_image + self.focal_real) / 2.
                a_hyper = (self.focal_real-foci_dist)
                b_hyper = np.sqrt(self.focal_image*self.focal_real)
                m_angle = aperture_angle_x/2.
                x_inter_max = self.intersectionlinehyperbol(a_hyper, b_hyper, m_angle, foci_dist)

                if aperture_angle_y is None:
                    rim_angle_y = aperture_angle_x/2.
                else:
                    rim_angle_y = abs(aperture_angle_y)/2.

                m_angle = secref_offset - rim_angle_y
                y_inter_min = self.intersectionlinehyperbol(a_hyper, b_hyper, m_angle, foci_dist)
                m_angle = secref_offset + rim_angle_y
                y_inter_max = self.intersectionlinehyperbol(a_hyper, b_hyper, m_angle, foci_dist)

                if aperture_angle_y is None:
                    secref_vert = np.array([x_inter_max])
                    if secref_offset == 0.:
                        y_center = 0.
                    else:
                        y_center = (y_inter_max + y_inter_min)/2.
                    center=np.array([0., y_center])
                    y_max=max(x_inter_max,abs(y_inter_max),abs(y_inter_min))
                else:
                    secref_vert = np.array([[0.,y_inter_max],[-x_inter_max,y_inter_max],[-x_inter_max,0.],[-x_inter_max,y_inter_min],[0.,y_inter_min],[x_inter_max,y_inter_min],[x_inter_max,0.],[x_inter_max,y_inter_max]])
                    y_max=max(abs(y_inter_max),abs(y_inter_min))

            else:
                x_inter_max=max(abs(secref_vert[:,0]))
                y_max=max(abs(secref_vert[:,1]))

            ## Vitual surface above secondary reflector
            virt_vert = np.array([ [-x_inter_max, y_max], [-x_inter_max, -y_max], [x_inter_max, -y_max], [x_inter_max, y_max] ])*50
            z_max = self.focal_real + self.hyperboloid(x_inter_max, y_max)
            asb_secref_tilt = abs(self.tilt_secref * np.pi/180.)
            z_max = z_max * np.cos(asb_secref_tilt) + y_max * np.sin(asb_secref_tilt)
            z_max += (self.h_CPC+self.rec_z)

            return secref_vert, virt_vert, z_max, center


        def cpcfaceparameters(self):
            '''
            Return the parameters used to create each CPC face with YAML files for the Solstice simulation:
            cpc_face_polygon: polygon clipping of an individual CPC face
            cpc_faces_trans: (x,y,z) for each face of the CPC
            '''
            # Lowest and highest point of the parabola curve of the CPC face given in local coordinates of the parabola
            self.focal_length = self.rec_radius * (np.sin(self.theta) + 1)
            p_B = 4. * self.focal_length * np.tan(self.theta)
            p_C = -4. * self.focal_length * self.focal_length
            y_min_para, z_min_para = self.intersectionpoint(1., p_B, p_C)
            p_B = -4. * self.focal_length / np.tan(2.*self.theta)
            _, z_max_para = self.intersectionpoint(1., p_B, p_C)

            # For each face, calculate the x,y,z translations and the clipping polygon
            cpc_face_polygon = self.cpcfaceclipping(y_min_para, z_min_para, z_max_para)

            cpc_faces_pos = self.regularpolygon(self.rec_x, self.rec_y, 0)
            cpc_faces_trans=np.zeros((self.cpc_nfaces,3))
            for i in range(self.cpc_nfaces):
                gamma = i*2.*np.pi/self.cpc_nfaces
                x_trans, y_trans, z_trans = self.cpcfacetranslation(y_min_para, z_min_para, gamma, cpc_faces_pos[i][0], cpc_faces_pos[i][1], self.rec_z)
                cpc_faces_trans[i] = [x_trans, y_trans, z_trans]

            return cpc_face_polygon, cpc_faces_trans


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
                (7) cpc_h_ratio     : float, ratio of critical CPC height calculated with cpc_theta_deg, [0,1]
                (8) cpc_nZ     : int, number of number of incrementation for the clipping polygon for the construction of each CPC face
                # Secondary Reflector (hyperboloid)
                (9) aim_z   : float, z (vertical) location of the heliostats' aiming point (m)
                (10) secref_inv_eccen    : float, hyperboloid inverse eccentricity: ratio of the apex distance over the foci distance to the origin, must be between 0 and 1
                (11) tilt_secref    : float, angle (in degree) of the tilted axis of the hyperboloid, from the vertical to the North (+) or South (-), [-180,180]
                (12) rho_secref	: float, secondary mirror reflectivity property, [0,1]
                (13) rho_cpc	: float, CPC reflectivity property, [0,1]
                (14) slope_error	: float, slope error of secondary mirror and CPC refelctivity (?)
                (15) secref_vert	: array, (x_i,y_i) clipping polygon of the secondary reflector
                # if length of array is = 2, rim angles along X and Y axis used to calculate the clipping polygon of the secondary reflector
            '''

            self.dependantparameters(rec_param)

            iyaml='\n'
            #
            #    Materials
            #
            # CREATE a specular material for the secondary reflector
            iyaml+='- material: &%s\n' % 'material_sec_mirror'
            iyaml+='   front:\n'
            iyaml+='     mirror: {reflectivity: %6.4f, slope_error: %15.8e }\n' % (self.rho_secref, self.slope_error)
            iyaml+='   back:\n'
            iyaml+='     matte: {reflectivity: 0.0 }\n'
            iyaml+='\n'
            # CREATE a specular material for the CPC
            iyaml+='- material: &%s\n' % 'material_cpc'
            iyaml+='   front:\n'
            iyaml+='     mirror: {reflectivity: %6.4f, slope_error: %15.8e }\n' % (self.rho_cpc, self.slope_error)
            iyaml+='   back:\n'
            iyaml+='     matte: {reflectivity: 0.0 }\n'
            iyaml+='\n'

            #
            #    Entities
            #
            # Receiver Polygon
            rec_vert=self.regularpolygon(self.rec_x, self.rec_y, 1)
            iyaml+='- entity:\n'
            iyaml+='    name: %s\n' % 'receiver'
            iyaml+='    primary: 0\n'
            iyaml+='    transform: { translation: %s, rotation: %s }\n' % ([0, 0, self.rec_z], [0, 0, 0])
            iyaml+='    geometry:\n'
            iyaml+='    - material: *%s\n' % 'material_target'
            iyaml+='      plane:\n'
            iyaml+='        slices:  %d\n' % self.rec_grid
            iyaml+='        clip: \n'
            iyaml+='        - operation: AND \n'
            iyaml+='          vertices: \n'
            for i in range(self.cpc_nfaces):
                iyaml+='            - %s\n' % ([rec_vert[i][0],rec_vert[i][1]])
            iyaml+='\n'
            #
            # Virtual target entity at the receiver level
            virt_vert = np.array([ [-1.0, 1.0], [-1.0, -1.0], [1.0, -1.0], [1.0, 1.0] ]) * (self.rec_radius * 50.0 * self.h_CPC)
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
                iyaml+='            - %s\n' % ([virt_vert[i][0], virt_vert[i][1]])
            iyaml+='\n'
            #
            # Secondary Reflector
            secref_vert, virt_vert, virt_vert_z, center=self.secrefpolygonsclipping(rec_param[14])
            iyaml+='- entity:\n'
            iyaml+='    name: %s\n' % 'secondary_reflector'
            iyaml+='    primary: 0\n'
            iyaml+='    transform: { translation: %s, rotation: %s }\n' % ([0, self.secref_y, self.secref_z], [-self.tilt_secref, 0., 0.])
            iyaml+='    geometry:\n'
            iyaml+='    - material: *%s\n' % 'material_sec_mirror'
            iyaml+='      hyperbol:\n'
            iyaml+='        focals: &hyperbol_focals { real: %s, image: %s }\n' % (self.focal_real, self.focal_image)
            iyaml+='        slices: 100\n'
            iyaml+='        clip: \n'
            if len(secref_vert)>1:
                iyaml+='        - operation: AND \n'
                iyaml+='          vertices: \n'
                for i in range(len(secref_vert)):
                    iyaml+='            - %s\n' % ([secref_vert[i][0],secref_vert[i][1]])
                iyaml+='\n'
            else:
                iyaml+='        - operation: AND \n'
                iyaml+='          circle: { radius: %s, center: %s} \n' % (secref_vert[0], [center[0], center[1]])
            #
            # Virtual target entity above the secondary reflector
            if virt_vert_z < self.aim_z:
                virt_vert_z = self.aim_z
            iyaml+='\n- entity:\n'
            iyaml+='    name: virtual_sec_ref\n'
            iyaml+='    primary: 0\n'
            iyaml+='    transform: { translation: %s, rotation: %s }\n' % ([0, self.secref_y, virt_vert_z], [0., 180., 0])
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
            # CPC Facets: implement individual faces of the CPC
            cpc_face_polygon, cpc_faces_trans = self.cpcfaceparameters()
            rec_vert = np.insert(rec_vert,0,rec_vert[self.cpc_nfaces-1],axis=0)
            for i in range(self.cpc_nfaces):
                gamma = i*2.*np.pi/self.cpc_nfaces
                iyaml+='- entity:\n'
                iyaml+='    name: %s\n' % ('CPC_face_'+str(i))
                iyaml+='    transform: { translation: %s }\n' % ([cpc_faces_trans[i][0],cpc_faces_trans[i][1],cpc_faces_trans[i][2]])
                iyaml+='    children: \n'
                iyaml+='    - name: RotationGamma\n'
                iyaml+='      transform: { rotation: %s } \n' % ([0, 0, gamma*180.0/np.pi])
                iyaml+='      children: \n'
                iyaml+='      - name: RotationTheta\n'
                iyaml+='        primary: 0\n'
                iyaml+='        transform: {rotation: %s } \n' % ([self.theta*180.0/np.pi, 0, 0.])
                iyaml+='        geometry:\n'
                iyaml+='        - material: *%s\n' % 'material_cpc'
                iyaml+='          parabolic-cylinder:\n'
                iyaml+='            slices: 30\n'
                iyaml+='            focal:  %s\n' % self.focal_length
                iyaml+='            clip: \n'
                iyaml+='            - operation: AND \n'
                iyaml+='              vertices: \n'
                increment = np.linalg.norm(rec_vert[i]-rec_vert[i+1])/2 - self.rec_radius
                midlength = int(len(cpc_face_polygon) / 2)
                for j in range(midlength):
                    iyaml+='               - %s\n' % ([cpc_face_polygon[j][0]+increment,cpc_face_polygon[j][1]])
                for j in range(midlength,midlength*2):
                    iyaml+='               - %s\n' % ([cpc_face_polygon[j][0]-increment,cpc_face_polygon[j][1]])
                iyaml+='\n'
            #
            #
            # Receivers objects in receiver yaml (receiver polygon, CPC, secondary mirror and 2 virtuals surfaces)
            rcv = '\n'
            rcv+='- name: receiver \n'
            rcv+='  side: %s \n' % 'FRONT_AND_BACK'
            rcv+='  per_primitive: %s \n' % 'INCOMING_AND_ABSORBED'
            rcv+='- name: secondary_reflector \n'
            rcv+='  side: %s \n' % 'FRONT'
            rcv+='  per_primitive: %s \n' % 'INCOMING_AND_ABSORBED'
            for i in range(self.cpc_nfaces):
                rcv+='- name: %s\n' % ('CPC_face_'+str(i)+'.RotationGamma.RotationTheta')
                rcv+='  side: %s \n' % 'FRONT'
                rcv+='  per_primitive: %s \n' % 'INCOMING_AND_ABSORBED'
            rcv+='- name: virtual_sec_ref \n'
            rcv+='  side: %s \n' % 'FRONT'
            rcv+='  per_primitive: %s \n' % 'INCOMING'
            rcv+='- name: virtual_target \n'
            rcv+='  side: %s \n' % 'FRONT'
            rcv+='  per_primitive: %s \n' % 'INCOMING'

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
        cpc_h_ratio=1.
        n_Z=20
        aim_z = 30.
        secref_inv_eccen = 0.72
        tilt_secref = 20.
        rho_bd = 1. # front
        slope_error = 0.
        secref_vert = np.array([[-15,25],[-15,-10],[15,-10],[15,25]])

        bd_param=np.array([rec_w, rec_l, rec_z, rec_grid, n_CPC_faces, theta_deg, cpc_h_ratio, n_Z, aim_z,
        secref_inv_eccen, tilt_secref, rho_bd, rho_bd, slope_error, secref_vert])

        cpc=CPC()
        cpc.beamdowncomponents(bd_param)

        eta=1.
        print('total efficiency:', eta)

        end=time.time()
        print('total time %.2f'%((end-start)/60.), 'min')
