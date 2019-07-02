""" A test function to try and build an additional scoring term using the dipole maps """

import biobox as bb
import numpy as np
import sys
import JabberDock as jd
from copy import deepcopy
import scipy.spatial as spaz

class Dipole(object):
    '''
    A Dipole object contains an ensemble of points in 3D space with an associated dipole value.
    '''

    def __init__(self, d=np.array([[]]), c = np.array([[]])):
        '''
        Here we initiate our Dipole object, which must include information on the coordinates of each voxel, or at least the means to 
        build it, as well as the dipoles from those points

        :param d: The dipole held within each voxel
        :param c: The coordinate information for the initial point of each dipole (each voxel)
        ''' 
        super(Dipole, self).__init__()

        if d.ndim == 2:
             self.dipole = d
             self.coordinates = c
        else:
            raise Exception("ERROR: Your supplied dipole information that does is not n x 3 (n = number of voxels).")

    def import_dipole(self, dip_name):
        '''
        Imports a dipole map written in the tcl format. This can be written by biobox.molecule.get_dipole_map.

        :param dip_name: Name of dipole tcl file
        '''

        map_raw = np.loadtxt(dip_name, dtype=str)
        # Get the origin of the dipole, and the end points of the dipole to recalculate the dipole in each voxel
        crd_base = np.asarray(map_raw[:, 3:6]).astype(float)
        crd_end = np.asarray(map_raw[:, 8:11]).astype(float)
 
        dipole_val = crd_end - crd_base

        self.dipole = deepcopy(dipole_val)
        self.coordinates = deepcopy(crd_base)     

    def translate_dipole(self, trans = np.array((0., 0., 0.,))):
        '''
        translate a Dipole map by a given amount in place. 

        :param x: The desired x translation amounts
        :param y: The desired y translation amounts
        :param z: The desired z translation amounts
        '''
        
        self.coordinates += trans
        #self.dipole += trans

    def rotate_dipole(self, COM = np.array((0., 0., 0.)), angle = 0., axisx = 0., axisy = 0., axisz =  0., R = np.identity(3)):
        '''
        Method to rotate electron density maps about the an axis of rotation from the origin in place

        :param COM: Centre of Mass of corresponding PDB (if using PDB as basis for roto-translations)
        :param angle: Angle to shift by (must be in deg)
        :param axisx: Point in x for the axis of rotation
        :param axisy: Point in y for the axis of rotation
        :param axisz: Point in z for the axis of rotation
        :returns: The rotation matrix used to rotate the map
        '''

        if np.array_equal(R, np.identity(3)):
            axis = [axisx, axisy, axisz]

            # Convert angle to rads
            angle = angle * np.pi / 180.

            # Normalise the column vector axis of rotation
            norm = axis[0]**2 + axis[1]**2 + axis[2]**2
            axis /= np.sqrt(norm)

            # First we need to define our rotation matrix
            # We begin by defining the components of said matrix
            c = np.cos(angle)
            s = np.sin(angle)
            t = 1 - c
            x = axis[0]
            y = axis[1]
            z = axis[2]

            # We then build our rotation matrix 
            R = [[t*x*x + c, t*x*y - z*s, t*x*z + y*s], [t*x*y + z*s, t*y*y + c , t*y*z - x*s], [t*x*z - y*s, t*y*z + x*s, t*z*z + c]]

        # Move coordinates with input PDB Centre of Mass
        self.translate_dipole(-COM)

        # Multiply our coordinate system matrix by our rotation matrix - rotating the dipole structure
        vals_c = np.matmul(self.coordinates, R)
        vals_d = np.matmul(self.dipole, R)

        self.coordinates = vals_c
        self.dipole = vals_d

        # Move back to original position
        self.translate_dipole(COM)

        # Return the rotation m1= self.data.M.atomselect("*", "*", self.params.atoms)
        return np.asarray(R, dtype=float)

    def dist_cutoff_list(self, verts, cutoff_low=2, cutoff_high=5):
        '''
        Return a list of distances that are between two cutoffs

        :param verts1: coordinate for first protein (the stationary one to be docked into)
        :param verts2: coordinate for second protein (the roto-translated structure)
        :param cutoff: A cutoff to exclude interactions from - typically between 1 Ang. and 3 Ang.
        :returns: Numpy array containing four separate arrays. The first contains the distances between the two specificed cutoffs for the first set of input dipoles. The second is a boolean array where the length corresponds to the number of dipoles in the first dipole map. Where, True, the corresponding dipole has a matching dipole on the second map (the ligand). The third and fourth arrays are the same but from the second dipole map's perspective.  
        '''

        dist1 = spaz.distance.cdist(self.coordinates, verts, 'euclidean')
        dist2 = spaz.distance.cdist(verts, self.coordinates, 'euclidean') # get verts2 nearest neighbours

        #########
        # Try just getting the minimum distance
        ########
        #print dist1[dist1 < 1.0]
        cutoff_index1 = np.logical_and(cutoff_low < dist1, dist1 < cutoff_high) # Get the distances less than the cutoff   
        cutoff_index2 = np.logical_and(cutoff_low < dist2, dist2 < cutoff_high)
        
        dist_cutoff1 = dist1[cutoff_index1]
        dist_cutoff2 = dist2[cutoff_index2]
    
        return np.array((dist_cutoff1, cutoff_index1, dist_cutoff2, cutoff_index2))

    def _score_maths(self, dipole_map2, cutoff_low=0., cutoff_high=2.):
        '''
        Internal function to perform some maths that a) finds a list of nearest neighbours within some cutoff and then b) calculates the
        dot product between them

        :param dipole_map2: Dipole object for dipole map 2
        :param cutoff_low: Lower cutoff threashold to consider interactions from
        :param cutoff_high: Upper threashold to consider interactions within
        :returns: The alignment score of the first dipole map with respect to the second
        :returns: The alignment score of the second dipole map with respect to the first
        '''
            
        distance_list = self.dist_cutoff_list(dipole_map2.coordinates, cutoff_low = cutoff_low, cutoff_high=cutoff_high)
            
        score1 = []
        for i in range(np.shape(distance_list[1])[0]):
            Sc_tmp = np.dot(self.dipole[i], dipole_map2.dipole[distance_list[1][i]].T)
            score1.append(np.sum(Sc_tmp))
        score1 = filter(None, score1)

        score2 = []
        for i in range(np.shape(distance_list[3])[0]):
            Sc_tmp = np.dot(dipole_map2.dipole[i], self.dipole[distance_list[3][i]].T)
            score2.append(np.sum(Sc_tmp))
        score2 = filter(None, score2)
           
        return score1, score2


    def dipole_score(self, dipole_map2, cutoff_low = 0., cutoff_high = 2.):
        '''
        Return a score based on the alignment of the dipoles

        :param dipole_map2: The Dipole map (in jd.Dipole() form) that is being docked into self.
        :returns: An overall score for the dipole alignment.
        '''

        score1, score2 = self._score_maths(dipole_map2, cutoff_low = 0.0, cutoff_high = cutoff_high)
        
        Sc = - (np.median(score1) +  np.median(score2)) / 2
        return Sc

    def write_dipole(self, fname, radius = 0.3):
        '''
        Write a dipole map in a tcl format that can be read in by VMD
        
        :param fname: Name of output tcl file
        :param radius: Radius of drawn cone if desired (visualisation purposes only)
        '''

        write_file = open(fname, "w") # open a file for writing to

        # Get the end points for the dipole w.r.t. the base coordinates
        end_crd = self.coordinates + self.dipole

        for i in range(np.shape(self.coordinates)[0]):
            write_file.write("draw cone { %f %f %f } { %f %f %f } radius %f\n"%(self.coordinates[i][0], self.coordinates[i][1], self.coordinates[i][2], end_crd[i][0], end_crd[i][1], end_crd[i][2], radius))

    def total_dipole_calc(self):
        '''
        Calculate the total value of all the dipoles summed in a dipole map, a test function to see if this is minimised
        during binding

        :returns: An 'energy' corresponding to the overall sum
        '''

        dip_sum = np.sum(self.dipole, axis=0)
        energy = np.sqrt(dip_sum[0]**2 + dip_sum[1]**2 + dip_sum[2]**2)
        return energy

    def combine_dipole(self, dipole2):
        '''
        Combine two dipole maps, with vectors closer than 1 Ang adding, though scaled by distance

        :param dipole2: Second Dipole() structure
        :returns: A new dipole map with the combined properties of the two inputs.
        '''
        D = Dipole(d = np.vstack((self.dipole, dipole2.dipole)), c = np.vstack((self.coordinates, dipole2.coordinates)))
        dist = D.dist_cutoff_list(D.coordinates, cutoff_low = 0., cutoff_high = 1.)

        # Don't need to remove self-interacting terms as its == 0. exactly
        clash_index = np.where(np.sum(dist[1], axis=1) > 0.)

        d_cnt = 0
        dip_sum = []
        for c in clash_index[0]:
            contacts = np.sum(dist[1][c])
            #Reassign dipole value depending on contact
            for i in range(contacts):
                D.dipole[c] = (1 - dist[0][d_cnt]) * (D.dipole[c] + D.dipole[dist[1][c]][i])
                d_cnt += 1

        return D #score_after - score_before
