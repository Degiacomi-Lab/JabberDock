# Copyright (c) 2018 Lucas Rudden
#
# JabberDock is free software ;
# you can redistribute it and/or modify it under the terms of the GNU General Public License as published by the Free Software Foundation ;
# either version 2 of the License, or (at your option) any later version.
# JabberDock is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY ;
# without even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.
# See the GNU General Public License for more details.
# You should have received a copy of the GNU General Public License along with JabberDock;
# if not, write to the Free Software Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA 02111-1307 USA.
#
# Please reference this software and its authors in any works that make use of it
# Author : Lucas Rudden, lucas.rudden1@btinternet.com

import biobox as bb
import numpy as np
cimport numpy as np
cimport cython
from cpython cimport bool
import sys
from skimage import measure
from copy import deepcopy
import scipy.spatial as spaz
import time
from libc.math cimport sqrt
import JabberDock as jd
#import gpu_pairwise as gpupdist

class Jabber(object):
    '''
    A Jabber object contains information pertaining to an isometric surface that described an input density at a desired threashold - in the literature it is known as the Spatial and Temporal Influence Density (STID) map

    It contains information on the vertices, normals to said vertices, faces and the max value of the triangle used to define the isosurface using the marching cubes algorithm.
    The vertices returned are 3D coordinates
    '''

    def __init__(self, str dx, float cutoff = 0.5):
        '''
        We initalise our Jabber object, where information regarding the coordinates of an isometric surface, their faces, normals and maximum values for local topology
        are stored. We require a density object (bb.Density()) in the form of a dx map to interpret this.

        :param dx: The name of the dx file required for the generation of the Jabber isosurface
        '''

        super(Jabber, self).__init__()

        S = bb.Density()
        S._import_dx(dx)

        volume = S.properties['density']
        delta = S.properties['delta']
        origin = S.properties['origin']
        # get the complete transformation matrix based on the origin point of the coordinate system (as marching cubes assumes (0,0,0))

        # Define the step size in a normal coordinate system
        spacing = np.sqrt(delta[0, 0]**2 + delta[0, 1]**2 + delta[0, 2]**2)  # assuming the voxel is a cube (which it should be!)
        # Problem with the spacing is that it means we can't replicate the same score if we just load in our map!
        step_matrix = [spacing, spacing, spacing]

        # Define the coordinate points of the verties, the triangular faces (index of verts), the normal direction at each vertex, and the maximum value in the local triangle region
        verts, self.faces, self.norm, self.values = measure.marching_cubes(volume, level = cutoff, spacing = step_matrix)

        self.verts = verts + origin

    def rotate(self, np.ndarray COM = np.array((0., 0., 0.)), np.ndarray R = np.identity(3)):
        '''
        Rotate a jabber object (STID map) via a rotation matrix in place (presumably obtained through geometry.translate_pdb or geometry.translate_map)

        :param COM: Centre of Mass of system
        :param R: Rotation matrix
        '''

        self.translate(-COM)
        self.verts = np.matmul(self.verts, R)
        self.translate(COM)
        self.norm = np.matmul(self.norm, R)

    def translate(self, np.ndarray trans = np.array((0., 0., 0.))):
        '''
        Translate a jabber object (STID map) via a x, y, z vector in place

        :param trans: Transformation vector
        '''

        self.verts += trans

    def distance_list(self, jabber2, float cutoff=1.6):
        '''
        Get a list of indicies for the vertices that are closest to one another, and use a cutoff to remove points that we're not interested in
        (i.e. beyond a typical non-bonded cutoff distance)

        :param jabber2: A Jabber() object to get a distance list between (we're interested in jabber.verts here)
        :param cutoff: A cutoff to exclude interactions from - typically between 1 Ang. and 3 Ang. (non-bonded cutoff is what we're interested in)
        :returns: Numpy array containing two arrays - the first contains the distances to the closest point on points2 from points1 for each point on points1 where said distance is less than the cutoff. The second contains the minimum distances from points1 to points2 for each point in points2 where that distance is less than the cutoff. Note that this means that the two don't have to have the same shape.
        :returns: Numpy array containing two arrays. The first has shape A, the second shape B. Each element in an array contains a boolean referencing whether each point is able to make contact with a point on the other surface within the cutoff.
        :returns: Numpy array containing two arrays. The first has shape A, the second shape B. Each element in an array references the index to an element in the other array to which the point in the first array is closest to. 
        
        :Example: Say you have two sets of jabber objects, jabber_A and jabber_B, with jabber_A.verts and jabber_B.verts with shape 3 x 3 and 6 x 3 respectively (3 corresponding to x, y, z). Plugging them into here you might get the following arrays:
        dist = array([array([1.09]), array([0.58, 1.47])]).
    
        bool_index = array([array([True, False, False]), array([True, True, False, False, False, False])]).
    
        min_index = array([array([4, 4, 0]), array([2, 2, 2, 0, 0, 2])])
    
        dist[0] says that there is only 1 point in A within the cutoff of B. It has a distance 1.09 to that point. dist[1] shows that there are two points in contact with A.
    
        bool_index[0] essentially gives us the index on A to which dist[0] applies to. So only the first point is in range of B, equally with B only the first two points are in range with A.
    
        Finally, min_index[0] provides us with the elements on B that each point is closest to. So the first point in A is closest to point 4 (or index 5) on B, so is point 2. Point 3 on A is closest to the first point on B. 
    
        This information together gives us each set of points that satisfy our cutoff, and the distance between them.
        '''

        if isinstance(self, jd.Jabber) and isinstance(jabber2, jd.Jabber):
            points1 = self.verts
            points2 = jabber2.verts
        else:
            raise Exception("ERROR: One of your inputs is not a Jabber input, please check and try again.")
            exit(0)

        # Include GPU functionality
        #dist1 = gpupdist.pairwise_distance(verts1, verts2)
        dist1 = spaz.distance.cdist(points1, points2, 'euclidean')
    
        cdef np.ndarray[np.int64_t, ndim=1] min_index2 = jd.geometry.fast_argmin_axis_0(dist1)
    
        dist2 = spaz.distance.cdist(points2, points1, 'euclidean') # get verts2 nearest neighbours
        #dist2 = gpupdist.pairwise_distance(verts2, verts1)
        cdef np.ndarray[np.int64_t, ndim=1] min_index1 = jd.geometry.fast_argmin_axis_0(dist2)
        dist_min1 = []

        for i in range(dist1.shape[0]):
            dist_min1.append(dist1[i, min_index1[i]])
    
        dist_min2 = []
        for i in range(dist2.shape[0]):
            dist_min2.append(dist2[i, min_index2[i]])

        dist_min1 = np.array(dist_min1)
        dist_min2 = np.array(dist_min2)

        cutoff_index1 = dist_min1 < cutoff # Get the distances less than the cutoff
    
        cutoff_index2 = dist_min2 < cutoff

        # Get the smallest point to point distance AND those within the cutoff range
        cdef np.ndarray[np.npy_bool, ndim=1, cast=True] bool_index1 = np.logical_and(min_index1, cutoff_index1)
        # This gives us the boolean values for each index as to whether it satifies the criteria of closness and within the cutoff
        cdef np.ndarray[dtype=np.npy_bool, ndim=1, cast=True] bool_index2 = np.logical_and(min_index2, cutoff_index2)

        cdef np.ndarray[double, ndim=1] dist_cutoff1 = dist_min1[bool_index1]
        cdef np.ndarray[double, ndim=1] dist_cutoff2 = dist_min2[bool_index2]
    
        return np.array((dist_cutoff1, dist_cutoff2)), np.array((bool_index1, bool_index2)), np.array((min_index1, min_index2))
    
    #cpdef float scoring_function(np.ndarray protein_1_iso, np.ndarray protein_2_iso, np.ndarray dist, np.ndarray bool_index, np.ndarray index, float weight = 0.5):
    def scoring_function(self, jabber2, np.ndarray dist, np.ndarray bool_index, np.ndarray index, float weight = 0.5):
        '''
        Define a scoring function to obtain a value to obtain the 'goodness of a fit', or surface complementarity, between two isosurfaces at a specific 
        region of the surface. This relies on you first obtaining the set of nearest neighbours within a sensible cutoff with distance_list.

        :param jabber2: Second jabber object (STID map) which we're treating as the ligand
        :param dist: distance list containing two lists with distances to their nearest neighbours for proteins 1 and 2 respectively in the cutoff region. Can be created with jabber.distance_list()
        :param bool_index: Series of booleans list (for proteins 1 and 2) which are true when we have an index for an atom that satisfies the cutoff - we need this for access to normals etc. Can be created with jabber.distance_list()
        :param index: Minimum indicies for atoms that are within the cutoff distance (corresponds to the distance lists). Can be created with jabber.distance_list()
        :param weight: An arbiturary weighting, can be used to modify the distance weighting in the scoring function, and tends to give worse scores for less well fitting models. Taken from Lawrence & Colman 1993
        :returns: A singular score that judges the fit of the two STID maps. The larger, the better
        '''

        norm1_keep = deepcopy(self.norm) # Keep a copy as we need to use it later to find the corresponding nearest normals from protein 2

        # Now we need to calculate our score for the two isosurfaces given the input data
        # First we find the complementarity of protein 1 into protein 2
        norm1 = self.norm[bool_index[0]]
        correspond_norm2 = np.array(list(jabber2.norm[index[0]]))
        correspond_norm2 = correspond_norm2[bool_index[0]]

        # Now we flip the normals for protein 2 so they're pointing inward
        correspond_norm2 *= -1

        # Finally we calculate the surface complementarity using the two normals and the distance for Protein 1 into 2
        S12 = (np.einsum('ij,ij->i', norm1, correspond_norm2)) * np.exp(-1 * weight * (np.absolute(dist[0]))**2) * (float(np.sum(bool_index[0])) / float(np.shape(self.verts)[0]))

        # Next we find the surface complementary value for protein 2 into protein 1
        norm2 = jabber2.norm[bool_index[1]]
        correspond_norm1 = np.array(list(norm1_keep[index[1]]))
        correspond_norm1 = correspond_norm1[bool_index[1]]

        # Flip protein 1 normals
        correspond_norm1 *= -1

        ######## 
        # We could also modify the score by molecular overlap - a simple solution is if a point is closer to the centre of the opposite proteins molecular surface
        #######

        # Calculate the surface complementarity using the two normals and the distance for Protein 2 into 1
        S21 = (np.einsum('ij,ij->i', norm2, correspond_norm1)) * np.exp(-1 * weight * (np.absolute(dist[1]))**2) * (float(np.sum(bool_index[1])) / float(np.shape(jabber2.verts)[0]))
        
        # Finally we calculate the total surface complementarity given the two
        Sc = (np.median(S12) + np.median(S21)) / 2.

        return Sc

    def write_jabber(self, fname):
        '''
        Write the STID map object so it can be visualised. (Is written as a PDB file with each triangular vertex corner defined as an atom)

        :param fname: Name of output file (will be a .pdb file)
        '''

        S = bb.Structure(p=self.verts)
        S.write_pdb(fname)
