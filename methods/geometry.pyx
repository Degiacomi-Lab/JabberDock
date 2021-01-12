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
from libc.math cimport sqrt
from skimage import measure
from copy import deepcopy
import scipy.spatial as spaz

cpdef float vdw_energy(np.ndarray m1, np.ndarray m2):
    '''
    Calculate the vdw energy between two PDB structures (typically bb.Molecule()) - a check for clashes

    :param m1: Coordinates for structure 1
    :param m2: Coordinates for structure 2
    :returns: The calculated VdW scalar energy
    '''
    cdef float epsilon=1.0
    cdef float sigma=2.7
    cdef float cutoff=12.0
    cdef float energy=0.0
    cdef np.ndarray[float, ndim=1] d = np.empty()
    for i in xrange(0,len(m1),1):
        d.append(np.sqrt(np.sum((m2-m1[i])**2,axis=1)))
    dist=np.array(d)
    #detect interfacing atoms (atom couples at less than a certain cutoff distance
    couples=np.array(np.where(dist<cutoff)) #detect couples of clashing residues
    for i in xrange(0,len(couples[0]),1):
        d=dist[couples[0,i],couples[1,i]]
        energy+=4*epsilon*((sigma/d)**9-(sigma/d)**6)
    return energy

cpdef fast_argmin_axis_0(np.ndarray a):
    '''
    A function that basically takes the place of numpy.argmin(), but is much faster
            
    :param a: The array that we need the index for the minimum value
    :returns: An array of indicies into the array (with len(flattened a))
    '''

    matches = np.nonzero((a == np.min(a, axis=0)).ravel())[0]
    rows, cols = np.unravel_index(matches, (a.shape[0], a.shape[1]))
    argmin_array = np.empty(a.shape[1], dtype=np.intp)
    argmin_array[cols] = rows
    return np.asarray(argmin_array, dtype=int)

cpdef translate_pdb(pdb, float x = 0., float y = 0., float z = 0.):
    '''
    translate a whole PDB structure by a given amount in place. Function is defunct relative to biobox.molecule.translate() 

    :param pdb: PDB structure (i.e. bb.Molecule() object)
    :param x: The desired x translation amounts
    :param y: The desired y translation amounts
    :param z: The desired z translation amounts
    '''

    if 'center' not in pdb.properties:
        pdb.get_center()

    pdb.properties['center'][0] += x
    pdb.properties['center'][1] += y
    pdb.properties['center'][2] += z

    pdb.points[:, 0] += x
    pdb.points[:, 1] += y
    pdb.points[:, 2] += z 

cpdef rotate_pdb(pdb, float angle = 0., float axisx = 0., float axisy = 0., float axisz = 0., R = np.identity(3)):
    '''
    Method to rotate pdb structures about an axis of rotation from the origin in place

    :param pdb: PDB structure (i.e. bb.Molecule() object)
    :param angle: Scalar for angle to rotate (units of degrees)
    :param axisx: Point in x for the axis of rotation
    :param axisy: Point in y for the axis of rotation
    :param axisz: Point in z for the axis of rotation
    :param R: Rotation matrix if this function has already been applied
    :returns: The rotation matrix used to rotate the protein 
    '''
    
    # First center molecule before applying rotation
    if 'center' not in pdb.properties:
        pdb.get_center()
    pdb_center = deepcopy(pdb.properties['center'])
    pdb.translate(x = - pdb_center[0], y = - pdb_center[1], z = -pdb_center[2])

    # Find R if R not given
    if np.array_equal(R, np.identity(3)):
        # Convert angle to rads
        angle = angle * np.pi / 180.
        axis = [axisx, axisy, axisz]
   
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
        R = [[t*x*x + c, t*x*y - z*s, t*x*z + y*s], [t*x*y + z*s, 	t*y*y + c , t*y*z - x*s], [t*x*z - y*s, t*y*z + x*s, t*z*z + c]]

    # Finally, we multiply our coordinate system matrix by our rotation matrix - rotating the pdb structure
    vals = np.matmul(pdb.coordinates[0], R)
    pdb.set_xyz(vals)
 
    # Move back to original position
    pdb.translate(x = pdb_center[0], y = pdb_center[1], z = pdb_center[2])

    # Return the rotation m1= self.data.M.atomselect("*", "*", self.params.atoms)
    return np.asarray(R, dtype=float)

cpdef translate_map(dx, float x = 0., float y = 0., float z = 0.):
    '''
    Method to translate electron density maps based on origin point of grid in place

    :param dx: Density map previously imported from dx format (i.e. bb.Density() object)
    :param x: x coordianate to shift by (units in Ang.)
    :param y: y coordianate to shift by    
    :param z: z coordianate to shift by
    '''
    trans = [x, y, z]
    dx.properties['origin'] = dx.properties['origin'] + trans

cpdef rotate_map(dx, np.ndarray COM = np.array((0., 0., 0.)), float angle = 0., float axisx = 0., float axisy = 0., float axisz =  0., R = np.identity(3)):
    '''
    Method to rotate electron density maps about the an axis of rotation from the origin in place

    :param dx: Density map previously imported from dx format (i.e. bb.Density() object)
    :param COM: Centre of Mass of corresponding PDB (if using PDB as basis for roto-translations)
    :param angle: Angle to shift by (must be in deg)
    :param axisx: Point in x for the axis of rotation
    :param axisy: Point in y for the axis of rotation
    :param axisz: Point in z for the axis of rotation
    :param R: Rotation matrix if already defined
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

    # Now we need to shift the entire box so that the centre of the box is at the origin of the system, and we can
    # apply our rotation matrix to the origin
    dx_center = dx.properties['origin']
    delta = dx.properties['delta']

    # Move origin to coordinate with input PDB Centre of Mass
    translate_map(dx, -COM[0], -COM[1], -COM[2])
    rot = np.matmul(delta, R)

    dx.properties['delta'] = rot
    new_origin = dx.properties['origin']

    origin_rot = np.matmul(new_origin, R)
    dx.properties['origin'] = origin_rot

    # Reshift the origin back to its origional position
    translate_map(dx, COM[0], COM[1], COM[2])   
 
    # Return the rotation matrix
    return np.asarray(R, dtype=float)

cpdef redefine_coord(dx, float cutoff=0.5):
    '''
    Define a coordinate system given a cutoff for the isosurface based on the input information of the dx map. The function is also 
    avaliable to rotate / translate maps as necessary, providing a map of points at a spectific isovalue cutoff. 
    The rotation already occurs given a specific delta in the input map.

    :param dx: Density map previously imported from dx format  (i.e. bb.Density() object)
    :param cutoff: Isolevel cutoff for map
    :returns: An array which contains the coordinates of every voxel point greater than the isovalue cutoff given
    '''

    delta = dx.properties['delta']
    origin = dx.properties['origin']
    dens = dx.properties['density']

    # Define a new coordinate system based on points and delta
    pts = []
    for x in range(dens.shape[0]):
       for y in range(dens.shape[1]):
          for z in range(dens.shape[2]):

              if dens[x,y,z]>=cutoff:
                  newx = delta[0, :]*x
                  newy = delta[1, :]*y
                  newz = delta[2, :]*z
                  p = newx + newy + newz + origin
                  pts.append(p)

    points = np.array(pts)
    return points

cpdef define_isosurface(dx, np.ndarray COM = np.array((0., 0., 0.)), float cutoff = 0.5, np.ndarray R = np.identity(3)):
    '''
    Obtain an isosurface of a density object that contains the vertices, normals to said vertices,
    faces and the max value of the triangle used to define the isosurface using the marching cubes algorithm. The vertices returned are 3D coordinates.
    This is defunct now, as when defining a jabber object this is performed automatically, and the isosurface is an associated property of the jabber object.

    :param dx: Density map previously imported from _dx format
    :param COM: Centre of mass of the density map. Can be the same as get_center() of corresponding PDB structure.
    :param cutoff: Isolevel cutoff for map
    :param R: Rotation matrix used to define rotation (can be output from rotate_map)
    :returns: Four arrays - verts, faces, normals and values. Verts contains the vertices of the triangular mesh that makes up the isosurface (V, 3) is the shape, which match the order of the input volume. V corresponds to the unique vertices. Faces contains the triangular faces that refer to verts, (F, 3) is the shape - where F is the number of faces. Normals are the normal direction to each vertex with shape (V, 3). The values contains a measure for the maximum value of the STID quantity nead each vertex (useful for visualisation), with shape (V).
    '''
    
    volume = dx.properties['density']
    delta = dx.properties['delta']
    origin = dx.properties['origin']
    # get the complete transformation matrix based on the origin point of the coordinate system (as marching cubes assumes (0,0,0))
    
    # Define the step size in a normal coordinate system
    spacing = np.sqrt(delta[0, 0]**2 + delta[0, 1]**2 + delta[0, 2]**2)  # assuming the voxel is a cube (which it should be!)
    # Problem with the spacing is that it means we can't replicate the same score if we just load in our map!
    step_matrix = [spacing, spacing, spacing]
    
    # Define the coordinate points of the verties, the triangular faces (index of verts), the normal direction at each vertex, and the maximum value in the local triangle region 
    verts, faces, normals, values = measure.marching_cubes_lewiner(volume, level = cutoff, spacing = step_matrix)
    newverts = np.matmul(verts, R) + origin # roto-translate our coordinate points as necessary
    return np.asarray(newverts), np.asarray(faces), np.asarray(normals), np.asarray(values)

cpdef get_mol_center(np.ndarray points):
    '''
    Get the geometrix center of a molecular surface, as defined by the vertix coordinates, primarily an internal function

    :param points: Input coordinates of molecular surface (in n by 3 shape), where n is the number of points
    :returns: A Cartesian vector containing the coordinates of the geometric centre.
    '''

    return np.mean(points, axis = 0)

cpdef get_center_distance(np.ndarray points):
    '''
    Return the distance (in angstroems) from a set of points from the center of a molecular surface

    :param points: Input coordinates of molecular surface (n by 3 in shape), where n is the number of points
    :returns: A Cartesian vector containing the distance of between the series of points input and the geometric centre of the STID map.
    '''

    center = get_mol_center(points)

    return points - center

cpdef distance_list(points1, points2, float cutoff=1.5):
    '''
    Get a list of indicies for the vertices that are closest to one another, and use a cutoff to remove points that we're not interested in
    (i.e. beyond a typical non-bonded cutoff distance)
    This function has the possibility of using the gpu pairwaise distance function (if uncommented) - though this is unstable for large proteins)

    :param points1: First set of coordinates, with shape A x 3
    :param points2: Second set of coordinates, with shape B x 3
    :param cutoff: A cutoff to exclude interactions from
    :returns: Numpy array containing two arrays - the first contains the distances to the closest point on points2 from points1 for each point on points1 where said distance is less than the cutoff. The second contains the minimum distances from points1 to points2 for each point in points2 where that distance is less than the cutoff. Note that this means that the two don't have to have the same shape.
    :returns: Numpy array containing two arrays. The first has shape A, the second shape B. Each element in an array contains a boolean referencing whether each point is able to make contact with a point on the other surface within the cutoff.
    :returns: Numpy array containing two arrays. The first has shape A, the second shape B. Each element in an array references the index to an element in the other array to which the point in the first array is closest to. 
    
    :Example: Say you have two sets of coordinates, A and B, with shape 3 x 3 and 6 x 3 (3 corresponding to x, y, z). Plugging them into here you might get the following arrays:
    dist = array([array([1.09]), array([0.58, 1.47])]).

    bool_index = array([array([True, False, False]), array([True, True, False, False, False, False])]).

    min_index = array([array([4, 4, 0]), array([2, 2, 2, 0, 0, 2])])

    dist[0] says that there is only 1 point in A within the cutoff of B. It has a distance 1.09 to that point. dist[1] shows that there are two points in contact with A.

    bool_index[0] essentially gives us the index on A to which dist[0] applies to. So only the first point is in range of B, equally with B only the first two points are in range with A.

    Finally, min_index[0] provides us with the elements on B that each point is closest to. So the first point in A is closest to point 4 (or index 5) on B, so is point 2. Point 3 on A is closest to the first point on B. 

    This information together gives us each set of points that satisfy our cutoff, and the distance between them.
    '''

    # Include GPU functionality
    #dist1 = gpupdist.pairwise_distance(verts1, verts2)
    dist1 = spaz.distance.cdist(points1, points2, 'euclidean')
    
    cdef np.ndarray[np.int64_t, ndim=1] min_index2 = fast_argmin_axis_0(dist1)
    
    dist2 = spaz.distance.cdist(points2, points1, 'euclidean') # get verts2 nearest neighbours
    #dist2 = gpupdist.pairwise_distance(verts2, verts1)
    cdef np.ndarray[np.int64_t, ndim=1] min_index1 = fast_argmin_axis_0(dist2)
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

cpdef get_minmax_crd(fname_rec, fname_times = 2, fname_buff = 0.0):
    ''' 
    Obtain the min max coordinates from a receptors conformational space to generate sensible structures for ensembles in JabberDock.
    We move the ligand around this space to represent the size the receptor occupies, but we do want it to be slightly bigger.

    :param fname_rec: PDB file name of receptor
    :param fname_times: Number greater than 1 to represent extra space added to search space (default is 2, so twice as big as the receptor)
    :param fname_buff: Buffer region at the edges to include (default is 0.0 ang.)
    :returns: numpy array size 3, containing the x, y and z lengths that provide the sample space about which the ligand will move in.
    '''
 
    receptor = bb.Molecule()
    receptor.import_pdb(fname_rec)

    rec_minx = np.min(receptor.points[:, 0])
    rec_miny = np.min(receptor.points[:, 1])
    rec_minz = np.min(receptor.points[:, 2])

    rec_maxx = np.max(receptor.points[:, 0])
    rec_maxy = np.max(receptor.points[:, 1])
    rec_maxz = np.max(receptor.points[:, 2])

    rec_x = rec_maxx - rec_minx
    rec_y = rec_maxy - rec_miny
    rec_z = rec_maxz - rec_minz

    rec_x *= fname_times
    rec_y *= fname_times
    rec_z *= fname_times

    rec_x += fname_buff
    rec_y += fname_buff
    rec_z += fname_buff

    return np.array((rec_x, rec_y, rec_z))

