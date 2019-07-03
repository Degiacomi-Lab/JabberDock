# Copyright (c) 2014-2017 Matteo Degiacomi
#
# BiobOx is free software ;
# you can redistribute it and/or modify it under the terms of the GNU General Public License as published by the Free Software Foundation ;
# either version 2 of the License, or (at your option) any later version.
# BiobOx is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY ;
# without even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.
# See the GNU General Public License for more details.
# You should have received a copy of the GNU General Public License along with BiobOx ;
# if not, write to the Free Software Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA 02111-1307 USA.
#
# Author : Matteo Degiacomi, matteothomas.degiacomi@gmail.com

import os

import scipy.ndimage.filters
from sklearn.cluster import DBSCAN
import numpy as np
import pandas as pd

from biobox.classes.structure import Structure

class Density(Structure):
    '''
    Subclass of :func:`Structure <structure.Structure>`, allows importing density map, and transform them in a PDB file containing a collection of spheres placed on the map's high density regions.
    '''

    def __init__(self):
        '''
        A density map is fully described by the following attributes, stored in the self.properties dictionary:
        
        :param density: density map
        :param delta:   scaling factor for voxels (default is [1, 1, 1] Angstrom)
        :param size:    dimensions in voxels
        :param origin:  bottom-left-front corner of the cube
        :param radius: radius of points composing the density map
        :param format: format name (only dx supported at the moment)
        '''

        super(Density, self).__init__()
        self._reset_info()

    def _reset_info(self, r=1.9):
        '''
        reset all properties related to a density map (clean)
        '''

        self.properties['density'] = np.array([])  # density map

        self.properties['radius'] = r
        # scaling factor for voxels (default is [1, 1, 1] Angstrom)
        self.properties['delta'] = np.array([])
        self.properties['size'] = np.array([])  # dimensions in voxels
        # bottom-left-front corner of the cube
        self.properties['origin'] = np.array([])
        self.properties['format'] = ""  # file format
        self.properties['scan'] = np.array([])
        self.properties['filename'] = ""

        self.clear()

    def return_density_map(self):
        '''
        :returns: density map as 3D numpy array
        '''
        return self.properties['density']

    def import_map(self, filename, fileformat='dx'):
        '''
        Import density map and fill up the points and properties data structures.

        :param  filename: name of density file to load
        :param  fileformat: at the moment supports dx, ccp4, mrc and imod
        '''

        if not os.path.exists(filename):
            raise Exception("%s not found!" % filename)

        # call format-specific loading functions.
        # function should fill up all required properties in _reset_info(), and
        # load the map as a 3D array, containing intensity values.
        try:
            if fileformat == 'dx':
                self._import_dx(filename)
            elif fileformat == 'ccp4':
                self._import_mrc(filename, 'ccp4')
            elif fileformat == 'mrc':
                self._import_mrc(filename, 'mrc')
            elif fileformat == 'imod':
                self._import_mrc(filename, 'imod')
            else:
                raise Exception("sorry, format %s is not supported" % fileformat)

        except Exception as e:
            self._reset_info()
            Exception("ERROR: %s" % e)

        # if any error went undetected during loading (missing information),
        # data structures may be inconsistent. Call cleaning procedure!
        if len(self.properties['density']) == 0:
            print("density map could not be correctly loaded!")
            self._reset_info(self)
        elif len(self.properties['size']) == 0:
            print("density map information missing!")
            self._reset_info(self)
        elif len(self.properties['origin']) == 0:
            print("map origin information missing!")
            self._reset_info(self)
        elif len(self.properties['delta']) == 0:
            print("voxel size information missing!")
            self._reset_info(self)

        # if all required information is present, place points instead of
        # voxels
        else:
            try:
                self.place_points()
            except Exception as e:
                pass

            self.properties['format'] = format
            self.properties['filename'] = filename

        self.properties["sigma"] = np.std(self.properties['density'])

    def import_numpy(self, data, origin=[0, 0, 0], delta=np.identity(3)):
        '''
        import a numpy 3D array to allow manipulation as a density map

        :param data: numpy 3D array
        :param origin: coordinates of bottom left corner of the map
        :param delta: voxels' shape (default is a cubic voxel of 1 Angstrom-long sides).
        '''
                
        if len(data.shape) != 3:
            raise Exception("ERROR: a 3D numpy array is expected")
        
        self.properties['density'] = data
        self.properties['origin'] = np.array(origin)
        self.properties['size'] = np.array(data.shape)
        self.properties['delta'] = delta

        # sphere size corresponding to the volume of one voxel
        voxel_volume = self.properties['delta'][0, 0] * self.properties['delta'][1, 1] * self.properties['delta'][2, 2]
        self.properties['radius'] = (voxel_volume * 3 / (4 * np.pi))**(1 / 3.0)

    def place_points(self, sigma=0, noise_filter=0.01):
        '''
        given density information, place points on every voxel above a given threshold.

        :param sigma: intensity threshold value.
        :param noise_filter: launch DBSCAN clustering algorithm to detect connected regions in density map. Regions representing less than noise_filter of the total will be removed. This is a ratio, value should be between 0 and 1.
        '''

        thresh = self.get_thresh_from_sigma(sigma)

        if not np.any(self.properties['density'] > thresh):
            raise IOError("selected threshold leads to empty point ensemble")

        if noise_filter >= 1 or noise_filter < 0:
            raise IOError("noise_filter should be between 0 and 1")

        # define scaling to shrink everything by a size equal to spheres radius
        size = np.max(np.transpose(np.where(self.properties['density'] > thresh)), axis=0) - np.min(np.transpose(np.where(self.properties['density'] > thresh)), axis=0)
        scaled_size = np.max(np.transpose(np.where(self.properties['density'] > thresh)) * np.diag(self.properties['delta']), axis=0)-np.min(np.transpose(np.where(self.properties['density'] > thresh)) * np.diag(self.properties['delta']), axis=0)

        if self.properties['delta'][0, 0] != 1:
            scale_x = (self.properties['delta'][0, 0] - 1.0) / (scaled_size[0] - size[0]) * (scaled_size[0] - self.properties['radius'])
        else:
            scale_x = 1.0

        if self.properties['delta'][1, 1] != 1:
            scale_y = (self.properties['delta'][1, 1] - 1.0) / (scaled_size[1] - size[1]) * (scaled_size[1] - self.properties['radius'])
        else:
            scale_y = 1.0

        if self.properties['delta'][2, 2] != 1:
            scale_z = (self.properties['delta'][2, 2] - 1.0) / (scaled_size[2] - size[2]) * (scaled_size[2] - self.properties['radius'])
        else:
            scale_z = 1
        # create structure data (ensemble of points and their center, and store
        # its geometric center)
        points = np.transpose(np.where(self.properties['density'] > thresh)) * np.array([scale_x, scale_y, scale_z]) + self.properties['origin'] + np.ones(3) * self.properties['radius']

        # remove previous points arrangement, and create new one (necessary,
        # since the amount of points will change, and cannot therefore be
        # considered a new conformation)
        self.clear()

        # apply DBSCAN noise filter
        if noise_filter != 0:
            step = self.properties['delta'][0, 0] * np.sqrt(3)
            db = DBSCAN(eps=step, min_samples=10).fit(points)
            pts2 = []
            for i in np.unique(db.labels_):
                if np.sum(db.labels_ == i) / float(len(points)) > 0.01 and i != -1:
                    if len(pts2) == 0:
                        pts2 = points[db.labels_ == i]
                    else:
                        pts2 = np.concatenate((pts2, points[db.labels_ == i]))

            self.add_xyz(pts2)


        else:
            self.add_xyz(points)

        #update radii list (useful for instance for CCS calculation)
        idx = np.arange(len(self.points))
        self.data = pd.DataFrame(idx, index=idx, columns=["radius"])

        # self.get_center()

    def write_dx(self, fname='dens.dx'):
        '''
        Write a density map in DX format

        :param fname: output file name
        '''

        dens = self.properties['density']
        origin = self.properties['origin']
        delta = self.properties['delta']

        fout = open(fname, "w")
        fout.write("# density generated with SBT\n#\n#\n#\n")
        fout.write("object 1 class gridpositions counts %s %s %s\n" % (dens.shape[0], dens.shape[1], dens.shape[2]))
        fout.write("origin %s %s %s\n"%(origin[0], origin[1], origin[2]))
        fout.write("delta %s %s %s\n"%(delta[0, 0], delta[0, 1], delta[0, 2]))
        fout.write("delta %s %s %s\n" % (delta[1, 0], delta[1, 1], delta[1, 2]))
        fout.write("delta %s %s %s\n" %(delta[2, 0], delta[2, 1], delta[2, 2]))
        fout.write("object 2 class gridconnections counts %s %s %s\n"%(dens.shape[0], dens.shape[1], dens.shape[2]))
        fout.write("object 3 class array type double rank 0 items %i data follows\n"%(dens.shape[0] * dens.shape[1] * dens.shape[2]))

        cnt = 0
        for xpos in range(0, dens.shape[0], 1):
            for ypos in range(0, dens.shape[1], 1):
                for zpos in range(0, dens.shape[2], 1):
                    fout.write("%s " % dens[xpos, ypos, zpos])
                    cnt += 1
                    if cnt % 3 == 0:
                        fout.write("\n")

        fout.close()

    def export_as_pdb(self, outname, step, threshold=0.1):
        '''
        Write a pdb file with points where the density exceeds a threshold

        :param outname: output file name
        :param step: stepsize used to generate the density map
        :param threshold: density to be exceeded to generate a point in pdb
        '''

        # @todo assign spheres beta factor to associated density value

        dens = self.properties['density']
        origin = self.properties['origin']

        fout = open(outname, 'w')

        cnt = 1
        identifier = 'ATOM'  # atom to be used to mimick density
        symbol = 'H'  # element for atom

        print('exporting density greater than %s to pdb' % threshold)

        for xpos in range(0, dens.shape[0], 1):
            for ypos in range(0, dens.shape[1], 1):
                for zpos in range(0, dens.shape[2], 1):
                    if dens[xpos, ypos, zpos] > threshold:
                        x_coord = float(xpos / step) + origin[0]
                        y_coord = float(ypos / step) + origin[1]
                        z_coord = float(zpos / step) + origin[2]
                        L = '%-6s%5s  %-4s%-4s  DIS    %8.3f%8.3f%8.3f  1.00  0.00            \n'%(identifier, cnt, symbol, symbol, x_coord, y_coord, z_coord)
                        fout.write(L)
                        cnt += 1

        fout.close()

    def _import_dx(self, filename):
        '''
        import density map and fill up the points and properties data structures.

        :param filename: name of dx file to load
        '''

        try:
            fin = open(filename, "r")
        except Exception as e:
            raise Exception('opening of file %s failed!' % filename)

        d = []
        dlt = []
        for line in fin:
            w = line.split()
            if w[0] == "#":
                continue
            elif len(w) <= 3: #and np.array(list(w)).dtype == ('float'):
                try:
                    w = np.array(w).astype(float)
                    for i in range(len(w)):
                        d.append(w[i])
                except Exception as e:
                    continue
                # get coordinates
            elif len(w) > 2 and w[0] == "object" and w[3] == "gridpositions":
                self.properties['size'] = np.array([w[-3], w[-2], w[-1]]).astype(int)
            # get scaling factor
            elif len(w) > 2 and w[0] == "delta":
                dlt.append(w[1:4])
            # get position of origin
            elif len(w) > 2 and w[0] == "origin":
                self.properties['origin'] = np.array(w[1:4]).astype(float)


        # scaling factor with respect of unit cell voxels
        self.properties['delta'] = np.array(dlt).astype(float)

        try:
            self.properties['density'] = np.reshape(
                np.array(d).astype(float),
                (self.properties['size'][0],
                 self.properties['size'][1],
                 self.properties['size'][2]))
        except Exception as ex:
            raise Exception("reshaping of dx data failed! Dimensions and dataset size are inconsistent!")

    def _import_mrc(self, filename, fileformat):
        '''
        import density map in MRC, CCP4 or IMOD format.

        MRC format here: www2.mrc-lmb.cam.ac.uk/image2000.html

        CCP4 format here: www.ccp4.ac.uk/html/maplib.html

        :param filename name of MRC or CCP4 file to load
        :param fileformat: can be mrc, imod or ccp4
        '''

        import biobox.classes.density_MRC as MRC
        try:
            [density, data] = MRC.read_density(filename, fileformat)
            self.properties['density'] = density
            self.properties['origin'] = np.array(data.origin)
            self.properties['size'] = np.array(density.shape)
            self.properties['delta'] = np.identity(3) * data.mrc_data.data_step

            # sphere size corresponding to the volume of one voxel
            voxel_volume = self.properties['delta'][0, 0] * self.properties['delta'][1, 1] * self.properties['delta'][2, 2]
            self.properties['radius'] = (voxel_volume * 3 / (4 * np.pi))**(1 / 3.0)

        except Exception as e:
            raise Exception("%s" % e)

    def get_volume(self):
        '''
        compute density map volume. This is done by counting the points, and multiplying that by voxels' volume.

        .. warning:: can be called only after :func:`place_points <density.Density.place_points>` has been called.

        .. warning:: supposes unskewed voxels.
        '''
        return self.properties['delta'][0, 0] * self.properties['delta'][1, 1] * self.properties['delta'][2, 2] * len(self.points)


