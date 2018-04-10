from Result import Result

import os.path as op
import glob
import subprocess
import time
import matplotlib.pyplot as plt

from collections import OrderedDict as odict
from shapely.geometry import Polygon
from descartes.patch import PolygonPatch
from matplotlib.collections import PatchCollection

from tools import *

import mesh as meshmod

class Marginal(Result):
    def __init__(self, io, nodename, vn=100, wn=100):
        super(Marginal, self).__init__(io, nodename)
        self.path = op.join(self.io.output_directory,
                            self.nodename + '_marginal_density')
        self.data_path = op.join(self.io.output_directory, 'marginal_density.npz')
        self.projfname = self.modelpath.replace('.model', '.projection')
        self.vn, self.wn = vn, wn

    def __getitem__(self, name):
        return self.density[name]

    @property
    def density(self):
        # Get the projection file
        self.read_projection()
        # If there's no new projection and we've already calculated everything,
        # just return the existing calculated marginals
        if op.exists(self.data_path) and not self.new_projection:
            load_data = np.load(self.data_path)['data'][()]
            if self.modelname in load_data:
                return load_data[self.modelname]

        # Initialise the marginal densities for each frame of the simulation
        v = np.zeros((len(self.times), self.projection['N_V']))
        w = np.zeros((len(self.times), self.projection['N_W']))

        # Load in the masses from the densities
        masses, coords_ = [], None
        for ii, fname in enumerate(self.fnames):
            density, coords = read_density(fname)
            if coords_ is None:
                coords_ = coords
            else:
                assert coords == coords_
            masses.append(calc_mass(self.mesh, density, coords))
        masses = np.vstack(masses)
        assert masses.shape[0] == len(self.fnames)

        # Calculate the merginals for each frame and store in 'data'
        v, w, bins_v, bins_w = self.calc_marginal_density(
            v, w, masses, coords, self.projection)
        data = {'v': v, 'w': w, 'bins_v': bins_v,
                'bins_w': bins_w, 'times': self.times}

        # Save 'data' into a compressed file
        if op.exists(self.data_path):
            other = np.load(self.data_path)['data'][()]
            other.update({self.modelname: data})
            save_data = other
        else:
            save_data = {self.modelname: data}
        np.savez(self.data_path, data=save_data)
        return data

    # Using the projection file, the mass from in each cell of the mesh
    # (for each frame of the simulation)
    # can be summed into the two marginal densities
    # note, v and w are 2D matrices of the marginal density bin values for each
    # frame of the simulation
    def calc_marginal_density(self, v, w, masses, coords, proj):
        # temp function 'scale' to parse a transition row
        # in the projection file
        def scale(var, proj, mass):
            bins = [marg.split(',') for marg in proj.split(';')
                    if len(marg) > 0]
            for jj, dd in bins:
                var[:, int(jj)] += mass * float(dd)
            return var

        # Each cell in the mesh has a transition row in the projection file
        for trans in proj['transitions'].findall('cell'):
            # Get the coordinates of this cell and its mass each time
            i, j = [int(a) for a in trans.find('coordinates').text.split(',')]
            cell_mass = masses[:, coords.index((i, j))]

            # Calculate and add the density values for each bin
            v = scale(v, trans.find('vbins').text, cell_mass)
            w = scale(w, trans.find('wbins').text, cell_mass)

        # Generate the linspace values for plotting
        bins_v = np.linspace(proj['V_min'], proj['V_max'], proj['N_V'])
        bins_w = np.linspace(proj['W_min'], proj['W_max'], proj['N_W'])
        return v, w, bins_v, bins_w

    def make_projection_file(self):
        # Run the projection app to analyse the model file and get
        # dimensions
        projection_exe = op.join(getMiindAppsPath(), 'Projection', 'Projection')
        out = subprocess.check_output(
          [projection_exe, self.modelfname], cwd=self.io.xml_location)

        # Parse the output
        vmax, wmax = np.array(out.split('\n')[3].split(' ')[2:], dtype=float)
        vmin, wmin = np.array(out.split('\n')[4].split(' ')[2:], dtype=float)
        # assert that we bound the range
        inc = lambda x: x * 1.01 if x > 0 else x * 0.99
        vmax, wmax = inc(vmax), inc(wmax)
        vmin, wmin = -inc(-vmin) , -inc(-wmin)

        # Run the projection app to generate the .projection file
        cmd = [projection_exe, self.modelfname, vmin, vmax,
               self.vn, wmin, wmax, self.wn]
        subprocess.call([str(c) for c in cmd], cwd=self.io.xml_location)

    def read_projection(self):
        proj_pathname = op.join(self.io.xml_location, self.projfname)
        self.new_projection = False
        # Does the pojection file exist? If not, generate it.
        if not op.exists(proj_pathname):
            print('No projection file found, generating...')
            self.make_projection_file()
            self.new_projection = True

        # If we don't have projection data loaded, load it!
        if not hasattr(self, 'projection'):
            with open(proj_pathname) as proj_file:
                self._proj = ET.fromstring(proj_file.read())

        # Has the projection file changed since we last loaded? If so, reload.
        if (int(self._proj.find('W_limit/N_W').text) != self.wn or
            int(self._proj.find('V_limit/N_V').text) != self.vn):
            print('New N in bins, generating projection file...')
            self.make_projection_file()
            with open(proj_pathname) as proj_file:
                self._proj = ET.fromstring(proj_file.read())
            self.new_projection = True

        self.projection =  {
            'transitions': self._proj.find('transitions'),
            'V_min': float(self._proj.find('V_limit/V_min').text),
            'V_max': float(self._proj.find('V_limit/V_max').text),
            'N_V': int(self._proj.find('V_limit/N_V').text),
            'W_min': float(self._proj.find('W_limit/W_min').text),
            'W_max': float(self._proj.find('W_limit/W_max').text),
            'N_W': int(self._proj.find('W_limit/N_W').text),
        }

    def plotV(self, time, ax=None):
        if not ax:
            fig, ax = plt.subplots()
            ax.plot(self['bins_v'], self['v'][self['times'].index(time), :])
            fig.show()
        else:
            ax.plot(self['bins_v'], self['v'][self['times'].index(time), :])

    def plotW(self, time, ax=None):
        if not ax:
            fig, ax = plt.subplots()
            ax.plot(self['bins_w'], self['w'][self['times'].index(time), :])
            fig.show()
        else:
            ax.plot(self['bins_w'], self['w'][self['times'].index(time), :])

    def generatePlotImages(self, image_size=300):
        if not op.exists(self.path):
            os.mkdir(self.path)
        for ii in range(len(self['times'])):
            fig, axs = plt.subplots(1, 2)
            plt.suptitle('time = {}'.format(self['times'][ii]))
            self.plotV(self['times'][ii], axs[0])
            self.plotW(self['times'][ii], axs[1])
            figname = op.join(self.path,
                              '{}_'.format(ii) +
                              '{}.png'.format(self['times'][ii]))
            fig.savefig(figname, res=image_size, bbox_inches='tight')
            plt.close(fig)

    def generateMarginalAnimation(self, filename, image_size=300, generate_images=True, time_scale=1.0):
        # Generate the density image files
        if generate_images:
            self.generatePlotImages(image_size)

        try:
            # grab all the filenames
            files = glob.glob(op.join(self.path, '*.png'))
            files.sort()

            # calculate duration of each frame - this is the time between each
            # image
            durations = [((self.times[0])*time_scale)/1000.0]
            for t in range(len(self.times)-1):
                durations.append(((self.times[t+1] - self.times[t])*time_scale)/1000.0)

            # Generate an image list file with the calculated durations
            with open(op.join(self.path, 'filelist.txt'), 'w') as lst:
                d = 0
                for f in files:
                    lst.write('file \'' + f + '\'\n')
                    lst.write('duration ' + str(durations[d]) + '\n')
                    d += 1

            # note ffmpeg must be installed
            process = ['ffmpeg',
                '-f', 'concat',
                '-safe', '0',
                '-i', op.join(self.path, 'filelist.txt')]

            process.append(filename + '.mp4')

            subprocess.call(process)
        except OSError as e:
            if e.errno == 2:
                print "MIIND API Error : generateMarginalAnimation() requires ffmpeg to be installed."
            else:
                print "MIIND API Error : Unknown Error"
                print e