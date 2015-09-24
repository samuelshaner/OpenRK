__author__ = 'Sam Shaner'
__email__ = 'shaner@mit.edu'

# Import modules
import matplotlib
import openrk as rk
import checkvalue as cv

# force headless backend, or set 'backend' to 'Agg'
# in your ~/.matplotlib/matplotlibrc
matplotlib.use('Agg')

import matplotlib.pyplot as plt

# Force non-interactive mode, or set 'interactive' to False
# in your ~/.matplotlib/matplotlibrc
plt.ioff()

import numpy as np
import numpy.random
import os

## A static variable for the output directory in which to save plots
SUBDIRECTORY = "plots/"
TINY_MOVE = 1.0e-8

def plot_flux(solver, plane='xy', offset=0., time=3, energy_groups=[0], 
              gridsize=250, name='flux', xlim=None, ylim=None, zlim=None):
    
    global SUBDIRECTORY

    # Make directory if it does not exist
    if not os.path.exists(SUBDIRECTORY):
        os.makedirs(SUBDIRECTORY)

    # Error checking
    if not isinstance(solver, rk.SolverDiffusion):
        msg = 'Unable to plot the flux since input was not ' \
              'a SolverDiffusion class object'
        raise ValueError(msg)

    if not cv.is_integer(gridsize):
        msg = 'Unable to plot the flux since the gridsize {0} ' \
              'is not an integer'.format(gridsize)
        raise ValueError(msg)

    if plane not in ['xy', 'xz', 'yz']:
        msg = 'Unable to plot the flux with an invalid ' \
            'plane {0}. Plane options xy, xz, yz'.format(plane)
        raise ValueError(msg)

    if gridsize <= 0:
        msg = 'Unable to plot the flux with a negative ' \
              'gridsize {0}'.format(gridsize)
        raise ValueError(msg)

    if not isinstance(energy_groups, list):
        energy_groups = [energy_groups]

    if not cv.is_float(offset):
        msg = 'Unable to plot the flux since the offset {0} ' \
            'is not a float'.format(offset)
        raise ValueError(msg)

    print('Plotting the fluxes...')

    # Initialize a numpy array for the groupwise scalar fluxes
    fluxes = numpy.zeros((len(energy_groups), gridsize, gridsize))

    geometry = solver.getGeometryDiffusion()
    coords = get_pixel_coords(geometry, plane, offset, gridsize, xlim, ylim, zlim)

    for i in range(gridsize):
        for j in range(gridsize):

            if plane == 'xy':
                cell = geometry.findCell(x=coords['x'][i], y=coords['y'][j], z=offset)
            elif plane == 'xz':
                cell = geometry.findCell(x=coords['x'][i], y=offset, z=coords['z'][j])
            else:
                cell = geometry.findCell(x=offset, y=coords['y'][i], z=coords['z'][j])

            for index, group in enumerate(energy_groups):
                fluxes[index][j][i] = solver.getFluxByValue(cell, group, time)

    # Loop over all energy group and create a plot
    for index, group in enumerate(energy_groups):
        # Plot a 2D color map of the flat source regions
        fig = plt.figure()
        plt.imshow(np.flipud(fluxes[index, :, :]), extent=coords['bounds'])
        plt.colorbar()
        plt.title('Mesh Cell Scalar Flux in Group ' + str(group))
        filename = SUBDIRECTORY + name + '-plane-' + plane + '-group-' + str(group) + '.png'
        fig.savefig(filename, bbox_inches='tight')


def plot_precursors(solver, plane='xy', offset=0., time=3, energy_groups=[0], 
                    gridsize=250, name='precursors', xlim=None, ylim=None, zlim=None):
    
    global SUBDIRECTORY

    # Make directory if it does not exist
    if not os.path.exists(SUBDIRECTORY):
        os.makedirs(SUBDIRECTORY)

    # Error checking
    if not isinstance(solver, rk.SolverDiffusion):
        msg = 'Unable to plot the precursors since input was not ' \
              'a SolverDiffusion class object'
        raise ValueError(msg)

    if not cv.is_integer(gridsize):
        msg = 'Unable to plot the precursors since the gridsize {0} ' \
              'is not an integer'.format(gridsize)
        raise ValueError(msg)

    if plane not in ['xy', 'xz', 'yz']:
        msg = 'Unable to plot the precursors with an invalid ' \
            'plane {0}. Plane options xy, xz, yz'.format(plane)
        raise ValueError(msg)

    if gridsize <= 0:
        msg = 'Unable to plot the precursors with a negative ' \
              'gridsize {0}'.format(gridsize)
        raise ValueError(msg)

    if not isinstance(energy_groups, list):
        energy_groups = [energy_groups]

    if not cv.is_float(offset):
        msg = 'Unable to plot the precursors since the offset {0} ' \
            'is not a float'.format(offset)
        raise ValueError(msg)

    print('Plotting the precursors...')

    # Initialize a numpy array for the groupwise scalar fluxes
    precursors = numpy.zeros((len(energy_groups), gridsize, gridsize))

    geometry = solver.getGeometryDiffusion()
    coords = get_pixel_coords(geometry, plane, offset, gridsize, xlim, ylim, zlim)

    for i in range(gridsize):
        for j in range(gridsize):

            if plane == 'xy':
                cell = geometry.findCell(x=coords['x'][i], y=coords['y'][j], z=offset)
            elif plane == 'xz':
                cell = geometry.findCell(x=coords['x'][i], y=offset, z=coords['z'][j])
            else:
                cell = geometry.findCell(x=offset, y=coords['y'][i], z=coords['z'][j])

            for index, group in enumerate(energy_groups):
              if geometry.getMaterial(cell).isFissionable() is True:
                precursors[index][j][i] = geometry.getMaterial(cell)\
                                                  .getPrecursorConcByGroup(group, time)

    # Loop over all energy group and create a plot
    for index, group in enumerate(energy_groups):
        # Plot a 2D color map of the flat source regions
        fig = plt.figure()
        plt.imshow(np.flipud(precursors[index, :, :]), extent=coords['bounds'])
        plt.colorbar()
        plt.title('Mesh Cell Precursor Conc in Group ' + str(group))
        filename = SUBDIRECTORY + name + '-plane-' + plane + '-group-' + str(group) + '.png'
        fig.savefig(filename, bbox_inches='tight')

        
def plot_frequency(solver, plane='xy', offset=0., time=3, energy_groups=[0], 
                   gridsize=250, name='frequency', xlim=None, ylim=None, zlim=None):
    
    global SUBDIRECTORY

    # Make directory if it does not exist
    if not os.path.exists(SUBDIRECTORY):
        os.makedirs(SUBDIRECTORY)

    # Error checking
    if not isinstance(solver, rk.SolverDiffusion):
        msg = 'Unable to plot the frequency since input was not ' \
              'a SolverDiffusion class object'
        raise ValueError(msg)

    if not cv.is_integer(gridsize):
        msg = 'Unable to plot the frequency since the gridsize {0} ' \
              'is not an integer'.format(gridsize)
        raise ValueError(msg)

    if plane not in ['xy', 'xz', 'yz']:
        msg = 'Unable to plot the frequency with an invalid ' \
            'plane {0}. Plane options xy, xz, yz'.format(plane)
        raise ValueError(msg)

    if gridsize <= 0:
        msg = 'Unable to plot the frequency with a negative ' \
              'gridsize {0}'.format(gridsize)
        raise ValueError(msg)

    if not isinstance(energy_groups, list):
        energy_groups = [energy_groups]

    if not cv.is_float(offset):
        msg = 'Unable to plot the frequency since the offset {0} ' \
            'is not a float'.format(offset)
        raise ValueError(msg)

    print('Plotting the frequencies...')

    # Initialize a numpy array for the groupwise scalar frequencyes
    frequencies = numpy.zeros((len(energy_groups), gridsize, gridsize))

    geometry = solver.getGeometryDiffusion()
    coords = get_pixel_coords(geometry, plane, offset, gridsize, xlim, ylim, zlim)

    for i in range(gridsize):
        for j in range(gridsize):

            if plane == 'xy':
                cell = geometry.findCell(x=coords['x'][i], y=coords['y'][j], z=offset)
            elif plane == 'xz':
                cell = geometry.findCell(x=coords['x'][i], y=offset, z=coords['z'][j])
            else:
                cell = geometry.findCell(x=offset, y=coords['y'][i], z=coords['z'][j])

            for index, group in enumerate(energy_groups):
                frequencies[index][j][i] = solver.getFrequencyByValue(cell, group, time)

    # Loop over all energy group and create a plot
    for index, group in enumerate(energy_groups):
        # Plot a 2D color map of the flat source regions
        fig = plt.figure()
        plt.imshow(np.flipud(frequencies[index, :, :]), extent=coords['bounds'])
        plt.colorbar()
        plt.title('Mesh Cell Frequency in Group ' + str(group))
        filename = SUBDIRECTORY + name + '-plane-' + plane + '-group-' + str(group) + '.png'
        fig.savefig(filename, bbox_inches='tight')

        
        
def get_pixel_coords(geometry, plane, offset, gridsize, xlim, ylim, zlim):

  # initialize variables to be returned
  bounds = [geometry.getXMin(), geometry.getXMax(), geometry.getYMin(), geometry.getYMax(), geometry.getZMin(), geometry.getZMax()] 
  xcoords = None
  ycoords = None
  zcoords = None
  coords = dict()

  if xlim is not None:
    bounds[0] = xlim[0]
    bounds[1] = xlim[1]

  if ylim is not None:
    bounds[2] = ylim[0]
    bounds[3] = ylim[1]

  if zlim is not None:
    bounds[4] = zlim[0]
    bounds[5] = zlim[1]

  if plane == 'xy':
    xcoords = np.linspace(bounds[0]+TINY_MOVE, bounds[1]-TINY_MOVE, gridsize)
    ycoords = np.linspace(bounds[2]+TINY_MOVE, bounds[3]-TINY_MOVE, gridsize)
    if offset < bounds[4] or offset > bounds[5]:
      msg = 'Unable to plot offset at z={0} as it must lie ' \
          'between the z bounds [{1},{2}]'.format(offset, \
                                                   bounds[4], bounds[5])
      raise ValueError(msg)
    del bounds[4:]
  elif plane == 'xz':
    xcoords = np.linspace(bounds[0]+TINY_MOVE, bounds[1]-TINY_MOVE, gridsize)
    zcoords = np.linspace(bounds[4]+TINY_MOVE, bounds[5]-TINY_MOVE, gridsize)
    if offset < bounds[2] or offset > bounds[3]:
      msg = 'Unable to plot offset at y={0} as it must lie ' \
          'between the y bounds [{1},{2}]'.format(offset, \
                                                   bounds[2], bounds[3])
      raise ValueError(msg)
    del bounds[2:4]
  else:
    ycoords = np.linspace(bounds[2]+TINY_MOVE, bounds[3]-TINY_MOVE, gridsize)
    zcoords = np.linspace(bounds[4]+TINY_MOVE, bounds[5]-TINY_MOVE, gridsize)
    if offset < bounds[0] or offset > bounds[1]:
      msg = 'Unable to plot offset at x={0} as it must lie ' \
          'between the x bounds [{1},{2}]'.format(offset, \
                                                   bounds[0], bounds[1])
      raise ValueError(msg)
    del bounds[:2]
  
  # add attributes to coords dictionary
  coords['x'] = xcoords
  coords['y'] = ycoords
  coords['z'] = zcoords
  coords['bounds'] = bounds

  return coords
