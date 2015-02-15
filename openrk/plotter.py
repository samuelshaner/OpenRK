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

def plot_flux(mesh, plane='xy', offset=0., time=rk.CURRENT, energy_groups=[0], 
              gridsize=250, name='mesh-flux', xlim=None, ylim=None, zlim=None):
    
    global SUBDIRECTORY

    # Make directory if it does not exist
    if not os.path.exists(SUBDIRECTORY):
        os.makedirs(SUBDIRECTORY)

    # Error checking
    if not isinstance(mesh, rk.Mesh):
        msg = 'Unable to plot the flux since input was not ' \
              'a Mesh class object'
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

    coords = get_pixel_coords(mesh, plane, offset, gridsize, xlim, ylim, zlim)

    for i in range(gridsize):
        for j in range(gridsize):

            if plane == 'xy':
                cell = mesh.findCell(x=coords['x'][i], y=coords['y'][j], z=offset)
            elif plane == 'xz':
                cell = mesh.findCell(x=coords['x'][i], y=offset, z=coords['z'][j])
            else:
                cell = mesh.findCell(x=offset, y=coords['y'][i], z=coords['z'][j])

            for index, group in enumerate(energy_groups):
                fluxes[index][j][i] = mesh.getFluxByValue(cell, group, time)

    # Loop over all energy group and create a plot
    for index, group in enumerate(energy_groups):
        # Plot a 2D color map of the flat source regions
        fig = plt.figure()
        plt.imshow(np.flipud(fluxes[index, :, :]), extent=coords['bounds'])
        plt.colorbar()
        plt.title('Mesh Cell Scalar Flux in Group ' + str(group))
        filename = SUBDIRECTORY + name + '-plane-' + plane + '-group-' + str(group) + '.png'
        fig.savefig(filename, bbox_inches='tight')


def plot_power(mesh, plane='xy', offset=0., time=rk.CURRENT, gridsize=250, name='mesh-power',
               xlim=None, ylim=None, zlim=None):
    global SUBDIRECTORY

    # Make directory if it does not exist
    if not os.path.exists(SUBDIRECTORY):
        os.makedirs(SUBDIRECTORY)

    # Error checking
    if not isinstance(mesh, rk.Mesh):
        msg = 'Unable to plot the power since input was not ' \
              'a Mesh class object'
        raise ValueError(msg)

    if not cv.is_integer(gridsize):
        msg = 'Unable to plot the power since the gridsize {0} ' \
              'is not an integer'.format(gridsize)
        raise ValueError(msg)

    if plane not in ['xy', 'xz', 'yz']:
        msg = 'Unable to plot the power with an invalid ' \
            'plane {0}. Plane options xy, xz, yz'.format(plane)
        raise ValueError(msg)

    if gridsize <= 0:
        msg = 'Unable to plot the power with a negative ' \
              'gridsize {0}'.format(gridsize)
        raise ValueError(msg)

    if not cv.is_float(offset):
        msg = 'Unable to plot the power since the offset {0} ' \
            'is not a float'.format(offset)
        raise ValueError(msg)

    print('Plotting the power...')

    # Compute the power
    mesh.computePower(time)

    # Initialize a numpy array for the groupwise scalar fluxes
    power = numpy.zeros((gridsize, gridsize))

    coords = get_pixel_coords(mesh, plane, offset, gridsize, xlim, ylim, zlim)

    for i in range(gridsize):
        for j in range(gridsize):

            if plane == 'xy':
                cell = mesh.findCell(x=coords['x'][i], y=coords['y'][j], z=offset)
            elif plane == 'xz':
                cell = mesh.findCell(x=coords['x'][i], y=offset, z=coords['z'][j])
            else:
                cell = mesh.findCell(x=offset, y=coords['y'][i], z=coords['z'][j])

            power[j][i] = mesh.getPowerByValue(cell, time)

    # Plot a 2D color map of the flat source regions
    fig = plt.figure()
    plt.imshow(np.flipud(power[:, :]), extent=coords['bounds'])
    plt.colorbar()
    plt.title('Mesh Cell Power')
    filename = SUBDIRECTORY + name + '-plane-' + plane + '.png'
    fig.savefig(filename, bbox_inches='tight')


def plot_temperature(mesh, plane='xy', offset=0., time=rk.CURRENT, gridsize=250, 
                     name='mesh-temperature', xlim=None, ylim=None, zlim=None):

    global SUBDIRECTORY

    # Make directory if it does not exist
    if not os.path.exists(SUBDIRECTORY):
        os.makedirs(SUBDIRECTORY)

    # Error checking
    if not isinstance(mesh, rk.Mesh):
        msg = 'Unable to plot the temperature since input was not ' \
              'a Mesh class object'
        raise ValueError(msg)

    if not cv.is_integer(gridsize):
        msg = 'Unable to plot the temperature since the gridsize {0} ' \
              'is not an integer'.format(gridsize)
        raise ValueError(msg)

    if plane not in ['xy', 'xz', 'yz']:
        msg = 'Unable to plot the temperature with an invalid ' \
            'plane {0}. Plane options xy, xz, yz'.format(plane)
        raise ValueError(msg)

    if gridsize <= 0:
        msg = 'Unable to plot the temperature with a negative ' \
              'gridsize {0}'.format(gridsize)
        raise ValueError(msg)

    if not cv.is_float(offset):
        msg = 'Unable to plot the temperature since the offset {0} ' \
            'is not a float'.format(offset)
        raise ValueError(msg)

    print('Plotting the temperature...')

    # Initialize a numpy array for the groupwise scalar fluxes
    power = numpy.zeros((gridsize, gridsize))

    coords = get_pixel_coords(mesh, plane, offset, gridsize, xlim, ylim, zlim)

    for i in range(gridsize):
        for j in range(gridsize):

            if plane == 'xy':
                cell = mesh.findCell(x=coords['x'][i], y=coords['y'][j], z=offset)
            elif plane == 'xz':
                cell = mesh.findCell(x=coords['x'][i], y=offset, z=coords['z'][j])
            else:
                cell = mesh.findCell(x=offset, y=coords['y'][i], z=coords['z'][j])

            power[j][i] = mesh.getTemperatureByValue(cell, time)

    # Plot a 2D color map of the flat source regions
    fig = plt.figure()
    plt.imshow(np.flipud(power[:, :]), extent=coords['bounds'])
    plt.colorbar()
    plt.title('Mesh Cell Temperature')
    filename = SUBDIRECTORY + name + '-plane-' + plane + '.png'
    fig.savefig(filename, bbox_inches='tight')


def plot_materials(mesh, plane='xy', offset=0., gridsize=250, name='mesh-materials',
                   xlim=None, ylim=None, zlim=None):

    global SUBDIRECTORY

    # Make directory if it does not exist
    if not os.path.exists(SUBDIRECTORY):
        os.makedirs(SUBDIRECTORY)

    # Error checking
    if not isinstance(mesh, rk.Mesh):
        msg = 'Unable to plot the materials since input was not ' \
              'a Mesh class object'
        raise ValueError(msg)

    if not cv.is_integer(gridsize):
        msg = 'Unable to plot the materials since the gridsize {0} ' \
              'is not an integer'.format(gridsize)
        raise ValueError(msg)

    if plane not in ['xy', 'xz', 'yz']:
        msg = 'Unable to plot the materials with an invalid ' \
            'plane {0}. Plane options xy, xz, yz'.format(plane)
        raise ValueError(msg)

    if gridsize <= 0:
        msg = 'Unable to plot the materials with a negative ' \
              'gridsize {0}'.format(gridsize)
        raise ValueError(msg)

    if not cv.is_float(offset):
        msg = 'Unable to plot the materials since the offset {0} ' \
            'is not a float'.format(offset)
        raise ValueError(msg)

    print('Plotting the materials...')

    # Initialize a numpy array for the groupwise scalar fluxes
    materials = numpy.zeros((gridsize, gridsize))

    coords = get_pixel_coords(mesh, plane, offset, gridsize, xlim, ylim, zlim)

    unique_materials = []

    for i in range(gridsize):
        for j in range(gridsize):

            if plane == 'xy':
                cell = mesh.findCell(x=coords['x'][i], y=coords['y'][j], z=offset)
            elif plane == 'xz':
                cell = mesh.findCell(x=coords['x'][i], y=offset, z=coords['z'][j])
            else:
                cell = mesh.findCell(x=offset, y=coords['y'][i], z=coords['z'][j])

            mat_name = mesh.getMaterial(cell).getId()

            if mat_name in unique_materials:
                materials[j][i] = unique_materials.index(mat_name)
            else:
                materials[j][i] = len(unique_materials)
                unique_materials.append(mat_name)

    # Plot a 2D color map of the flat source regions
    fig = plt.figure()
    plt.imshow(np.flipud(materials[:, :]), extent=coords['bounds'])
    plt.title('Mesh Cell Materials')
    filename = SUBDIRECTORY + name + '-plane-' + plane + '.png'
    fig.savefig(filename, bbox_inches='tight')


def plot_sigma_a(mesh, plane='xy', offset=0., energy_groups=[0], gridsize=250, time=rk.CURRENT, 
                 name='mesh-sigma-a', xlim=None, ylim=None, zlim=None):

    global SUBDIRECTORY

    # Make directory if it does not exist
    if not os.path.exists(SUBDIRECTORY):
        os.makedirs(SUBDIRECTORY)

    # Error checking
    if not isinstance(mesh, rk.Mesh):
        msg = 'Unable to plot sigma a since input was not ' \
              'a Mesh class object'
        raise ValueError(msg)

    if not cv.is_integer(gridsize):
        msg = 'Unable to plot the sigma a since the gridsize {0} ' \
              'is not an integer'.format(gridsize)
        raise ValueError(msg)

    if plane not in ['xy', 'xz', 'yz']:
        msg = 'Unable to plot the sigma a with an invalid ' \
            'plane {0}. Plane options xy, xz, yz'.format(plane)
        raise ValueError(msg)

    if gridsize <= 0:
        msg = 'Unable to plot the sigma a with a negative ' \
              'gridsize {0}'.format(gridsize)
        raise ValueError(msg)

    if not cv.is_float(offset):
        msg = 'Unable to plot the sigma a since the offset {0} ' \
            'is not a float'.format(offset)
        raise ValueError(msg)

    print('Plotting the material absorption xs...')

    # Initialize a numpy array for the groupwise scalar fluxes
    sigma_a = numpy.zeros((len(energy_groups), gridsize, gridsize))

    coords = get_pixel_coords(mesh, plane, offset, gridsize, xlim, ylim, zlim)

    for i in range(gridsize):
        for j in range(gridsize):

            if plane == 'xy':
                cell = mesh.findCell(x=coords['x'][i], y=coords['y'][j], z=offset)
            elif plane == 'xz':
                cell = mesh.findCell(x=coords['x'][i], y=offset, z=coords['z'][j])
            else:
                cell = mesh.findCell(x=offset, y=coords['y'][i], z=coords['z'][j])

            mat = mesh.getMaterial(cell)
            temp = mesh.getTemperatureByValue(cell, time)

            for index, group in enumerate(energy_groups):
                sigma_a[index][j][i] = mat.getSigmaAByGroup(group, time, temp)

    # Loop over all energy group and create a plot
    for index, group in enumerate(energy_groups):
        # Plot a 2D color map of the flat source regions
        fig = plt.figure()
        plt.imshow(np.flipud(sigma_a[index, :, :]), extent=coords['bounds'])
        plt.colorbar()
        plt.title('Mesh Cell Sigma A in Group ' + str(group))
        filename = SUBDIRECTORY + name + '-plane-' + plane + '-group-' + str(group) + '.png'
        fig.savefig(filename, bbox_inches='tight')


def plot_precursor_conc(mesh, plane='xy', offset=0., time=rk.CURRENT, delayed_groups=[0], 
                        gridsize=250, name='mesh-precursor-conc', xlim=None, ylim=None, zlim=None):

    global SUBDIRECTORY

    # Make directory if it does not exist
    if not os.path.exists(SUBDIRECTORY):
        os.makedirs(SUBDIRECTORY)

    # Error checking
    if not isinstance(mesh, rk.Mesh):
        msg = 'Unable to plot the delayed precursors since input was not ' \
              'a Mesh class object'
        raise ValueError(msg)

    if not cv.is_integer(gridsize):
        msg = 'Unable to plot the delayed precursors since the gridsize {0} ' \
              'is not an integer'.format(gridsize)
        raise ValueError(msg)

    if plane not in ['xy', 'xz', 'yz']:
        msg = 'Unable to plot the delayed precursors with an invalid ' \
            'plane {0}. Plane options xy, xz, yz'.format(plane)
        raise ValueError(msg)

    if gridsize <= 0:
        msg = 'Unable to plot the delayed precursors with a negative ' \
              'gridsize {0}'.format(gridsize)
        raise ValueError(msg)

    if not cv.is_float(offset):
        msg = 'Unable to plot the delayed precursors since the offset {0} ' \
            'is not a float'.format(offset)
        raise ValueError(msg)

    if not isinstance(delayed_groups, list):
        delayed_groups = [delayed_groups]

    print('Plotting the precursor conc...')

    # Initialize a numpy array for the groupwise precursor conc
    precursor_conc = numpy.zeros((len(delayed_groups), gridsize, gridsize))

    coords = get_pixel_coords(mesh, plane, offset, gridsize, xlim, ylim, zlim)

    for i in range(gridsize):
        for j in range(gridsize):

            if plane == 'xy':
                cell = mesh.findCell(x=coords['x'][i], y=coords['y'][j], z=offset)
            elif plane == 'xz':
                cell = mesh.findCell(x=coords['x'][i], y=offset, z=coords['z'][j])
            else:
                cell = mesh.findCell(x=offset, y=coords['y'][i], z=coords['z'][j])

            for index, group in enumerate(delayed_groups):
                precursor_conc[index][j][i] = mesh.getMaterial(cell).getPrecursorConcByGroup(group, time)

    # Loop over all energy group and create a plot
    for index, group in enumerate(delayed_groups):
        # Plot a 2D color map of the flat source regions
        fig = plt.figure()
        plt.imshow(np.flipud(precursor_conc[index, :, :]), extent=coords['bounds'])
        plt.colorbar()
        plt.title('Mesh Cell Precursor Conc in Group ' + str(group))
        filename = SUBDIRECTORY + name + '-plane-' + plane + '-group-' + str(group) + '.png'
        fig.savefig(filename, bbox_inches='tight')


def get_pixel_coords(mesh, plane, offset, gridsize, xlim, ylim, zlim):

  # initialize variables to be returned
  bounds = [mesh.getXMin(), mesh.getXMax(), mesh.getYMin(), mesh.getYMax(), mesh.getZMin(), mesh.getZMax()] 
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
