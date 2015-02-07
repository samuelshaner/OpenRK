__author__ = 'Sam Shaner'
__email__ = 'shaner@mit.edu'

# Import modules
import matplotlib
from material import Material, FunctionalMaterial
from mesh import *
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


def plot_flux(mesh, time='CURRENT', energy_groups=[0], gridsize=250, name='mesh-flux'):
    global SUBDIRECTORY

    # Make directory if it does not exist
    if not os.path.exists(SUBDIRECTORY):
        os.makedirs(SUBDIRECTORY)

    # Error checking
    if not isinstance(mesh, Mesh):
        msg = 'Unable to plot the cells since input was not ' \
              'a Mesh class object'
        raise ValueError(msg)

    if not cv.is_integer(gridsize):
        msg = 'Unable to plot the cells since the gridsize {0} ' \
              'is not an integer'.format(gridsize)
        raise ValueError(msg)

    if gridsize <= 0:
        msg = 'Unable to plot the cells with a negative ' \
              'gridsize {0}'.format(gridsize)
        raise ValueError(msg)

    if not isinstance(energy_groups, list):
        energy_groups = [energy_groups]

    print('Plotting the fluxes...')

    # Initialize a numpy array for the groupwise scalar fluxes
    fluxes = numpy.zeros((len(energy_groups), gridsize, gridsize))

    tiny_move = 1.0e-8
    bounds = mesh.get_bounds()

    # Retrieve the bounding box for the geometry
    xmin = bounds[0] + tiny_move
    xmax = bounds[1] - tiny_move
    ymin = bounds[2] + tiny_move
    ymax = bounds[3] - tiny_move

    # Initialize numpy arrays for the grid points
    xcoords = np.linspace(xmin, xmax, gridsize)
    ycoords = np.linspace(ymin, ymax, gridsize)

    for i in range(gridsize):
        for j in range(gridsize):

            # Find the flat source region IDs for each grid point
            x = xcoords[i]
            y = ycoords[j]

            cell = mesh.find_cell(x, y)

            for index, group in enumerate(energy_groups):
                fluxes[index][j][i] = mesh.get_flux_by_value(cell, group, time)

    # Loop over all energy group and create a plot
    for index, group in enumerate(energy_groups):
        # Plot a 2D color map of the flat source regions
        fig = plt.figure()
        plt.imshow(np.flipud(fluxes[index, :, :]), extent=[xmin, xmax, ymin, ymax])
        plt.colorbar()
        plt.title('Mesh Cell Scalar Flux in Group ' + str(group))
        filename = SUBDIRECTORY + name + '-group-' + str(group) + '.png'
        fig.savefig(filename, bbox_inches='tight')


def plot_power(mesh, time='CURRENT', gridsize=250, name='mesh-power'):
    global SUBDIRECTORY

    # Make directory if it does not exist
    if not os.path.exists(SUBDIRECTORY):
        os.makedirs(SUBDIRECTORY)

    # Error checking
    if not isinstance(mesh, Mesh):
        msg = 'Unable to plot the cells since input was not ' \
              'a Mesh class object'
        raise ValueError(msg)

    if not cv.is_integer(gridsize):
        msg = 'Unable to plot the cells since the gridsize {0} ' \
              'is not an integer'.format(gridsize)
        raise ValueError(msg)

    if gridsize <= 0:
        msg = 'Unable to plot the cells with a negative ' \
              'gridsize {0}'.format(gridsize)
        raise ValueError(msg)

    print('Plotting the power...')

    # Compute the power
    mesh.compute_power(time)

    # Initialize a numpy array for the groupwise scalar fluxes
    power = numpy.zeros((gridsize, gridsize))

    tiny_move = 1.0e-8
    bounds = mesh.get_bounds()

    # Retrieve the bounding box for the geometry
    xmin = bounds[0] + tiny_move
    xmax = bounds[1] - tiny_move
    ymin = bounds[2] + tiny_move
    ymax = bounds[3] - tiny_move

    # Initialize numpy arrays for the grid points
    xcoords = np.linspace(xmin, xmax, gridsize)
    ycoords = np.linspace(ymin, ymax, gridsize)

    for i in range(gridsize):
        for j in range(gridsize):

            # Find the flat source region IDs for each grid point
            x = xcoords[i]
            y = ycoords[j]

            cell = mesh.find_cell(x, y)

            power[j][i] = mesh.get_power_by_value(cell, time)

    # Plot a 2D color map of the flat source regions
    fig = plt.figure()
    plt.imshow(np.flipud(power[:, :]), extent=[xmin, xmax, ymin, ymax])
    plt.colorbar()
    plt.title('Mesh Cell Power')
    filename = SUBDIRECTORY + name + '.png'
    fig.savefig(filename, bbox_inches='tight')


def plot_temperature(mesh, time='CURRENT', gridsize=250, name='mesh-temperature'):
    global SUBDIRECTORY

    # Make directory if it does not exist
    if not os.path.exists(SUBDIRECTORY):
        os.makedirs(SUBDIRECTORY)

    # Error checking
    if not isinstance(mesh, Mesh):
        msg = 'Unable to plot the cells since input was not ' \
              'a Mesh class object'
        raise ValueError(msg)

    if not cv.is_integer(gridsize):
        msg = 'Unable to plot the cells since the gridsize {0} ' \
              'is not an integer'.format(gridsize)
        raise ValueError(msg)

    if gridsize <= 0:
        msg = 'Unable to plot the cells with a negative ' \
              'gridsize {0}'.format(gridsize)
        raise ValueError(msg)

    print('Plotting the temperature...')

    # Initialize a numpy array for the groupwise scalar fluxes
    power = numpy.zeros((gridsize, gridsize))

    tiny_move = 1.0e-8
    bounds = mesh.get_bounds()

    # Retrieve the bounding box for the geometry
    xmin = bounds[0] + tiny_move
    xmax = bounds[1] - tiny_move
    ymin = bounds[2] + tiny_move
    ymax = bounds[3] - tiny_move

    # Initialize numpy arrays for the grid points
    xcoords = np.linspace(xmin, xmax, gridsize)
    ycoords = np.linspace(ymin, ymax, gridsize)

    for i in range(gridsize):
        for j in range(gridsize):

            # Find the flat source region IDs for each grid point
            x = xcoords[i]
            y = ycoords[j]

            cell = mesh.find_cell(x, y)

            power[j][i] = mesh.get_temperature_by_value(cell, time)

    # Plot a 2D color map of the flat source regions
    fig = plt.figure()
    plt.imshow(np.flipud(power[:, :]), extent=[xmin, xmax, ymin, ymax])
    plt.colorbar()
    plt.title('Mesh Cell Temperature')
    filename = SUBDIRECTORY + name + '.png'
    fig.savefig(filename, bbox_inches='tight')


def plot_materials(mesh, gridsize=250, name='mesh-materials'):
    global SUBDIRECTORY

    # Make directory if it does not exist
    if not os.path.exists(SUBDIRECTORY):
        os.makedirs(SUBDIRECTORY)

    # Error checking
    if not isinstance(mesh, Mesh):
        msg = 'Unable to plot the cells since input was not ' \
              'a Mesh class object'
        raise ValueError(msg)

    if not cv.is_integer(gridsize):
        msg = 'Unable to plot the cells since the gridsize {0} ' \
              'is not an integer'.format(gridsize)
        raise ValueError(msg)

    if gridsize <= 0:
        msg = 'Unable to plot the cells with a negative ' \
              'gridsize {0}'.format(gridsize)
        raise ValueError(msg)

    print('Plotting the materials...')

    # Initialize a numpy array for the groupwise scalar fluxes
    materials = numpy.zeros((gridsize, gridsize))

    tiny_move = 1.0e-8
    bounds = mesh.get_bounds()

    # Retrieve the bounding box for the geometry
    xmin = bounds[0] + tiny_move
    xmax = bounds[1] - tiny_move
    ymin = bounds[2] + tiny_move
    ymax = bounds[3] - tiny_move

    # Initialize numpy arrays for the grid points
    xcoords = np.linspace(xmin, xmax, gridsize)
    ycoords = np.linspace(ymin, ymax, gridsize)

    unique_materials = []

    for i in range(gridsize):
        for j in range(gridsize):

            # Find the flat source region IDs for each grid point
            x = xcoords[i]
            y = ycoords[j]

            mat_name = mesh.get_material(mesh.find_cell(x, y)).get_name()

            if mat_name in unique_materials:
                materials[j][i] = unique_materials.index(mat_name)
            else:
                materials[j][i] = len(unique_materials)
                unique_materials.append(mat_name)

    # Plot a 2D color map of the flat source regions
    fig = plt.figure()
    plt.imshow(np.flipud(materials[:, :]), extent=[xmin, xmax, ymin, ymax])
    plt.title('Mesh Cell Materials')
    filename = SUBDIRECTORY + name + '.png'
    fig.savefig(filename, bbox_inches='tight')


def plot_sigma_a(mesh, energy_groups=[0], gridsize=250, time='CURRENT', name='mesh-sigma-a'):
    global SUBDIRECTORY

    # Make directory if it does not exist
    if not os.path.exists(SUBDIRECTORY):
        os.makedirs(SUBDIRECTORY)

    # Error checking
    if not isinstance(mesh, Mesh):
        msg = 'Unable to plot the cells since input was not ' \
              'a Mesh class object'
        raise ValueError(msg)

    if not cv.is_integer(gridsize):
        msg = 'Unable to plot the cells since the gridsize {0} ' \
              'is not an integer'.format(gridsize)
        raise ValueError(msg)

    if gridsize <= 0:
        msg = 'Unable to plot the cells with a negative ' \
              'gridsize {0}'.format(gridsize)
        raise ValueError(msg)

    print('Plotting the material absorption xs...')

    # Initialize a numpy array for the groupwise scalar fluxes
    sigma_a = numpy.zeros((len(energy_groups), gridsize, gridsize))

    tiny_move = 1.0e-8
    bounds = mesh.get_bounds()

    # Retrieve the bounding box for the geometry
    xmin = bounds[0] + tiny_move
    xmax = bounds[1] - tiny_move
    ymin = bounds[2] + tiny_move
    ymax = bounds[3] - tiny_move

    # Initialize numpy arrays for the grid points
    xcoords = np.linspace(xmin, xmax, gridsize)
    ycoords = np.linspace(ymin, ymax, gridsize)

    for i in range(gridsize):
        for j in range(gridsize):

            # Find the flat source region IDs for each grid point
            x = xcoords[i]
            y = ycoords[j]

            cell = mesh.find_cell(x, y)
            mat = mesh.get_material(cell)
            temp = mesh.get_temperature(time)[cell]

            for index, group in enumerate(energy_groups):
                sigma_a[index][j][i] = mat.get_sigma_a_by_group(group, time, temp)

    # Loop over all energy group and create a plot
    for index, group in enumerate(energy_groups):
        # Plot a 2D color map of the flat source regions
        fig = plt.figure()
        plt.imshow(np.flipud(sigma_a[index, :, :]), extent=[xmin, xmax, ymin, ymax])
        plt.colorbar()
        plt.title('Mesh Cell Sigma A in Group ' + str(group))
        filename = SUBDIRECTORY + name + '-group-' + str(group) + '.png'
        fig.savefig(filename, bbox_inches='tight')


def plot_precursor_conc(mesh, time='CURRENT', delayed_groups=[0], gridsize=250, name='mesh-precursor-conc'):
    global SUBDIRECTORY

    # Make directory if it does not exist
    if not os.path.exists(SUBDIRECTORY):
        os.makedirs(SUBDIRECTORY)

    # Error checking
    if not isinstance(mesh, Mesh):
        msg = 'Unable to plot the cells since input was not ' \
              'a Mesh class object'
        raise ValueError(msg)

    if not cv.is_integer(gridsize):
        msg = 'Unable to plot the cells since the gridsize {0} ' \
              'is not an integer'.format(gridsize)
        raise ValueError(msg)

    if gridsize <= 0:
        msg = 'Unable to plot the cells with a negative ' \
              'gridsize {0}'.format(gridsize)
        raise ValueError(msg)

    if not isinstance(delayed_groups, list):
        delayed_groups = [delayed_groups]

    print('Plotting the precursor conc...')

    # Initialize a numpy array for the groupwise precursor conc
    precursor_conc = numpy.zeros((len(delayed_groups), gridsize, gridsize))

    tiny_move = 1.0e-8
    bounds = mesh.get_bounds()

    # Retrieve the bounding box for the geometry
    xmin = bounds[0] + tiny_move
    xmax = bounds[1] - tiny_move
    ymin = bounds[2] + tiny_move
    ymax = bounds[3] - tiny_move

    # Initialize numpy arrays for the grid points
    xcoords = np.linspace(xmin, xmax, gridsize)
    ycoords = np.linspace(ymin, ymax, gridsize)

    for i in range(gridsize):
        for j in range(gridsize):

            # Find the flat source region IDs for each grid point
            x = xcoords[i]
            y = ycoords[j]

            cell = mesh.find_cell(x, y)

            for index, group in enumerate(delayed_groups):
                if isinstance(mesh.get_material(cell), TransientMaterial):
                    precursor_conc[index][j][i] = mesh.get_material(cell).get_precursor_conc_by_group(group, time)

    # Loop over all energy group and create a plot
    for index, group in enumerate(delayed_groups):
        # Plot a 2D color map of the flat source regions
        fig = plt.figure()
        plt.imshow(np.flipud(precursor_conc[index, :, :]), extent=[xmin, xmax, ymin, ymax])
        plt.colorbar()
        plt.title('Mesh Cell Precursor Conc in Group ' + str(group))
        filename = SUBDIRECTORY + name + '-group-' + str(group) + '.png'
        fig.savefig(filename, bbox_inches='tight')