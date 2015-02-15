from openmoc.compatible.opencg_compatible import *
from geometry import geometry


###############################################################################
##################   Exporting to OpenMC geometry.xml File  ###################
###############################################################################

openmoc_geometry = get_openmoc_geometry(geometry)
openmoc_geometry.initializeFlatSourceRegions()

import openmoc.plotter as plotter
plotter.plot_cells(openmoc_geometry)