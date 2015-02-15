from opencg import Universe, Cell
from datasets.BEAVRS.surfaces import *
from datasets.BEAVRS.materials import opencg_materials

# NOTE - These are all taken at axial locations above the dashpot.
# NOTE - The guide tubes are assumed to be empty.
# NOTE - The instrument tubes are assumed to be empty.

# Keys are pin cell string names. Values are dictionaries with string name keys
# for each cell in the pin cell and OpenCG Cells as values.
pin_cells = dict()


################################################################################
##########################  1.6% Enriched Fuel Pin  ############################
################################################################################


universe = Universe(name='1.6% Fuel Pin')
pin_cells[universe._name] = universe

cell = Cell(name='Fuel', fill=opencg_materials['1.6% Fuel'])
cell.addSurface(surface=surfaces['Fuel Radius-1'], halfspace=-1)
universe.addCell(cell)

cell = Cell(name='Gap', fill=opencg_materials['Helium'])
cell.addSurface(surface=surfaces['Fuel Radius-1'], halfspace=+1)
cell.addSurface(surface=surfaces['Fuel Radius-2'], halfspace=-1)
universe.addCell(cell)

cell = Cell(name='Clad', fill=opencg_materials['Zircaloy'])
cell.addSurface(surface=surfaces['Fuel Radius-2'], halfspace=+1)
cell.addSurface(surface=surfaces['Fuel Radius-3'], halfspace=-1)
universe.addCell(cell)

cell = Cell(name='Moderator', fill=opencg_materials['Borated Water'])
cell.addSurface(surface=surfaces['Fuel Radius-3'], halfspace=+1)
universe.addCell(cell)


################################################################################
##########################  2.4% Enriched Fuel Pin  ############################
################################################################################

universe = Universe(name='2.4% Fuel Pin')
pin_cells[universe._name] = universe

cell = Cell(name='Fuel', fill=opencg_materials['2.4% Fuel'])
cell.addSurface(surface=surfaces['Fuel Radius-1'], halfspace=-1)
universe.addCell(cell)

cell = Cell(name='Gap', fill=opencg_materials['Helium'])
cell.addSurface(surface=surfaces['Fuel Radius-1'], halfspace=+1)
cell.addSurface(surface=surfaces['Fuel Radius-2'], halfspace=-1)
universe.addCell(cell)

cell = Cell(name='Clad', fill=opencg_materials['Zircaloy'])
cell.addSurface(surface=surfaces['Fuel Radius-2'], halfspace=+1)
cell.addSurface(surface=surfaces['Fuel Radius-3'], halfspace=-1)
universe.addCell(cell)

cell = Cell(name='Moderator', fill=opencg_materials['Borated Water'])
cell.addSurface(surface=surfaces['Fuel Radius-3'], halfspace=+1)
universe.addCell(cell)


################################################################################
##########################  3.1% Enriched Fuel Pin  ############################
################################################################################

universe = Universe(name='3.1% Fuel Pin')
pin_cells[universe._name] = universe

cell = Cell(name='Fuel', fill=opencg_materials['3.1% Fuel'])
cell.addSurface(surface=surfaces['Fuel Radius-1'], halfspace=-1)
universe.addCell(cell)

cell = Cell(name='Gap', fill=opencg_materials['Helium'])
cell.addSurface(surface=surfaces['Fuel Radius-1'], halfspace=+1)
cell.addSurface(surface=surfaces['Fuel Radius-2'], halfspace=-1)
universe.addCell(cell)

cell = Cell(name='Clad', fill=opencg_materials['Zircaloy'])
cell.addSurface(surface=surfaces['Fuel Radius-2'], halfspace=+1)
cell.addSurface(surface=surfaces['Fuel Radius-3'], halfspace=-1)
universe.addCell(cell)

cell = Cell(name='Moderator', fill=opencg_materials['Borated Water'])
cell.addSurface(surface=surfaces['Fuel Radius-3'], halfspace=+1)
universe.addCell(cell)


################################################################################
################################  Guide Tube  ##################################
################################################################################

universe = Universe(name='Guide Tube')
pin_cells[universe._name] = universe

cell = Cell(name='Inner Water', fill=opencg_materials['Borated Water'])
cell.addSurface(surface=surfaces['GT Radius-1'], halfspace=-1)
universe.addCell(cell)

cell = Cell(name='Clad', fill=opencg_materials['Zircaloy'])
cell.addSurface(surface=surfaces['GT Radius-1'], halfspace=+1)
cell.addSurface(surface=surfaces['GT Radius-2'], halfspace=-1)
universe.addCell(cell)

cell = Cell(name='Outer Water', fill=opencg_materials['Borated Water'])
cell.addSurface(surface=surfaces['GT Radius-2'], halfspace=+1)
universe.addCell(cell)


################################################################################
##############################  Instrument Tube  ###############################
################################################################################

universe = Universe(name='Instrument Tube')
pin_cells[universe._name] = universe

cell = Cell(name='Air', fill=opencg_materials['Air'])
cell.addSurface(surface=surfaces['IT Radius-1'], halfspace=-1)
universe.addCell(cell)

cell = Cell(name='Clad', fill=opencg_materials['Zircaloy'])
cell.addSurface(surface=surfaces['IT Radius-1'], halfspace=+1)
cell.addSurface(surface=surfaces['IT Radius-2'], halfspace=-1)
universe.addCell(cell)

cell = Cell(name='Moderator', fill=opencg_materials['Borated Water'])
cell.addSurface(surface=surfaces['IT Radius-2'], halfspace=+1)
universe.addCell(cell)


################################################################################
#############################  Burnable Absorber  ##############################
################################################################################

universe = Universe(name='Burnable Absorber')
pin_cells[universe._name] = universe

cell = Cell(name='Region-1', fill=opencg_materials['Air'])
cell.addSurface(surface=surfaces['BA Radius-1'], halfspace=-1)
universe.addCell(cell)

cell = Cell(name='Region-2', fill=opencg_materials['Steel'])
cell.addSurface(surface=surfaces['BA Radius-1'], halfspace=+1)
cell.addSurface(surface=surfaces['BA Radius-2'], halfspace=-1)
universe.addCell(cell)

cell = Cell(name='Region-3', fill=opencg_materials['Air'])
cell.addSurface(surface=surfaces['BA Radius-2'], halfspace=+1)
cell.addSurface(surface=surfaces['BA Radius-3'], halfspace=-1)
universe.addCell(cell)

cell = Cell(name='Region-4', fill=opencg_materials['Boro. Glass'])
cell.addSurface(surface=surfaces['BA Radius-3'], halfspace=+1)
cell.addSurface(surface=surfaces['BA Radius-4'], halfspace=-1)
universe.addCell(cell)

cell = Cell(name='Region-5', fill=opencg_materials['Air'])
cell.addSurface(surface=surfaces['BA Radius-4'], halfspace=+1)
cell.addSurface(surface=surfaces['BA Radius-5'], halfspace=-1)
universe.addCell(cell)

cell = Cell(name='Region-6', fill=opencg_materials['Steel'])
cell.addSurface(surface=surfaces['BA Radius-5'], halfspace=+1)
cell.addSurface(surface=surfaces['BA Radius-6'], halfspace=-1)
universe.addCell(cell)

cell = Cell(name='Region-7', fill=opencg_materials['Borated Water'])
cell.addSurface(surface=surfaces['BA Radius-6'], halfspace=+1)
cell.addSurface(surface=surfaces['BA Radius-7'], halfspace=-1)
universe.addCell(cell)

cell = Cell(name='Region-8', fill=opencg_materials['Zircaloy'])
cell.addSurface(surface=surfaces['BA Radius-7'], halfspace=+1)
cell.addSurface(surface=surfaces['BA Radius-8'], halfspace=-1)
universe.addCell(cell)

cell = Cell(name='Region-9', fill=opencg_materials['Borated Water'])
cell.addSurface(surface=surfaces['BA Radius-8'], halfspace=+1)
universe.addCell(cell)