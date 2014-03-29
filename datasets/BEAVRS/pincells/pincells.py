from surfaces import *
from materials import *


# NOTE - These are all taken at axial locations above the dashpot.
# NOTE - The guide tubes are assumed to be empty.
# NOTE - The instrument tubes are assumed to be empty.

# Keys are pin cell string names. Values are dictionaries with string name keys
# for each cell in the pin cell and OpenMOC CellBasics as values.
pincells = dict()

# Elements are all CellBasic
cells = list()


################################################################################
##########################  1.6% Enriched Fuel Pin  ############################
################################################################################

pin = {}
pin_id = universe_id()

pin['Universe'] = pin_id

pin['Fuel'] = CellBasic(pin_id, materials['1.6% Fuel'].getId())
pin['Fuel'].addSurface(halfspace=-1, surface=surfaces['Fuel Radius-1'])
cells.append(pin['Fuel'])

pin['Gap'] = CellBasic(pin_id, materials['Gap'].getId())
pin['Gap'].addSurface(halfspace=+1, surface=surfaces['Fuel Radius-1'])
pin['Gap'].addSurface(halfspace=-1, surface=surfaces['Fuel Radius-2'])
cells.append(pin['Gap'])

pin['Clad'] = CellBasic(pin_id, materials['Zircaloy'].getId())
pin['Clad'].addSurface(halfspace=+1, surface=surfaces['Fuel Radius-2'])
pin['Clad'].addSurface(halfspace=-1, surface=surfaces['Fuel Radius-3'])
cells.append(pin['Clad'])

pin['Moderator'] = CellBasic(pin_id, materials['Borated Water'].getId())
pin['Moderator'].addSurface(halfspace=+1, surface=surfaces['Fuel Radius-3'])
cells.append(pin['Moderator'])

pincells['1.6% Fuel'] = pin


################################################################################
##########################  2.4% Enriched Fuel Pin  ############################
################################################################################

pin = {}
pin_id = universe_id()

pin['Universe'] = pin_id

pin['Fuel'] = CellBasic(pin_id, materials['2.4% Fuel'].getId())
pin['Fuel'].addSurface(halfspace=-1, surface=surfaces['Fuel Radius-1'])
cells.append(pin['Fuel'])

pin['Gap'] = CellBasic(pin_id, materials['Gap'].getId())
pin['Gap'].addSurface(halfspace=+1, surface=surfaces['Fuel Radius-1'])
pin['Gap'].addSurface(halfspace=-1, surface=surfaces['Fuel Radius-2'])
cells.append(pin['Gap'])

pin['Clad'] = CellBasic(pin_id, materials['Zircaloy'].getId())
pin['Clad'].addSurface(halfspace=+1, surface=surfaces['Fuel Radius-2'])
pin['Clad'].addSurface(halfspace=-1, surface=surfaces['Fuel Radius-3'])
cells.append(pin['Clad'])

pin['Moderator'] = CellBasic(pin_id, materials['Borated Water'].getId())
pin['Moderator'].addSurface(halfspace=+1, surface=surfaces['Fuel Radius-3'])
cells.append(pin['Moderator'])

pincells['2.4% Fuel'] = pin


################################################################################
##########################  3.1% Enriched Fuel Pin  ############################
################################################################################

pin = {}
pin_id = universe_id()

pin['Universe'] = pin_id

pin['Fuel'] = CellBasic(pin_id, materials['3.1% Fuel'].getId())
pin['Fuel'].addSurface(halfspace=-1, surface=surfaces['Fuel Radius-1'])
cells.append(pin['Fuel'])

pin['Gap'] = CellBasic(pin_id, materials['Gap'].getId())
pin['Gap'].addSurface(halfspace=+1, surface=surfaces['Fuel Radius-1'])
pin['Gap'].addSurface(halfspace=-1, surface=surfaces['Fuel Radius-2'])
cells.append(pin['Gap'])

pin['Clad'] = CellBasic(pin_id, materials['Zircaloy'].getId())
pin['Clad'].addSurface(halfspace=+1, surface=surfaces['Fuel Radius-2'])
pin['Clad'].addSurface(halfspace=-1, surface=surfaces['Fuel Radius-3'])
cells.append(pin['Clad'])

pin['Moderator'] = CellBasic(pin_id, materials['Borated Water'].getId())
pin['Moderator'].addSurface(halfspace=+1, surface=surfaces['Fuel Radius-3'])
cells.append(pin['Moderator'])

pincells['3.1% Fuel'] = pin


################################################################################
################################  Guide Tube  ##################################
################################################################################

pin = {}
pin_id = universe_id()

pin['Universe'] = pin_id

pin['Inner Water'] = CellBasic(pin_id, materials['Borated Water'].getId())
pin['Inner Water'].addSurface(halfspace=-1, surface=surfaces['GT Radius-1'])
cells.append(pin['Inner Water'])

pin['Clad'] = CellBasic(pin_id, materials['Zircaloy'].getId())
pin['Clad'].addSurface(halfspace=+1, surface=surfaces['GT Radius-1'])
pin['Clad'].addSurface(halfspace=-1, surface=surfaces['GT Radius-2'])
cells.append(pin['Clad'])

pin['Outer Water'] = CellBasic(pin_id, materials['Borated Water'].getId())
pin['Outer Water'].addSurface(halfspace=+1, surface=surfaces['GT Radius-2'])
cells.append(pin['Outer Water'])

pincells['Guide Tube'] = pin


################################################################################
##############################  Instrument Tube  ###############################
################################################################################

pin = {}
pin_id = universe_id()

pin['Universe'] = pin_id

pin['Air'] = CellBasic(pin_id, materials['Air'].getId())
pin['Air'].addSurface(halfspace=-1, surface=surfaces['IT Radius-1'])
cells.append(pin['Air'])

pin['Clad'] = CellBasic(pin_id, materials['Zircaloy'].getId())
pin['Clad'].addSurface(halfspace=+1, surface=surfaces['IT Radius-1'])
pin['Clad'].addSurface(halfspace=-1, surface=surfaces['IT Radius-2'])
cells.append(pin['Clad'])

pin['Water'] = CellBasic(pin_id, materials['Borated Water'].getId())
pin['Water'].addSurface(halfspace=+1, surface=surfaces['IT Radius-2'])
cells.append(pin['Water'])

pincells['Instrument Tube'] = pin


################################################################################
#############################  Burnable Absorber  ##############################
################################################################################

pin = {}
pin_id = universe_id()

pin['Universe'] = pin_id

pin['Region-1'] = CellBasic(pin_id, materials['Air'].getId())
pin['Region-1'].addSurface(halfspace=-1, surface=surfaces['BA Radius-1'])
cells.append(pin['Region-1'])

pin['Region-2'] = CellBasic(pin_id, materials['Steel'].getId())
pin['Region-2'].addSurface(halfspace=-1, surface=surfaces['BA Radius-1'])
pin['Region-2'].addSurface(halfspace=+1, surface=surfaces['BA Radius-2'])
cells.append(pin['Region-2'])

pin['Region-3'] = CellBasic(pin_id, materials['Air'].getId())
pin['Region-3'].addSurface(halfspace=-1, surface=surfaces['BA Radius-2'])
pin['Region-3'].addSurface(halfspace=+1, surface=surfaces['BA Radius-3'])
cells.append(pin['Region-3'])

pin['Region-4'] = CellBasic(pin_id, materials['Boro. Glass'].getId())
pin['Region-4'].addSurface(halfspace=-1, surface=surfaces['BA Radius-3'])
pin['Region-4'].addSurface(halfspace=+1, surface=surfaces['BA Radius-4'])
cells.append(pin['Region-4'])

pin['Region-5'] = CellBasic(pin_id, materials['Air'].getId())
pin['Region-5'].addSurface(halfspace=-1, surface=surfaces['BA Radius-4'])
pin['Region-5'].addSurface(halfspace=+1, surface=surfaces['BA Radius-5'])
cells.append(pin['Region-5'])

pin['Region-6'] = CellBasic(pin_id, materials['Steel'].getId())
pin['Region-6'].addSurface(halfspace=-1, surface=surfaces['BA Radius-5'])
pin['Region-6'].addSurface(halfspace=+1, surface=surfaces['BA Radius-6'])
cells.append(pin['Region-6'])

pin['Region-7'] = CellBasic(pin_id, materials['Borated Water'].getId())
pin['Region-7'].addSurface(halfspace=-1, surface=surfaces['BA Radius-6'])
pin['Region-7'].addSurface(halfspace=+1, surface=surfaces['BA Radius-7'])
cells.append(pin['Region-7'])

pin['Region-8'] = CellBasic(pin_id, materials['Zircaloy'].getId())
pin['Region-8'].addSurface(halfspace=-1, surface=surfaces['BA Radius-7'])
pin['Region-8'].addSurface(halfspace=+1, surface=surfaces['BA Radius-8'])
cells.append(pin['Region-8'])

pin['Region-9'] = CellBasic(pin_id, materials['Borated Water'].getId())
pin['Region-9'].addSurface(halfspace=+1, surface=surfaces['BA Radius-8'])
cells.append(pin['Region-9'])

pincells['Burnable Absorber'] = pin