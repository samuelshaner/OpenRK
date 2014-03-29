from surfaces import *
from materials import *

# Keys are k
pincells = {}

################################################################################
##########################  1.6% Enriched Fuel Pin  ############################
################################################################################

pin = {}
pin_id = universe_id()

pin['Fuel'] = CellBasic(universe=pin_id,
                        material=materials['1.6% Fuel'].getId())
pin['Fuel'].addSurface(halfspace=-1, surface=surfaces['Fuel Radius'])

pin['Gap'] = CellBasic(universe=pin_id, material=materials['Gap'].getId())
pin['Gap'].addSurface(halfspace=+1, surface=surfaces['Fuel Radius'])
pin['Gap'].addSurface(halfspace=-1, surface=surfaces['Gap Radius'])

pin['Clad'] = CellBasic(universe=pin_id, material=materials['Zircaloy'].getId())
pin['Clad'].addSurface(halfspace=+1, surface=surfaces['Gap Radius'])
pin['Clad'].addSurface(halfspace=-1, surface=surfaces['Fuel Clad Radius'])

pin['Moderator'] = CellBasic(universe=pin_id,
                             material=materials['Borated Water'].getId())
pin['Moderator'].addSurface(halfspace=+1, surface=surfaces['Fuel Clad Radius'])

pincells['1.6% Fuel'] = pin


################################################################################
##########################  2.4% Enriched Fuel Pin  ############################
################################################################################

pin = {}
pin_id = universe_id()

pin['Fuel'] = CellBasic(universe=pin_id,
                        material=materials['2.4% Fuel'].getId())
pin['Fuel'].addSurface(halfspace=-1, surface=surfaces['Fuel Radius'])

pin['Gap'] = CellBasic(universe=pin_id, material=materials['Gap'].getId())
pin['Gap'].addSurface(halfspace=+1, surface=surfaces['Fuel Radius'])
pin['Gap'].addSurface(halfspace=-1, surface=surfaces['Gap Radius'])

pin['Clad'] = CellBasic(universe=pin_id, material=materials['Zircaloy'].getId())
pin['Clad'].addSurface(halfspace=+1, surface=surfaces['Gap Radius'])
pin['Clad'].addSurface(halfspace=-1, surface=surfaces['Fuel Clad Radius'])

pin['Moderator'] = CellBasic(universe=pin_id,
                             material=materials['Borated Water'].getId())
pin['Moderator'].addSurface(halfspace=+1, surface=surfaces['Fuel Clad Radius'])

pincells['2.4% Fuel'] = pin


################################################################################
##########################  3.1% Enriched Fuel Pin  ############################
################################################################################

pin = {}
pin_id = universe_id()

pin['Fuel'] = CellBasic(universe=pin_id,
                        material=materials['3.1% Fuel'].getId())
pin['Fuel'].addSurface(halfspace=-1, surface=surfaces['Fuel Radius'])

pin['Gap'] = CellBasic(universe=pin_id, material=materials['Gap'].getId())
pin['Gap'].addSurface(halfspace=+1, surface=surfaces['Fuel Radius'])
pin['Gap'].addSurface(halfspace=-1, surface=surfaces['Gap Radius'])

pin['Clad'] = CellBasic(universe=pin_id, material=materials['Zircaloy'].getId())
pin['Clad'].addSurface(halfspace=+1, surface=surfaces['Gap Radius'])
pin['Clad'].addSurface(halfspace=-1, surface=surfaces['Fuel Clad Radius'])

pin['Moderator'] = CellBasic(universe=pin_id,
                             material=materials['Borated Water'].getId())
pin['Moderator'].addSurface(halfspace=+1, surface=surfaces['Fuel Clad Radius'])

pincells['3.1% Fuel'] = pin


################################################################################
#############################  Empty Guide Tube  ###############################
################################################################################

# Empty guide tube
pin = {}
pin_id = universe_id()

pin['Inner Water'] = CellBasic(universe=pin_id,
                               material=materials['Borated Water'].getId())
pin['Inner Water'].addSurface(halfspace=-1, surface=surfaces['GT Radius'])

pin['Clad'] = CellBasic(universe=pin_id, material=materials['Zircaloy'].getId())
pin['Clad'].addSurface(halfspace=+1, surface=surfaces['GT Radius'])
pin['Clad'].addSurface(halfspace=-1, surface=surfaces['GT Clad Radius'])

pin['Outer Water'] = CellBasic(universe=pin_id,
                               material=materials['Borated Water'].getId())
pin['Outer Water'].addSurface(halfspace=+1, surface=surfaces['GT Clad Radius'])

pincells['Empty Guide Tube'] = pin


################################################################################
#############################  Burnable Absorber  ##############################
################################################################################

pin = {}
pin_id = universe_id()

pin['Region-1'] = CellBasic(universe=pin_id, material=materials['Air'].getId())
pin['Region-1'].addSurface(halfspace=-1, surface=surfaces['BP Radius 1'])

pin['Region-2'] = CellBasic(universe=pin_id,
                            material=materials['Steel'].getId())
pin['Region-2'].addSurface(halfspace=-1, surface=surfaces['BP Radius 2'])
pin['Region-2'].addSurface(halfspace=+1, surface=surfaces['BP Radius 1'])

pin['Region-3'] = CellBasic(universe=pin_id, material=materials['Air'].getId())
pin['Region-3'].addSurface(halfspace=-1, surface=surfaces['BP Radius 3'])
pin['Region-3'].addSurface(halfspace=+1, surface=surfaces['BP Radius 2'])

pin['Region-4'] = CellBasic(universe=pin_id,
                            material=materials['Boro. Glass'].getId())
pin['Region-4'].addSurface(halfspace=-1, surface=surfaces['BP Radius 4'])
pin['Region-4'].addSurface(halfspace=+1, surface=surfaces['BP Radius 3'])

pin['Region-5'] = CellBasic(universe=pin_id, material=materials['Air'].getId())
pin['Region-5'].addSurface(halfspace=-1, surface=surfaces['BP Radius 5'])
pin['Region-5'].addSurface(halfspace=+1, surface=surfaces['BP Radius 4'])

pin['Region-6'] = CellBasic(universe=pin_id,
                            material=materials['Steel'].getId())
pin['Region-6'].addSurface(halfspace=-1, surface=surfaces['BP Radius 6'])
pin['Region-6'].addSurface(halfspace=+1, surface=surfaces['BP Radius 5'])

pin['Region-7'] = CellBasic(universe=pin_id,
                            material=materials['Borated Water'].getId())
pin['Region-7'].addSurface(halfspace=-1, surface=surfaces['BP Radius 7'])
pin['Region-7'].addSurface(halfspace=+1, surface=surfaces['BP Radius 6'])

pin['Region-8'] = CellBasic(universe=pin_id,
                            material=materials['Zircaloy'].getId())
pin['Region-8'].addSurface(halfspace=-1, surface=surfaces['BP Radius 8'])
pin['Region-8'].addSurface(halfspace=+1, surface=surfaces['BP Radius 7'])

pin['Region-9'] = CellBasic(universe=pin_id,
                            material=materials['Borated Water'].getId())
pin['Region-9'].addSurface(halfspace=+1, surface=surfaces['BP Radius 8'])

pincells['Burnable Asborber'] = pin