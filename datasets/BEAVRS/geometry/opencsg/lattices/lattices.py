import copy

from pincells import *
from datasets.BEAVRS.opencsg.lattices.templates import *


# Keys are Lattice string names and values are OpenMOC Lattices.
lattices = dict()


###############################################################################
####################   BEAVRS Fuel Assembly Parameters   ######################
###############################################################################

pin_pitch = 1.25984                                        # centimeters
lattice_pin_dims = 17.                                     # 17 pins
lattice_width = lattice_pin_dims * pin_pitch               # Lattice width [cm]


###############################################################################
####################   1.6% Enriched Fuel Assembly   ##########################
###############################################################################

cells = copy.deepcopy(templates['0BA'])

for i in range(len(cells)):

  row = cells[i]

  # Fuel Pins
  univ_id = pincells['1.6% Fuel']['Universe']
  row = [univ_id if x == masks['Fuel Pin'] else x for x in row]

  # Guide Tubes
  univ_id = pincells['Guide Tube']['Universe']
  row = [univ_id if x == masks['Guide Tube'] else x for x in row]

  # Instrument Tubes
  univ_id = pincells['Instrument Tube']['Universe']
  row = [univ_id if x == masks['Instrument Tube'] else x for x in row]

  cells[i] = row

lattice = Lattice(id=universe_id(), width_x=pin_pitch, width_y=pin_pitch)
lattice.setLatticeCells(cells)

lattices['1.6% Fuel - 0BA'] = lattice


###############################################################################
####################   2.4% Enriched Fuel Assembly   ##########################
###############################################################################

cells = copy.deepcopy(templates['0BA'])

for i in range(len(cells)):

  row = cells[i]

  # Fuel Pins
  univ_id = pincells['2.4% Fuel']['Universe']
  row = [univ_id if x == masks['Fuel Pin'] else x for x in row]

  # Guide Tubes
  univ_id = pincells['Guide Tube']['Universe']
  row = [univ_id if x == masks['Guide Tube'] else x for x in row]

  # Instrument Tubes
  univ_id = pincells['Instrument Tube']['Universe']
  row = [univ_id if x == masks['Instrument Tube'] else x for x in row]

  cells[i] = row

lattice = Lattice(id=universe_id(), width_x=pin_pitch, width_y=pin_pitch)
lattice.setLatticeCells(cells)

lattices['2.4% Fuel - 0BA'] = lattice


###############################################################################
################   2.4% Enriched Fuel Assembly w/ 12 BAs  #####################
###############################################################################

cells = copy.deepcopy(templates['12BA'])

for i in range(len(cells)):

  row = cells[i]

  # Fuel Pins
  univ_id = pincells['2.4% Fuel']['Universe']
  row = [univ_id if x == masks['Fuel Pin'] else x for x in row]

  # Guide Tubes
  univ_id = pincells['Guide Tube']['Universe']
  row = [univ_id if x == masks['Guide Tube'] else x for x in row]

  # Instrument Tubes
  univ_id = pincells['Instrument Tube']['Universe']
  row = [univ_id if x == masks['Instrument Tube'] else x for x in row]

  # Burnable Absorbers
  univ_id = pincells['Burnable Absorber']['Universe']
  row = [univ_id if x == masks['Burnable Absorber'] else x for x in row]

  cells[i] = row

lattice = Lattice(id=universe_id(), width_x=pin_pitch, width_y=pin_pitch)
lattice.setLatticeCells(cells)

lattices['2.4% Fuel - 12BA'] = lattice


###############################################################################
################   2.4% Enriched Fuel Assembly w/ 16 BAs  #####################
###############################################################################

cells = copy.deepcopy(templates['16BA'])

for i in range(len(cells)):

  row = cells[i]

  # Fuel Pins
  univ_id = pincells['2.4% Fuel']['Universe']
  row = [univ_id if x == masks['Fuel Pin'] else x for x in row]

  # Guide Tubes
  univ_id = pincells['Guide Tube']['Universe']
  row = [univ_id if x == masks['Guide Tube'] else x for x in row]

  # Instrument Tubes
  univ_id = pincells['Instrument Tube']['Universe']
  row = [univ_id if x == masks['Instrument Tube'] else x for x in row]

  # Burnable Absorbers
  univ_id = pincells['Burnable Absorber']['Universe']
  row = [univ_id if x == masks['Burnable Absorber'] else x for x in row]

  cells[i] = row

lattice = Lattice(id=universe_id(), width_x=pin_pitch, width_y=pin_pitch)
lattice.setLatticeCells(cells)

lattices['2.4% Fuel - 16BA'] = lattice


###############################################################################
#####################   3.1% Enriched Fuel Assembly  ##########################
###############################################################################

cells = copy.deepcopy(templates['0BA'])

for i in range(len(cells)):

  row = cells[i]

  # Fuel Pins
  univ_id = pincells['3.1% Fuel']['Universe']
  row = [univ_id if x == masks['Fuel Pin'] else x for x in row]

  # Guide Tubes
  univ_id = pincells['Guide Tube']['Universe']
  row = [univ_id if x == masks['Guide Tube'] else x for x in row]

  # Instrument Tubes
  univ_id = pincells['Instrument Tube']['Universe']
  row = [univ_id if x == masks['Instrument Tube'] else x for x in row]

  cells[i] = row

lattice = Lattice(id=universe_id(), width_x=pin_pitch, width_y=pin_pitch)
lattice.setLatticeCells(cells)

lattices['3.1% Fuel - 0BA'] = lattice


###############################################################################
################   3.1% Enriched Fuel Assembly w/ 6 BAs  ######################
###############################################################################

cells = copy.deepcopy(templates['6BA'])

for i in range(len(cells)):

  row = cells[i]

  # Fuel Pins
  univ_id = pincells['3.1% Fuel']['Universe']
  row = [univ_id if x == masks['Fuel Pin'] else x for x in row]

  # Guide Tubes
  univ_id = pincells['Guide Tube']['Universe']
  row = [univ_id if x == masks['Guide Tube'] else x for x in row]

  # Instrument Tubes
  univ_id = pincells['Instrument Tube']['Universe']
  row = [univ_id if x == masks['Instrument Tube'] else x for x in row]

  # Burnable Absorbers
  univ_id = pincells['Burnable Absorber']['Universe']
  row = [univ_id if x == masks['Burnable Absorber'] else x for x in row]

  cells[i] = row

lattice = Lattice(id=universe_id(), width_x=pin_pitch, width_y=pin_pitch)
lattice.setLatticeCells(cells)

lattices['3.1% Fuel - 6BA'] = lattice


###############################################################################
################   3.1% Enriched Fuel Assembly w/ 15 BAs  #####################
###############################################################################

cells = copy.deepcopy(templates['15BA'])

for i in range(len(cells)):

  row = cells[i]

  # Fuel Pins
  univ_id = pincells['3.1% Fuel']['Universe']
  row = [univ_id if x == masks['Fuel Pin'] else x for x in row]

  # Guide Tubes
  univ_id = pincells['Guide Tube']['Universe']
  row = [univ_id if x == masks['Guide Tube'] else x for x in row]

  # Instrument Tubes
  univ_id = pincells['Instrument Tube']['Universe']
  row = [univ_id if x == masks['Instrument Tube'] else x for x in row]

  # Burnable Absorbers
  univ_id = pincells['Burnable Absorber']['Universe']
  row = [univ_id if x == masks['Burnable Absorber'] else x for x in row]

  cells[i] = row

lattice = Lattice(id=universe_id(), width_x=pin_pitch, width_y=pin_pitch)
lattice.setLatticeCells(cells)

lattices['3.1% Fuel - 15BA'] = lattice


###############################################################################
################   3.1% Enriched Fuel Assembly w/ 16 BAs  #####################
###############################################################################

cells = copy.deepcopy(templates['16BA'])

for i in range(len(cells)):

  row = cells[i]

  # Fuel Pins
  univ_id = pincells['3.1% Fuel']['Universe']
  row = [univ_id if x == masks['Fuel Pin'] else x for x in row]

  # Guide Tubes
  univ_id = pincells['Guide Tube']['Universe']
  row = [univ_id if x == masks['Guide Tube'] else x for x in row]

  # Instrument Tubes
  univ_id = pincells['Instrument Tube']['Universe']
  row = [univ_id if x == masks['Instrument Tube'] else x for x in row]

  # Burnable Absorbers
  univ_id = pincells['Burnable Absorber']['Universe']
  row = [univ_id if x == masks['Burnable Absorber'] else x for x in row]

  cells[i] = row

lattice = Lattice(id=universe_id(), width_x=pin_pitch, width_y=pin_pitch)
lattice.setLatticeCells(cells)

lattices['3.1% Fuel - 16BA'] = lattice


###############################################################################
################   3.1% Enriched Fuel Assembly w/ 20 BAs  #####################
###############################################################################

cells = copy.deepcopy(templates['20BA'])

for i in range(len(cells)):

  row = cells[i]

  # Fuel Pins
  univ_id = pincells['3.1% Fuel']['Universe']
  row = [univ_id if x == masks['Fuel Pin'] else x for x in row]

  # Guide Tubes
  univ_id = pincells['Guide Tube']['Universe']
  row = [univ_id if x == masks['Guide Tube'] else x for x in row]

  # Instrument Tubes
  univ_id = pincells['Instrument Tube']['Universe']
  row = [univ_id if x == masks['Instrument Tube'] else x for x in row]

  # Burnable Absorbers
  univ_id = pincells['Burnable Absorber']['Universe']
  row = [univ_id if x == masks['Burnable Absorber'] else x for x in row]

  cells[i] = row

lattice = Lattice(id=universe_id(), width_x=pin_pitch, width_y=pin_pitch)
lattice.setLatticeCells(cells)

lattices['3.1% Fuel - 20BA'] = lattice