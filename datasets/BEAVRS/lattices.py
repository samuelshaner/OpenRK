import copy
from opencsg import Lattice
from datasets.BEAVRS.pin_cells import pin_cells
from datasets.BEAVRS.templates import *


# Keys are Lattice string names and values are OpenCSG Lattices.
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
  universe = pin_cells['1.6% Fuel Pin']
  row = [universe if x == masks['Fuel Pin'] else x for x in row]

  # Guide Tubes
  universe = pin_cells['Guide Tube']
  row = [universe if x == masks['Guide Tube'] else x for x in row]

  # Instrument Tubes
  universe = pin_cells['Instrument Tube']
  row = [universe if x == masks['Instrument Tube'] else x for x in row]

  cells[i] = row

lattice = Lattice(name='1.6% Fuel - 0BA')
lattices[lattice._name] = lattice

lattice.setDimension((17, 17))
lattice.setWidth((pin_pitch, pin_pitch))
lattice.setUniverses(cells)



###############################################################################
####################   2.4% Enriched Fuel Assembly   ##########################
###############################################################################

cells = copy.deepcopy(templates['0BA'])

for i in range(len(cells)):

  row = cells[i]

  # Fuel Pins
  universe = pin_cells['2.4% Fuel Pin']
  row = [universe if x == masks['Fuel Pin'] else x for x in row]

  # Guide Tubes
  universe = pin_cells['Guide Tube']
  row = [universe if x == masks['Guide Tube'] else x for x in row]

  # Instrument Tubes
  universe = pin_cells['Instrument Tube']
  row = [universe if x == masks['Instrument Tube'] else x for x in row]

  cells[i] = row

lattice = Lattice(name='2.4% Fuel - 0BA')
lattices[lattice._name] = lattice

lattice.setDimension((17, 17))
lattice.setWidth((pin_pitch, pin_pitch))
lattice.setUniverses(cells)



###############################################################################
################   2.4% Enriched Fuel Assembly w/ 12 BAs  #####################
###############################################################################

cells = copy.deepcopy(templates['12BA'])

for i in range(len(cells)):

  row = cells[i]

  # Fuel Pins
  universe = pin_cells['2.4% Fuel Pin']
  row = [universe if x == masks['Fuel Pin'] else x for x in row]

  # Guide Tubes
  universe = pin_cells['Guide Tube']
  row = [universe if x == masks['Guide Tube'] else x for x in row]

  # Instrument Tubes
  universe = pin_cells['Instrument Tube']
  row = [universe if x == masks['Instrument Tube'] else x for x in row]

  # Burnable Absorbers
  universe = pin_cells['Burnable Absorber']
  row = [universe if x == masks['Burnable Absorber'] else x for x in row]

  cells[i] = row

lattice = Lattice(name='2.4% Fuel - 12BA')
lattices[lattice._name] = lattice

lattice.setDimension((17, 17))
lattice.setWidth((pin_pitch, pin_pitch))
lattice.setUniverses(cells)



###############################################################################
################   2.4% Enriched Fuel Assembly w/ 16 BAs  #####################
###############################################################################

cells = copy.deepcopy(templates['16BA'])

for i in range(len(cells)):

  row = cells[i]

  # Fuel Pins
  universe = pin_cells['2.4% Fuel Pin']
  row = [universe if x == masks['Fuel Pin'] else x for x in row]

  # Guide Tubes
  universe = pin_cells['Guide Tube']
  row = [universe if x == masks['Guide Tube'] else x for x in row]

  # Instrument Tubes
  universe = pin_cells['Instrument Tube']
  row = [universe if x == masks['Instrument Tube'] else x for x in row]

  # Burnable Absorbers
  universe = pin_cells['Burnable Absorber']
  row = [universe if x == masks['Burnable Absorber'] else x for x in row]

  cells[i] = row

lattice = Lattice(name='2.4% Fuel - 16BA')
lattices[lattice._name] = lattice

lattice.setDimension((17, 17))
lattice.setWidth((pin_pitch, pin_pitch))
lattice.setUniverses(cells)



###############################################################################
#####################   3.1% Enriched Fuel Assembly  ##########################
###############################################################################

cells = copy.deepcopy(templates['0BA'])

for i in range(len(cells)):

  row = cells[i]

  # Fuel Pins
  universe = pin_cells['3.1% Fuel Pin']
  row = [universe if x == masks['Fuel Pin'] else x for x in row]

  # Guide Tubes
  univ_id = pin_cells['Guide Tube']
  row = [univ_id if x == masks['Guide Tube'] else x for x in row]

  # Instrument Tubes
  universe = pin_cells['Instrument Tube']
  row = [universe if x == masks['Instrument Tube'] else x for x in row]

  cells[i] = row

lattice = Lattice(name='3.1% Fuel - 0BA')
lattices[lattice._name] = lattice

lattice.setDimension((17, 17))
lattice.setWidth((pin_pitch, pin_pitch))
lattice.setUniverses(cells)



###############################################################################
################   3.1% Enriched Fuel Assembly w/ 6 BAs  ######################
###############################################################################

cells = copy.deepcopy(templates['6BA'])

for i in range(len(cells)):

  row = cells[i]

  # Fuel Pins
  universe = pin_cells['3.1% Fuel Pin']
  row = [universe if x == masks['Fuel Pin'] else x for x in row]

  # Guide Tubes
  universe = pin_cells['Guide Tube']
  row = [universe if x == masks['Guide Tube'] else x for x in row]

  # Instrument Tubes
  universe = pin_cells['Instrument Tube']
  row = [universe if x == masks['Instrument Tube'] else x for x in row]

  # Burnable Absorbers
  universe = pin_cells['Burnable Absorber']
  row = [universe if x == masks['Burnable Absorber'] else x for x in row]

  cells[i] = row

lattice = Lattice(name='3.1% Fuel - 6BA')
lattices[lattice._name] = lattice

lattice.setDimension((17, 17))
lattice.setWidth((pin_pitch, pin_pitch))
lattice.setUniverses(cells)



###############################################################################
################   3.1% Enriched Fuel Assembly w/ 15 BAs  #####################
###############################################################################

cells = copy.deepcopy(templates['15BA'])

for i in range(len(cells)):

  row = cells[i]

  # Fuel Pins
  universe = pin_cells['3.1% Fuel Pin']
  row = [universe if x == masks['Fuel Pin'] else x for x in row]

  # Guide Tubes
  universe = pin_cells['Guide Tube']
  row = [universe if x == masks['Guide Tube'] else x for x in row]

  # Instrument Tubes
  universe = pin_cells['Instrument Tube']
  row = [universe if x == masks['Instrument Tube'] else x for x in row]

  # Burnable Absorbers
  universe = pin_cells['Burnable Absorber']
  row = [universe if x == masks['Burnable Absorber'] else x for x in row]

  cells[i] = row

lattice = Lattice(name='3.1% Fuel - 15BA')
lattices[lattice._name] = lattice

lattice.setDimension((17, 17))
lattice.setWidth((pin_pitch, pin_pitch))
lattice.setUniverses(cells)



###############################################################################
################   3.1% Enriched Fuel Assembly w/ 16 BAs  #####################
###############################################################################

cells = copy.deepcopy(templates['16BA'])

for i in range(len(cells)):

  row = cells[i]

  # Fuel Pins
  universe = pin_cells['3.1% Fuel Pin']
  row = [universe if x == masks['Fuel Pin'] else x for x in row]

  # Guide Tubes
  universe = pin_cells['Guide Tube']
  row = [universe if x == masks['Guide Tube'] else x for x in row]

  # Instrument Tubes
  universe = pin_cells['Instrument Tube']
  row = [universe if x == masks['Instrument Tube'] else x for x in row]

  # Burnable Absorbers
  universe = pin_cells['Burnable Absorber']
  row = [universe if x == masks['Burnable Absorber'] else x for x in row]

  cells[i] = row

lattice = Lattice(name='3.1% Fuel - 16BA')
lattices[lattice._name] = lattice

lattice.setDimension((17, 17))
lattice.setWidth((pin_pitch, pin_pitch))
lattice.setUniverses(cells)



###############################################################################
################   3.1% Enriched Fuel Assembly w/ 20 BAs  #####################
###############################################################################

cells = copy.deepcopy(templates['20BA'])

for i in range(len(cells)):

  row = cells[i]

  # Fuel Pins
  universe = pin_cells['3.1% Fuel Pin']
  row = [universe if x == masks['Fuel Pin'] else x for x in row]

  # Guide Tubes
  universe = pin_cells['Guide Tube']
  row = [universe if x == masks['Guide Tube'] else x for x in row]

  # Instrument Tubes
  universe = pin_cells['Instrument Tube']
  row = [universe if x == masks['Instrument Tube'] else x for x in row]

  # Burnable Absorbers
  universe = pin_cells['Burnable Absorber']
  row = [universe if x == masks['Burnable Absorber'] else x for x in row]

  cells[i] = row

lattice = Lattice(name='3.1% Fuel - 20BA')
lattices[lattice._name] = lattice

lattice.setDimension((17, 17))
lattice.setWidth((pin_pitch, pin_pitch))
lattice.setUniverses(cells)