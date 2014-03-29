from pincells import *
from templates import *

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

cells = templates['0BA']

# Fuel Pins
univ_id = pincells['1.6% Fuel']['Universe']
cells = [univ_id if x is masks['Fuel Pin'] else x for x in cells]

# Guide Tubes
univ_id = pincells['Guide Tube']['Universe']
cells = [univ_id if x is masks['Guide Tube'] else x for x in cells]

# Instrument Tubes
univ_id = pincells['Instrument Tube']['Universe']
cells = [univ_id if x is masks['Instrument Tube'] else x for x in cells]

lattice = Lattice(id=universe_id(), width_x=pin_pitch, width_y=pin_pitch)
lattice.setLatticeCells(cells)

lattices['1.6% Fuel - 0BA'] = lattice


###############################################################################
####################   2.4% Enriched Fuel Assembly   ##########################
###############################################################################

cells = templates['0BA']

# Fuel Pins
univ_id = pincells['2.4% Fuel']['Universe']
cells = [univ_id if x is masks['Fuel Pin'] else x for x in cells]

# Guide Tubes
univ_id = pincells['Guide Tube']['Universe']
cells = [univ_id if x is masks['Guide Tube'] else x for x in cells]

# Instrument Tubes
univ_id = pincells['Instrument Tube']['Universe']
cells = [univ_id if x is masks['Instrument Tube'] else x for x in cells]

lattice = Lattice(id=universe_id(), width_x=pin_pitch, width_y=pin_pitch)
lattice.setLatticeCells(cells)

lattices['2.4% Fuel- 0BA'] = lattice


###############################################################################
################   2.4% Enriched Fuel Assembly w/ 12 BAs  #####################
###############################################################################

cells = templates['12BA']

# Fuel Pins
univ_id = pincells['2.4% Fuel']['Universe']
cells = [univ_id if x is masks['Fuel Pin'] else x for x in cells]

# Guide Tubes
univ_id = pincells['Guide Tube']['Universe']
cells = [univ_id if x is masks['Guide Tube'] else x for x in cells]

# Instrument Tubes
univ_id = pincells['Instrument Tube']['Universe']
cells = [univ_id if x is masks['Instrument Tube'] else x for x in cells]

# Burnable Absorbers
univ_id = pincells['Burnable Absorber']['Universe']
cells = [univ_id if x is masks['Burnable Absorber'] else x for x in cells]

lattice = Lattice(id=universe_id(), width_x=pin_pitch, width_y=pin_pitch)
lattice.setLatticeCells(cells)

lattices['2.4% Fuel - 12BA'] = lattice


###############################################################################
################   2.4% Enriched Fuel Assembly w/ 16 BAs  #####################
###############################################################################

cells = templates['16BA']

# Fuel Pins
univ_id = pincells['2.4% Fuel']['Universe']
cells = [univ_id if x is masks['Fuel Pin'] else x for x in cells]

# Guide Tubes
univ_id = pincells['Guide Tube']['Universe']
cells = [univ_id if x is masks['Guide Tube'] else x for x in cells]

# Instrument Tubes
univ_id = pincells['Instrument Tube']['Universe']
cells = [univ_id if x is masks['Instrument Tube'] else x for x in cells]

# Burnable Absorbers
univ_id = pincells['Burnable Absorber']['Universe']
cells = [univ_id if x is masks['Burnable Absorber'] else x for x in cells]

lattice = Lattice(id=universe_id(), width_x=pin_pitch, width_y=pin_pitch)
lattice.setLatticeCells(cells)

lattices['2.4% Fuel - 16BA'] = lattice


###############################################################################
#####################   3.1% Enriched Fuel Assembly  ##########################
###############################################################################

cells = templates['0BA']

# Fuel Pins
univ_id = pincells['3.1% Fuel']['Universe']
cells = [univ_id if x is masks['Fuel Pin'] else x for x in cells]

# Guide Tubes
univ_id = pincells['Guide Tube']['Universe']
cells = [univ_id if x is masks['Guide Tube'] else x for x in cells]

# Instrument Tubes
univ_id = pincells['Instrument Tube']['Universe']
cells = [univ_id if x is masks['Instrument Tube'] else x for x in cells]

lattice = Lattice(id=universe_id(), width_x=pin_pitch, width_y=pin_pitch)
lattice.setLatticeCells(cells)

lattices['3.1% Fuel- 0BA'] = lattice


###############################################################################
################   3.1% Enriched Fuel Assembly w/ 6 BAs  ######################
###############################################################################

cells = templates['6BA']

# Fuel Pins
univ_id = pincells['2.4% Fuel']['Universe']
cells = [univ_id if x is masks['Fuel Pin'] else x for x in cells]

# Guide Tubes
univ_id = pincells['Guide Tube']['Universe']
cells = [univ_id if x is masks['Guide Tube'] else x for x in cells]

# Instrument Tubes
univ_id = pincells['Instrument Tube']['Universe']
cells = [univ_id if x is masks['Instrument Tube'] else x for x in cells]

# Burnable Absorbers
univ_id = pincells['Burnable Absorber']['Universe']
cells = [univ_id if x is masks['Burnable Absorber'] else x for x in cells]

lattice = Lattice(id=universe_id(), width_x=pin_pitch, width_y=pin_pitch)
lattice.setLatticeCells(cells)

lattices['3.1% Fuel - 6BA'] = lattice


###############################################################################
################   3.1% Enriched Fuel Assembly w/ 15 BAs  #####################
###############################################################################

cells = templates['15BA']

# Fuel Pins
univ_id = pincells['3.1% Fuel']['Universe']
cells = [univ_id if x is masks['Fuel Pin'] else x for x in cells]

# Guide Tubes
univ_id = pincells['Guide Tube']['Clad']
cells = [univ_id if x is masks['Guide Tube'] else x for x in cells]

# Instrument Tubes
univ_id = pincells['Instrument Tube']['Universe']
cells = [univ_id if x is masks['Instrument Tube'] else x for x in cells]

# Burnable Absorbers
univ_id = pincells['Burnable Absorber']['Universe']
cells = [univ_id if x is masks['Burnable Absorber'] else x for x in cells]

lattice = Lattice(id=universe_id(), width_x=pin_pitch, width_y=pin_pitch)
lattice.setLatticeCells(cells)

lattices['3.1% Fuel - 15BA'] = lattice


###############################################################################
################   3.1% Enriched Fuel Assembly w/ 16 BAs  #####################
###############################################################################

cells = templates['16BA']

# Fuel Pins
univ_id = pincells['3.1% Fuel']['Universe']
cells = [univ_id if x is masks['Fuel Pin'] else x for x in cells]

# Guide Tubes
univ_id = pincells['Guide Tube']['Universe']
cells = [univ_id if x is masks['Guide Tube'] else x for x in cells]

# Instrument Tubes
univ_id = pincells['Instrument Tube']['Universe']
cells = [univ_id if x is masks['Instrument Tube'] else x for x in cells]

# Burnable Absorbers
univ_id = pincells['Burnable Absorber']['Universe']
cells = [univ_id if x is masks['Burnable Absorber'] else x for x in cells]

lattice = Lattice(id=universe_id(), width_x=pin_pitch, width_y=pin_pitch)
lattice.setLatticeCells(cells)

lattices['3.1% Fuel - 16BA'] = lattice


###############################################################################
################   3.1% Enriched Fuel Assembly w/ 20 BAs  #####################
###############################################################################

cells = templates['20BA']

# Fuel Pins
univ_id = pincells['3.1% Fuel']['Universe']
cells = [univ_id if x is masks['Fuel Pin'] else x for x in cells]

# Guide Tubes
univ_id = pincells['Guide Tube']['Universe']
cells = [univ_id if x is masks['Guide Tube'] else x for x in cells]

# Instrument Tubes
univ_id = pincells['Instrument Tube']['Universe']
cells = [univ_id if x is masks['Instrument Tube'] else x for x in cells]

# Burnable Absorbers
univ_id = pincells['Burnable Absorber']['Universe']
cells = [univ_id if x is masks['Burnable Absorber'] else x for x in cells]

lattice = Lattice(id=universe_id(), width_x=pin_pitch, width_y=pin_pitch)
lattice.setLatticeCells(cells)

lattices['3.1% Fuel - 20BA'] = lattice