# NOTE - These are for Cycle 1. The 12BA configuration is different for Cycle 2.
# NOTE - These are all taken at axial locations above the dashpot.


###############################################################################
####################   BEAVRS Fuel Assembly Parameters   ######################
###############################################################################

pin_pitch = 1.25984                                        # centimeters
lattice_pin_dims = 17.                                     # 17 pins
lattice_width = lattice_pin_dims * pin_pitch               # Lattice width [cm]


###############################################################################
#########################   Pin Cell Universe Masks   #########################
###############################################################################

# Keys are pin cell string names. Values are Universe IDs for each pin cell,
masks = dict()

masks['Fuel Pin'] = 1
masks['Guide Tube'] = 2
masks['Instrument Tube'] = 3
masks['Burnable Absorber'] = 4


###############################################################################
#################   BEAVRS Fuel Assembly Universe Templates   #################
###############################################################################

templates = dict()

templates['0BA'] = [[1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1],
                    [1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1],
                    [1,1,1,1,1,2,1,1,2,1,1,2,1,1,1,1,1],
                    [1,1,1,2,1,1,1,1,1,1,1,1,1,2,1,1,1],
                    [1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1],
                    [1,1,2,1,1,2,1,1,2,1,1,2,1,1,2,1,1],
                    [1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1],
                    [1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1],
                    [1,1,2,1,1,2,1,1,3,1,1,2,1,1,2,1,1],
                    [1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1],
                    [1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1],
                    [1,1,2,1,1,2,1,1,2,1,1,2,1,1,2,1,1],
                    [1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1],
                    [1,1,1,2,1,1,1,1,1,1,1,1,1,2,1,1,1],
                    [1,1,1,1,1,2,1,1,2,1,1,2,1,1,1,1,1],
                    [1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1],
                    [1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1]]

templates['4BA'] = [[1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1],
                    [1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1],
                    [1,1,1,1,1,2,1,1,2,1,1,2,1,1,1,1,1],
                    [1,1,1,4,1,1,1,1,1,1,1,1,1,4,1,1,1],
                    [1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1],
                    [1,1,2,1,1,2,1,1,2,1,1,2,1,1,2,1,1],
                    [1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1],
                    [1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1],
                    [1,1,2,1,1,2,1,1,3,1,1,2,1,1,2,1,1],
                    [1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1],
                    [1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1],
                    [1,1,2,1,1,2,1,1,2,1,1,2,1,1,2,1,1],
                    [1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1],
                    [1,1,1,4,1,1,1,1,1,1,1,1,1,4,1,1,1],
                    [1,1,1,1,1,2,1,1,2,1,1,2,1,1,1,1,1],
                    [1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1],
                    [1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1]]

templates['6BA'] = [[1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1],
                    [1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1],
                    [1,1,1,1,1,4,1,1,2,1,1,4,1,1,1,1,1],
                    [1,1,1,4,1,1,1,1,1,1,1,1,1,4,1,1,1],
                    [1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1],
                    [1,1,4,1,1,2,1,1,2,1,1,2,1,1,4,1,1],
                    [1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1],
                    [1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1],
                    [1,1,2,1,1,2,1,1,3,1,1,2,1,1,2,1,1],
                    [1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1],
                    [1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1],
                    [1,1,2,1,1,2,1,1,2,1,1,2,1,1,2,1,1],
                    [1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1],
                    [1,1,1,2,1,1,1,1,1,1,1,1,1,2,1,1,1],
                    [1,1,1,1,1,2,1,1,2,1,1,2,1,1,1,1,1],
                    [1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1],
                    [1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1]]

templates['8BA'] = [[1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1],
                    [1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1],
                    [1,1,1,1,1,2,1,1,2,1,1,2,1,1,1,1,1],
                    [1,1,1,4,1,1,1,1,1,1,1,1,1,4,1,1,1],
                    [1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1],
                    [1,1,2,1,1,2,1,1,4,1,1,2,1,1,2,1,1],
                    [1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1],
                    [1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1],
                    [1,1,2,1,1,4,1,1,3,1,1,4,1,1,2,1,1],
                    [1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1],
                    [1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1],
                    [1,1,2,1,1,2,1,1,4,1,1,2,1,1,2,1,1],
                    [1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1],
                    [1,1,1,4,1,1,1,1,1,1,1,1,1,4,1,1,1],
                    [1,1,1,1,1,2,1,1,2,1,1,2,1,1,1,1,1],
                    [1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1],
                    [1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1]]

templates['12BA'] = [[1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1],
                     [1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1],
                     [1,1,1,1,1,4,1,1,2,1,1,4,1,1,1,1,1],
                     [1,1,1,4,1,1,1,1,1,1,1,1,1,4,1,1,1],
                     [1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1],
                     [1,1,4,1,1,2,1,1,2,1,1,2,1,1,4,1,1],
                     [1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1],
                     [1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1],
                     [1,1,2,1,1,2,1,1,3,1,1,2,1,1,2,1,1],
                     [1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1],
                     [1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1],
                     [1,1,4,1,1,2,1,1,2,1,1,2,1,1,4,1,1],
                     [1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1],
                     [1,1,1,4,1,1,1,1,1,1,1,1,1,4,1,1,1],
                     [1,1,1,1,1,4,1,1,2,1,1,4,1,1,1,1,1],
                     [1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1],
                     [1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1]]

templates['15BA'] = [[1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1],
                     [1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1],
                     [1,1,1,1,1,4,1,1,4,1,1,4,1,1,1,1,1],
                     [1,1,1,4,1,1,1,1,1,1,1,1,1,2,1,1,1],
                     [1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1],
                     [1,1,4,1,1,4,1,1,4,1,1,4,1,1,2,1,1],
                     [1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1],
                     [1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1],
                     [1,1,4,1,1,4,1,1,3,1,1,4,1,1,2,1,1],
                     [1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1],
                     [1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1],
                     [1,1,4,1,1,4,1,1,4,1,1,4,1,1,2,1,1],
                     [1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1],
                     [1,1,1,2,1,1,1,1,1,1,1,1,1,2,1,1,1],
                     [1,1,1,1,1,2,1,1,2,1,1,2,1,1,1,1,1],
                     [1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1],
                     [1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1]]

templates['16BA'] = [[1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1],
                     [1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1],
                     [1,1,1,1,1,4,1,1,4,1,1,4,1,1,1,1,1],
                     [1,1,1,4,1,1,1,1,1,1,1,1,1,4,1,1,1],
                     [1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1],
                     [1,1,4,1,1,2,1,1,2,1,1,2,1,1,4,1,1],
                     [1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1],
                     [1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1],
                     [1,1,4,1,1,2,1,1,3,1,1,2,1,1,4,1,1],
                     [1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1],
                     [1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1],
                     [1,1,4,1,1,2,1,1,2,1,1,2,1,1,4,1,1],
                     [1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1],
                     [1,1,1,4,1,1,1,1,1,1,1,1,1,4,1,1,1],
                     [1,1,1,1,1,4,1,1,4,1,1,4,1,1,1,1,1],
                     [1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1],
                     [1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1]]

templates['20BA'] = [[1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1],
                     [1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1],
                     [1,1,1,1,1,4,1,1,4,1,1,4,1,1,1,1,1],
                     [1,1,1,4,1,1,1,1,1,1,1,1,1,4,1,1,1],
                     [1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1],
                     [1,1,4,1,1,4,1,1,2,1,1,4,1,1,4,1,1],
                     [1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1],
                     [1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1],
                     [1,1,4,1,1,2,1,1,3,1,1,2,1,1,4,1,1],
                     [1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1],
                     [1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1],
                     [1,1,4,1,1,4,1,1,2,1,1,4,1,1,4,1,1],
                     [1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1],
                     [1,1,1,4,1,1,1,1,1,1,1,1,1,4,1,1,1],
                     [1,1,1,1,1,4,1,1,4,1,1,4,1,1,1,1,1],
                     [1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1],
                     [1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1]]