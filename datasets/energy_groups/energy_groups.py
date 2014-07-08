from infermc.build import EnergyGroups
import numpy as np

group_structures = dict()


################################################################################
######################################  CASMO  #################################
################################################################################

# Create a sub-dictionary for the CASMO energy group structures
casmo = dict()


# 2-group structure
casmo['2-group'] = EnergyGroups()
group_edges = np.array([1.e-11, 0.625e-6, 10.])
casmo['2-group'].setGroupEdges(group_edges)


# 8-group structure
casmo['8-group'] = EnergyGroups()
group_edges = np.array([1.e-11, 0.03e-6, 0.1e-6, 0.3e-6,
                        0.625e-6, 4.e-6, 5.5308e-3, 1., 20.])
casmo['8-group'].setGroupEdges(group_edges)


# 20-group structure
casmo['20-group'] = EnergyGroups()
group_edges = np.array([1.e-11, 0.03e-6, 0.1e-6, 0.19e-6,
                        0.3e-6, 0.35e-6, 0.625e-6, 2.e-6,
                        4.e-6, 6.525e-6, 6.825e-6, 16.e-6,
                        27.7e-6, 47.9e-6, 148.7e-6, 5.5308e-3,
                        19.305e-3, 1., 4.9659, 20.])
casmo['20-group'].setGroupEdges(group_edges)


# Store the sub-dictionary in the global dictionary
group_structures['CASMO'] = casmo