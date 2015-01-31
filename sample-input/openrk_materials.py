"""
This file writes all of the materials data (multi-group nuclear 
cross-sections) for the LRA diffusion
benchmark problem to an HDF5 file. The script uses the h5py Python package
to interact with the HDF5 file format. This may be a good example for those
wishing ot write their nuclear data to an HDF5 file to import using the
OpenMOC 'materialize' Python module.
"""

from openrk import FunctionalMaterial
from openrk import Material

# Create a Python dictionary to store LRA multi-group cross-sections
dataset = {}

materials = []

dataset['Energy Groups'] = 2
dataset['Delayed Groups'] = 2
dataset['Materials'] = {}

lra_materials = dataset['Materials']

###############################################################################
################################   region 1    ################################
###############################################################################

# Create material
m = FunctionalMaterial(name='region 1')

# Set general material properties
m.setNumTimeSteps(3)
m.setTimeSteps([0, 2, 3])
m.setNumEnergyGroups(2)
m.setNumDelayedGroups(2)

# Set specific material properties
m.setSigmaA([[0.008252, 0.1003], [0.008252, 0.1003], [0.008252, 0.1003]])
m.setSigmaT([[0.2656, 1.5798], [0.2656, 1.5798], [0.2656, 1.5798]])
m.setSigmaF([[0.002, 0.05], [0.002, 0.05], [0.002, 0.05]])
m.setNuSigmaF([[0.004602, 0.1091], [0.004602, 0.1091], [0.004602, 0.1091]])
m.setSigmaS([[0.231892, 0.02533, 0.00, 1.47948],[0.231892, 0.02533, 0.00, 1.47948], [0.231892, 0.02533, 0.00, 1.47948]])
m.setChi([[1.0, 0.0], [1.0, 0.0], [1.0, 0.0]])
m.setDecayConstants([0.00654, 1.35])
m.setDelayedFractions([0.0054, 0.001087])
m.setVelocity([3e7, 3e5])
m.setDopplerCoefficients([1.e-3, 0.0])

# add material to materials list
materials.append(m)

###############################################################################
################################   region 2    ################################
###############################################################################

# Create material
m = FunctionalMaterial(name='region 2')

# Set general material properties
m.setNumTimeSteps(3)
m.setTimeSteps([0, 2, 3])
m.setNumEnergyGroups(2)
m.setNumDelayedGroups(2)

# Set specific material properties
m.setSigmaA([[0.007181, 0.07047], [0.007181, 0.07047], [0.007181, 0.07047]])
m.setSigmaT([[0.2629, 1.7525], [0.2629, 1.7525], [0.2629, 1.7525]])
m.setSigmaF([[0.002, 0.045], [0.002, 0.045], [0.002, 0.045]])
m.setNuSigmaF([[0.004609, 0.08675], [0.004609, 0.08675], [0.004609, 0.08675]])
m.setSigmaS([[0.22792, 0.02767, 0.00, 1.68201], [0.22792, 0.02767, 0.00, 1.68201], [0.22792, 0.02767, 0.00, 1.68201]])
m.setChi([[1.0, 0.0], [1.0, 0.0], [1.0, 0.0]])
m.setDecayConstants([0.00654, 1.35])
m.setDelayedFractions([0.0054, 0.001087])
m.setVelocity([3e7, 3e5])
m.setDopplerCoefficients([1.e-3, 0.0])

# add material to materials list
materials.append(m)

###############################################################################
################################   region 3    ################################
###############################################################################

# Create material
m = FunctionalMaterial(name='region 3')

# Set general material properties
m.setNumTimeSteps(3)
m.setTimeSteps([0, 2, 3])
m.setNumEnergyGroups(2)
m.setNumDelayedGroups(2)

# Set specific material properties
m.setSigmaA([[0.008002, 0.08344], [0.008002, 0.08344], [0.008002, 0.08344]])
m.setSigmaT([[0.2648, 1.5941], [0.2648, 1.5941], [0.2648, 1.5941]])
m.setSigmaF([[0.002, 0.045], [0.002, 0.045], [0.002, 0.045]])
m.setNuSigmaF([[0.004663, 0.1021], [0.004663, 0.1021], [0.004663, 0.1021]])
m.setSigmaS([[0.230502, 0.02617, 0.00, 1.510639], [0.230502, 0.02617, 0.00, 1.510639], [0.230502, 0.02617, 0.00, 1.510639]])
m.setChi([[1.0, 0.0], [1.0, 0.0], [1.0, 0.0]])
m.setDecayConstants([0.00654, 1.35])
m.setDelayedFractions([0.0054, 0.001087])
m.setVelocity([3e7, 3e5])
m.setDopplerCoefficients([1.e-3, 0.0])

# add material to materials list
materials.append(m)

###############################################################################
################################   region 4    ################################
###############################################################################

# Create material
m = FunctionalMaterial(name='region 4')

# Set general material properties
m.setNumTimeSteps(3)
m.setTimeSteps([0, 2, 3])
m.setNumEnergyGroups(2)
m.setNumDelayedGroups(2)

# Set specific material properties
m.setSigmaA([[0.008002, 0.073324], [0.008002, 0.073324], [0.008002, 0.073324]])
m.setSigmaT([[0.2648, 1.5941], [0.2648, 1.5941], [0.2648, 1.5941]])
m.setSigmaF([[0.002, 0.045], [0.002, 0.045], [0.002, 0.045]])
m.setNuSigmaF([[0.004663, 0.1021], [0.004663, 0.1021], [0.004663, 0.1021]])
m.setSigmaS([[0.230462, 0.02617, 0.00, 1.520789], [0.230462, 0.02617, 0.00, 1.520789], [0.230462, 0.02617, 0.00, 1.520789]])
m.setChi([[1.0, 0.0], [1.0, 0.0], [1.0, 0.0]])
m.setDecayConstants([0.00654, 1.35])
m.setDelayedFractions([0.0054, 0.001087])
m.setVelocity([3e7, 3e5])
m.setDopplerCoefficients([1.e-3, 0.0])

# add material to materials list
materials.append(m)

###############################################################################
################################   region 5    ################################
###############################################################################

# Create material
m = FunctionalMaterial(name='region 5')

# Set general material properties
m.setNumTimeSteps(3)
m.setTimeSteps([0, 2, 3])
m.setNumEnergyGroups(2)
m.setNumDelayedGroups(2)

# Set specific material properties
m.setSigmaA([[0.008002, 0.08344], [0.008002, 0.08344], [0.008002, 0.08344]])
m.setSigmaT([[0.2648, 1.5941], [0.2648, 1.5941], [0.2648, 1.5941]])
m.setSigmaF([[0.002, 0.045], [0.002, 0.045], [0.002, 0.045]])
m.setNuSigmaF([[0.004663, 0.1021], [0.004663, 0.1021], [0.004663, 0.1021]])
m.setSigmaS([[0.230462, 0.02617, 0.00, 1.510672], [0.230462, 0.02617, 0.00, 1.510672], [0.230462, 0.02617, 0.00, 1.510672]])
m.setChi([[1.0, 0.0], [1.0, 0.0], [1.0, 0.0]])
m.setDecayConstants([0.00654, 1.35])
m.setDelayedFractions([0.0054, 0.001087])
m.setVelocity([3e7, 3e5])
m.setDopplerCoefficients([1.e-3, 0.0])

# add material to materials list
materials.append(m)

###############################################################################
################################   region 6    ################################
###############################################################################

# Create material
m = FunctionalMaterial(name='region 6')

# Set general material properties
m.setNumTimeSteps(3)
m.setTimeSteps([0, 2, 3])
m.setNumEnergyGroups(2)
m.setNumDelayedGroups(2)

# Set specific material properties
m.setSigmaA([[0.0006034, 0.01911], [0.0006034, 0.01911], [0.0006034, 0.01911]])
m.setSigmaT([[0.2652, 2.0938], [0.2652, 2.0938], [0.2652, 2.0938]])
m.setSigmaF([[0.0, 0.0], [0.0, 0.0], [0.0, 0.0]])
m.setNuSigmaF([[0.0, 0.0], [0.0, 0.0], [0.0, 0.0]])
m.setSigmaS([[0.216931, 0.04754, 0.00, 2.074676], [0.216931, 0.04754, 0.00, 2.074676], [0.216931, 0.04754, 0.00, 2.074676]])
m.setChi([[1.0, 0.0], [1.0, 0.0], [1.0, 0.0]])
m.setDecayConstants([0.00654, 1.35])
m.setDelayedFractions([0.0054, 0.001087])
m.setVelocity([3e7, 3e5])
m.setDopplerCoefficients([1.e-3, 0.0])

# add material to materials list
materials.append(m)
