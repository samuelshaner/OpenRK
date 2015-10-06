import openrk as rk
import openrk.plotter as plotter

#################################################
###############  Create Materials ###############
#################################################

# create fuel 1 blade in
fuel1bin = rk.FunctionalMaterial()
fuel1bin.setTimeSteps([0.0, 2.0, 3.0])
fuel1bin.setNumEnergyGroups(2)
fuel1bin.setNumDelayedGroups(2)
fuel1bin.setSigmaA([[0.008252, 0.1003], [0.008252, 0.1003], [0.008252, 0.1003]])
fuel1bin.setDifCoef([[1.255, 0.211], [1.255, 0.211], [1.255, 0.211]])
fuel1bin.setSigmaT([[1.0/(3.0*1.255), 1.0/(3.0*0.211)], [1.0/(3.0*1.255), 1.0/(3.0*0.211)], [1.0/(3.0*1.255), 1.0/(3.0*0.211)]])
fuel1bin.setNuSigmaF([[0.004602, 0.1091], [0.004602, 0.1091], [0.004602, 0.1091]])
fuel1bin.setSigmaF([[0.004602/2.43, 0.1091/2.43], [0.004602/2.43, 0.1091/2.43], [0.004602/2.43, 0.1091/2.43]])
fuel1bin.setChi([[1.0, 0.0], [1.0, 0.0], [1.0, 0.0]])
fuel1bin.setSigmaS([[0.0, 0.02533, 0.0, 0.0], [0.0, 0.02533, 0.0, 0.0], [0.0, 0.02533, 0.0, 0.0]])
fuel1bin.setDopplerCoefficients([3.034e-3, 0.0])
fuel1bin.setEnergyPerFission(3.204e-11)
fuel1bin.setVelocity([[3.e7, 3.e5], [3.e7, 3.e5], [3.e7, 3.e5]])
fuel1bin.setTemperatureConversionFactor(3.83e-11)
fuel1bin.setDecayConstant([0.00654, 1.35])
fuel1bin.setDelayedFraction([[0.0054, 0.001087], [0.0054, 0.001087], [0.0054, 0.001087]])

# create fuel 1 blade out
fuel1bo = rk.FunctionalMaterial()
fuel1bo.setTimeSteps([0.0, 2.0, 3.0])
fuel1bo.setNumEnergyGroups(2)
fuel1bo.setNumDelayedGroups(2)
fuel1bo.setSigmaA([[0.007181, 0.07047], [0.007181, 0.07047], [0.007181, 0.07047]])
fuel1bo.setDifCoef([[1.268, 0.1902], [1.268, 0.1902], [1.268, 0.1902]])
fuel1bo.setSigmaT([[1.0/(3.0*1.268), 1.0/(3.0*0.1902)], [1.0/(3.0*1.268), 1.0/(3.0*0.1902)], [1.0/(3.0*1.268), 1.0/(3.0*0.1902)]])
fuel1bo.setNuSigmaF([[0.004609, 0.08675], [0.004609, 0.08675], [0.004609, 0.08675]])
fuel1bo.setSigmaF([[0.004609/2.43, 0.08675/2.43], [0.004609/2.43, 0.08675/2.43], [0.004609/2.43, 0.08675/2.43]])
fuel1bo.setChi([[1.0, 0.0], [1.0, 0.0], [1.0, 0.0]])
fuel1bo.setSigmaS([[0.0, 0.02767, 0.0, 0.0], [0.0, 0.02767, 0.0, 0.0], [0.0, 0.02767, 0.0, 0.0]])
fuel1bo.setDopplerCoefficients([3.034e-3, 0.0])
fuel1bo.setEnergyPerFission(3.204e-11)
fuel1bo.setVelocity([[3.e7, 3.e5], [3.e7, 3.e5], [3.e7, 3.e5]])
fuel1bo.setTemperatureConversionFactor(3.83e-11)
fuel1bo.setDecayConstant([0.00654, 1.35])
fuel1bo.setDelayedFraction([[0.0054, 0.001087], [0.0054, 0.001087], [0.0054, 0.001087]])

# create fuel 2 blade in
fuel2bin = rk.FunctionalMaterial()
fuel2bin.setTimeSteps([0.0, 2.0, 3.0])
fuel2bin.setNumEnergyGroups(2)
fuel2bin.setNumDelayedGroups(2)
fuel2bin.setSigmaA([[0.008002, 0.08344], [0.008002, 0.08344], [0.008002, 0.08344]])
fuel2bin.setDifCoef([[1.259, 0.2091], [1.259, 0.2091], [1.259, 0.2091]])
fuel2bin.setSigmaT([[1.0/(3.0*1.259), 1.0/(3.0*0.2091)], [1.0/(3.0*1.259), 1.0/(3.0*0.2091)], [1.0/(3.0*1.259), 1.0/(3.0*0.2091)]])
fuel2bin.setNuSigmaF([[0.004663, 0.1021], [0.004663, 0.1021], [0.004663, 0.1021]])
fuel2bin.setSigmaF([[0.004663/2.43, 0.1021/2.43], [0.004663/2.43, 0.1021/2.43], [0.004663/2.43, 0.1021/2.43]])
fuel2bin.setChi([[1.0, 0.0], [1.0, 0.0], [1.0, 0.0]])
fuel2bin.setSigmaS([[0.0, 0.02617, 0.0, 0.0], [0.0, 0.02617, 0.0, 0.0], [0.0, 0.02617, 0.0, 0.0]])
fuel2bin.setDopplerCoefficients([3.034e-3, 0.0])
fuel2bin.setEnergyPerFission(3.204e-11)
fuel2bin.setVelocity([[3.e7, 3.e5], [3.e7, 3.e5], [3.e7, 3.e5]])
fuel2bin.setTemperatureConversionFactor(3.83e-11)
fuel2bin.setDecayConstant([0.00654, 1.35])
fuel2bin.setDelayedFraction([[0.0054, 0.001087], [0.0054, 0.001087], [0.0054, 0.001087]])

# create fuel 2 blade out
fuel2bo = rk.FunctionalMaterial()
fuel2bo.setTimeSteps([0.0, 2.0, 3.0])
fuel2bo.setNumEnergyGroups(2)
fuel2bo.setNumDelayedGroups(2)
fuel2bo.setSigmaA([[0.008002, 0.073324], [0.008002, 0.073324], [0.008002, 0.073324]])
fuel2bo.setDifCoef([[1.259, 0.2091], [1.259, 0.2091], [1.259, 0.2091]])
fuel2bo.setSigmaT([[1.0/(3.0*1.259), 1.0/(3.0*0.2091)], [1.0/(3.0*1.259), 1.0/(3.0*0.2091)], [1.0/(3.0*1.259), 1.0/(3.0*0.2091)]])
fuel2bo.setNuSigmaF([[0.004663, 0.1021], [0.004663, 0.1021], [0.004663, 0.1021]])
fuel2bo.setSigmaF([[0.004663/2.43, 0.1021/2.43], [0.004663/2.43, 0.1021/2.43], [0.004663/2.43, 0.1021/2.43]])
fuel2bo.setChi([[1.0, 0.0], [1.0, 0.0], [1.0, 0.0]])
fuel2bo.setSigmaS([[0.0, 0.02617, 0.0, 0.0], [0.0, 0.02617, 0.0, 0.0], [0.0, 0.02617, 0.0, 0.0]])
fuel2bo.setDopplerCoefficients([3.034e-3, 0.0])
fuel2bo.setEnergyPerFission(3.204e-11)
fuel2bo.setVelocity([[3.e7, 3.e5], [3.e7, 3.e5], [3.e7, 3.e5]])
fuel2bo.setTemperatureConversionFactor(3.83e-11)
fuel2bo.setDecayConstant([0.00654, 1.35])
fuel2bo.setDelayedFraction([[0.0054, 0.001087], [0.0054, 0.001087], [0.0054, 0.001087]])

# create fuel 2 blade in then out
fuel2bino = rk.FunctionalMaterial()
fuel2bino.setTimeSteps([0.0, 2.0, 3.0])
fuel2bino.setNumEnergyGroups(2)
fuel2bino.setNumDelayedGroups(2)
#fuel2bino.setSigmaA([[0.008002, 0.08344], [0.008002, 0.08344], [0.008002, 0.08344]])
fuel2bino.setSigmaA([[0.008002, 0.08344], [0.008002, 0.073324], [0.008002, 0.073324]])
fuel2bino.setDifCoef([[1.259, 0.2091], [1.259, 0.2091], [1.259, 0.2091]])
fuel2bino.setSigmaT([[1.0/(3.0*1.259), 1.0/(3.0*0.2091)], [1.0/(3.0*1.259), 1.0/(3.0*0.2091)], [1.0/(3.0*1.259), 1.0/(3.0*0.2091)]])
fuel2bino.setNuSigmaF([[0.004663, 0.1021], [0.004663, 0.1021], [0.004663, 0.1021]])
fuel2bino.setSigmaF([[0.004663/2.43, 0.1021/2.43], [0.004663/2.43, 0.1021/2.43], [0.004663/2.43, 0.1021/2.43]])
fuel2bino.setChi([[1.0, 0.0], [1.0, 0.0], [1.0, 0.0]])
fuel2bino.setSigmaS([[0.0, 0.02617, 0.0, 0.0], [0.0, 0.02617, 0.0, 0.0], [0.0, 0.02617, 0.0, 0.0]])
fuel2bino.setDopplerCoefficients([3.034e-3, 0.0])
fuel2bino.setEnergyPerFission(3.204e-11)
fuel2bino.setVelocity([[3.e7, 3.e5], [3.e7, 3.e5], [3.e7, 3.e5]])
fuel2bino.setTemperatureConversionFactor(3.83e-11)
fuel2bino.setDecayConstant([0.00654, 1.35])
fuel2bino.setDelayedFraction([[0.0054, 0.001087], [0.0054, 0.001087], [0.0054, 0.001087]])


# create reflector
reflector = rk.Material()
reflector.setNumEnergyGroups(2)
reflector.setSigmaA([0.000603, 0.01911])
reflector.setDifCoef([1.257, 0.1592])
reflector.setSigmaT([1.0/(3.0*1.257), 1.0/(3.0*0.1592)])
reflector.setNuSigmaF([0.0, 0.0])
reflector.setSigmaF([0.0, 0.0])
reflector.setChi([1.0, 0.0])
reflector.setSigmaS([0.0, 0.04754, 0.0, 0.0])
reflector.setVelocity([3.e7, 3.e5])

#################################################
###############  Create Geometry  ###############
#################################################

geometry = rk.GeometryDiffusion(165.0, 165.0)
geometry.setBoundary(0, 1)
geometry.setBoundary(1, 0)
geometry.setBoundary(2, 1)
geometry.setBoundary(3, 0)
geometry.setAmpMeshDimensions(11,11)
geometry.setShapeMeshDimensions(11,11)

nx = 11

for i in xrange(9, 11):
    for cell_id in xrange(i*nx, (i+1)*nx):
        geometry.setMaterial(reflector, cell_id)

for i in xrange(9):
    for cell_id in xrange(i*nx+7, (i+1)*nx):
        geometry.setMaterial(reflector, cell_id)

for i in xrange(7, 9):
    for cell_id in xrange(i*nx, i*nx+7):
        geometry.setMaterial(fuel2bin, cell_id)

for i in xrange(5, 7):
    for cell_id in xrange(i*nx+7, i*nx+9):
        geometry.setMaterial(fuel2bino, cell_id)

for i in xrange(5):
    for cell_id in xrange(i*nx+7, i*nx+9):
        geometry.setMaterial(fuel2bin, cell_id)

geometry.setMaterial(fuel1bo, 6*nx)
geometry.setMaterial(fuel1bo, 5*nx)

for i in xrange(5, 7):
    for cell_id in xrange(i*nx+5, i*nx+7):
        geometry.setMaterial(fuel1bo, cell_id)

geometry.setMaterial(fuel1bo, 0)
for cell_id in xrange(5, 7):
    geometry.setMaterial(fuel1bo, cell_id)

geometry.setMaterial(fuel2bo, 7*nx+7)

for i in xrange(5, 7):
    for cell_id in xrange(i*nx+1, i*nx+5):
        geometry.setMaterial(fuel1bin, cell_id)

for i in xrange(1, 5):
    for cell_id in xrange(i*nx, i*nx+7):
        geometry.setMaterial(fuel1bin, cell_id)

for cell_id in xrange(1, 5):
    geometry.setMaterial(fuel1bin, cell_id)

geometry = geometry.uniformRefine(5,5,1)
geometry.uniquifyMaterials()
geometry.generateCellMap()

#################################################
###############   Create Solver   ###############
#################################################

solver = rk.SolverDiffusionTheta(geometry)
solver.setBuckling(1.e-4)
rk.setNumThreads(4)
solver.setStepSize(1.e-2)
solver.setFluxSolveTolerance(1.e-7)
solver.computeInitialFlux(1.e-7)

for i in range(300):
  solver.takeStep()

#solver.computeFrequency()
  
#plotter.plot_flux(solver)
#plotter.plot_precursors(solver)
#plotter.plot_frequency(solver)
