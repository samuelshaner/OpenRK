import openrk.openrk as rk
import openrk.plotter as plotter
from time import sleep

# Create materials

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

# create fuel 2 blade in then out
fuel2bino = rk.FunctionalMaterial()
fuel2bino.setTimeSteps([0.0, 2.0, 3.0])
fuel2bino.setNumEnergyGroups(2)
fuel2bino.setNumDelayedGroups(2)
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

shape_mesh = rk.StructuredShapeMesh(width=165.0, height=165.0, depth=360.0, num_x=11, num_y=11, num_z=12)
shape_mesh.setBoundary(2, 1)
shape_mesh.setBoundary(3, 1)
shape_mesh.setBoundary(4, 1)
shape_mesh.setBoundary(5, 1)
shape_mesh.setNumAmpEnergyGroups(2)
shape_mesh.setNumShapeEnergyGroups(2)
shape_mesh.setNumDelayedGroups(2)
shape_mesh.setBuckling(0.0)
shape_mesh.setDelayedFractions([0.0054, 0.001087])
shape_mesh.setDecayConstants([0.00654, 1.35])
shape_mesh.initialize()
shape_mesh.setTemperature(300)

nx = shape_mesh.getNumX()
ny = shape_mesh.getNumY()
nz = shape_mesh.getNumZ()

# Set upper and lower reflector cells
for i in xrange(nx*ny):
    shape_mesh.setMaterial(reflector, i)
    shape_mesh.setMaterial(reflector, 11*nx*ny+i)

# set core cells
for k in xrange(1,11):
    for i in xrange(9, 11):
        for cell_id in xrange(i*nx, (i+1)*nx):
            cell = k*nx*ny + cell_id
            shape_mesh.setMaterial(reflector, cell)

    for i in xrange(9):
        for cell_id in xrange(i*nx+7, (i+1)*nx):
            cell = k*nx*ny + cell_id
            shape_mesh.setMaterial(reflector, cell)

    for i in xrange(7, 9):
        for cell_id in xrange(i*nx, i*nx+7):
            cell = k*nx*ny + cell_id
            shape_mesh.setMaterial(fuel2bin, cell)

    for i in xrange(5, 7):
        for cell_id in xrange(i*nx+7, i*nx+9):
            cell = k*nx*ny + cell_id
            shape_mesh.setMaterial(fuel2bino, cell)

    for i in xrange(5):
        for cell_id in xrange(i*nx+7, i*nx+9):
            cell = k*nx*ny + cell_id
            shape_mesh.setMaterial(fuel2bin, cell)

    shape_mesh.setMaterial(fuel1bo, k*nx*ny + 6*nx)
    shape_mesh.setMaterial(fuel1bo, k*nx*ny + 5*nx)

    for i in xrange(5, 7):
        for cell_id in xrange(i*nx+5, i*nx+7):
            cell = k*nx*ny + cell_id
            shape_mesh.setMaterial(fuel1bo, cell)

    shape_mesh.setMaterial(fuel1bo, k*nx*ny + 0)
    for cell_id in xrange(5, 7):
        cell = k*nx*ny + cell_id
        shape_mesh.setMaterial(fuel1bo, cell)

    shape_mesh.setMaterial(fuel2bo, k*nx*ny + 7*nx+7)

    for i in xrange(5, 7):
        for cell_id in xrange(i*nx+1, i*nx+5):
            cell = k*nx*ny + cell_id
            shape_mesh.setMaterial(fuel1bin, cell)

    for i in xrange(1, 5):
        for cell_id in xrange(i*nx, i*nx+7):
            cell = k*nx*ny + cell_id
            shape_mesh.setMaterial(fuel1bin, cell)

    for cell_id in xrange(1, 5):
        cell = k*nx*ny + cell_id
        shape_mesh.setMaterial(fuel1bin, cell)

#shape_mesh = shape_mesh.uniformRefine(3,3,4)
shape_mesh.uniquifyMaterials()

# Create and initialize the amplitude mesh
amp_mesh = rk.AmpMesh(width=165.0, height=165.0, depth=360.0, num_x=11, num_y=11, num_z=12)
amp_mesh.setNumAmpEnergyGroups(2)
amp_mesh.setNumShapeEnergyGroups(2)
amp_mesh.setNumDelayedGroups(2)
amp_mesh.setOpticallyThick(False)
amp_mesh.setBoundary(2, 1)
amp_mesh.setBoundary(3, 1)
amp_mesh.setBoundary(4, 1)
amp_mesh.setBoundary(5, 1)
amp_mesh.setBuckling(0.0)
amp_mesh.setDelayedFractions([0.0054, 0.001087])
amp_mesh.setDecayConstants([0.00654, 1.35])
amp_mesh.initialize()
amp_mesh.setShapeMesh(shape_mesh)
amp_mesh.setGroupStructure([0, 1, 2])
shape_mesh.setAmpMesh(amp_mesh)
shape_mesh.setGroupStructure([0, 1])


# Solve diffusion problem
solver = rk.Solver(shape_mesh, amp_mesh)
rk.setNumThreads(2)

transient = rk.Transient()
clock = rk.Clock(dt_inner=1.e-3, dt_outer=1.e-1)
transient.setClock(clock)
#transient.setOuterMethod(rk.CRANK_NICOLSON)
#transient.setInnerMethod(rk.CRANK_NICOLSON)
transient.setShapeMesh(shape_mesh)
transient.setAmpMesh(amp_mesh)
transient.setSolver(solver)
transient.setInitialPower(1.e-6)
transient.computeInitialShape()

#for i in xrange(300):
#    transient.takeOuterStep()
    #rk.plotter.plot_power(mesh, name='mesh-power-{:.4f}s'.format(mesh.get_clock().get_time('CURRENT')))
    #rk.plotter.plot_temperature(mesh, name='mesh-temp-{:.4f}s'.format(mesh.get_clock().get_time('CURRENT')))


plotter.plot_flux(shape_mesh, plane='xy')
plotter.plot_power(shape_mesh, plane='xy')
plotter.plot_temperature(shape_mesh, plane='yz')
plotter.plot_materials(shape_mesh, plane='xz')
plotter.plot_sigma_a(shape_mesh, plane='yz')
plotter.plot_precursor_conc(shape_mesh, plane='xz')
