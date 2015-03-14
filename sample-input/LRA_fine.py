import openrk as rk
import openrk.mesh
import openrk.material
import openrk.solver
import openrk.plotter
import numpy as np

mesh = rk.StructuredMesh(165.0, 165.0, 165, 165)
mesh.setBoundary(2, 1)
mesh.setBoundary(3, 1)

# create fuel 1 blade in
fuel1bin = rk.material.Material(name='fuel 1 blade in', num_energy_groups=2)
fuel1bin.setSigmaA([0.0083775, 0.1003211])
fuel1bin.setDifCoef([1.255, 0.211])
fuel1bin.setNuSigmaF([0.004602, 0.1091])
fuel1bin.setChi([1.0, 0.0])
fuel1bin.setSigmaS(np.array([[0.0, 0.02533],[0.0, 0.0]]))

# create fuel 1 blade out
fuel1bo = rk.material.Material(name='fuel 1 blade out', num_energy_groups=2)
fuel1bo.setSigmaA([0.0073078, 0.07048902])
fuel1bo.setDifCoef([1.268, 0.1902])
fuel1bo.setNuSigmaF([0.004609, 0.08675])
fuel1bo.setChi([1.0, 0.0])
fuel1bo.setSigmaS(np.array([[0.0, 0.02767],[0.0, 0.0]]))

# create fuel 2 blade in
fuel2bin = rk.material.Material(name='fuel 2 blade in', num_energy_groups=2)
fuel2bin.setSigmaA([0.0081279, 0.08346091])
fuel2bin.setDifCoef([1.259, 0.2091])
fuel2bin.setNuSigmaF([0.004663, 0.1021])
fuel2bin.setChi([1.0, 0.0])
fuel2bin.setSigmaS(np.array([[0.0, 0.02617],[0.0, 0.0]]))

# create fuel 2 blade out
fuel2bo = rk.material.Material(name='fuel 2 blade out', num_energy_groups=2)
fuel2bo.setSigmaA([0.0081279, 0.07334491])
fuel2bo.setDifCoef([1.259, 0.2091])
fuel2bo.setNuSigmaF([0.004663, 0.1021])
fuel2bo.setChi([1.0, 0.0])
fuel2bo.setSigmaS(np.array([[0.0, 0.02617],[0.0, 0.0]]))

# create reflector
reflector = rk.material.Material(name='reflector', num_energy_groups=2)
reflector.setSigmaA([0.0007291, 0.01912592])
reflector.setDifCoef([1.257, 0.1592])
reflector.setNuSigmaF([0.0, 0.0])
reflector.setChi([1.0, 0.0])
reflector.setSigmaS(np.array([[0.0, 0.04754],[0.0, 0.0]]))

for j in xrange(9,11):
  for k in xrange(15):
    for cell in mesh._cells[j*15+k][:]: cell.setMaterial(reflector)

for j in xrange(9):
  for k in xrange(15):
    for cell in mesh._cells[j*15+k][7*15:]: cell.setMaterial(reflector)

for j in xrange(7,9):
  for k in xrange(15):
    for cell in mesh._cells[j*15+k][:7*15]: cell.setMaterial(fuel2bin)

for j in xrange(7):
  for k in xrange(15):
    for cell in mesh._cells[j*15+k][7*15:9*15]: cell.setMaterial(fuel2bin)

for k in xrange(15):
  for cell in mesh._cells[6*15+k][0*15:1*15]: cell.setMaterial(fuel1bo)
  for cell in mesh._cells[5*15+k][0*15:1*15]: cell.setMaterial(fuel1bo)

for j in xrange(5,7):
  for k in xrange(15):
    for cell in mesh._cells[j*15+k][5*15:7*15]: cell.setMaterial(fuel1bo)

for k in xrange(15):
  for cell in mesh._cells[0*15+k][0*15:1*15]: cell.setMaterial(fuel1bo)

for k in xrange(15):
  for cell in mesh._cells[0*15+k][5*15:7*15]: cell.setMaterial(fuel1bo)

for k in xrange(15):
  for cell in mesh._cells[7*15+k][7*15:8*15]: cell.setMaterial(fuel2bo)

for j in xrange(5,7):
  for k in xrange(15):
    for cell in mesh._cells[j*15+k][1*15:5*15]: cell.setMaterial(fuel1bin)

for j in xrange(1,5):
  for k in xrange(15):
    for cell in mesh._cells[j*15+k][:7*15]: cell.setMaterial(fuel1bin)

for k in xrange(15):
  for cell in mesh._cells[0*15+k][1*15:5*15]: cell.setMaterial(fuel1bin)

mesh.initializeSurfaces()

solver = rk.solver.Solver(mesh)
solver.setNumThreads(4)
solver.solve(1.e-6, 1000)

openrk.plotter.plot_flux(solver, [0,1], 500)
