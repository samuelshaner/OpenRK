import openrk as rk
mesh = rk.StructuredMesh(1.0, 1.0, 1, 1)
mat1 = rk.Material(name='mat1')
mesh._cells[0][0].setMaterial(mat1)
mesh.initializeSurfaces()
solver = rk.Solver(mesh)
mat1.setSigmaA([1.0])
mat1.setSigmaF([1.0])
mat1.setNuSigmaF([1.0])
mat1.setDifCoef([1.0])
mat1.setChi([1.0])
mat1.setSigmaS([[1.0]])
solver.solve(1.e-6, 1000)


