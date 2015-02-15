from statepoint import StatePoint
import glob
from geometry import geometry
import infermc.plotter as plotter

# Get statepoint files
files = glob.glob('statepoint.*.h5')

for file in files:
  sp = StatePoint(file)
  sp.read_results()
  plotter.plot_fluxes(geometry, sp, energies=[0, 1], gridsize=200)
