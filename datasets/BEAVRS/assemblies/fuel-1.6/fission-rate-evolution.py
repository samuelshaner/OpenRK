import re
import numpy, h5py
import matplotlib.pyplot as plt


def get_pin_powers(time_key):
  '''Get the pin power array from an OpenMOC simulation state file'''

  # The time stamp in the simulation state file
  time_group = date_group[time_key]

  # Find the OpenMC batch number
  note = time_group.attrs['note']
  batch = int(note.split('-')[1])
  print('Batch: {0}'.format(batch))

  # Allocate memory for a 2D array of pin powers
  pin_powers = numpy.zeros((17,17), dtype=numpy.float64)

  # Retrieve the fission rates dictionary
  # Keys    - FSR "paths"
  # Values  - Fission rates
  fiss_rates = time_group['fission-rates']['fission-rates']

  # Loop over all keys (FSR paths)
  for rate in fiss_rates.attrs:

    # Find the x,y indices into the lattice from the FSR path
    index = re.findall('[0-9]*, [0-9]*', rate)
    index = index[0].split(',')
    lat_x = int(index[0])
    lat_y = int(index[1])

    pin_powers[lat_x,lat_y] = float(fiss_rates.attrs[rate])

  # Normalize the pin powers to the maximum
  pin_powers /= numpy.max(pin_powers)

  # Return both the OpenMOC pin powers and OpenMC at this timestamp
  return pin_powers, batch


# Get reference pin powers
f = h5py.File('statepoint.100.h5', 'r')
ref_pin_powers = f['tallies']['tally 10000']['results']['sum'][...]
ref_pin_powers.shape = (17,17)
ref_pin_powers /= numpy.max(ref_pin_powers)
f.close()

# Initialize arrays for the max, mean errors for each batch
max_err = numpy.zeros(91)
mean_err = numpy.zeros(91)

# Open the OpenMOC simulation state file
f = h5py.File('sim-state-2-group.h5', 'r')
date_key = f.keys()[0]
date_group = f[date_key]
time_keys = date_group.keys()

# Loop over batches
for i, time_key in enumerate(time_keys):


  # Calculate the pin power percent relative error
  pin_powers, batch = get_pin_powers(time_key)
  pin_power_err = (pin_powers - ref_pin_powers) / ref_pin_powers
  pin_power_err *= 100.
  pin_power_err = numpy.nan_to_num(pin_power_err)
  pin_power_err = numpy.fabs(pin_power_err)

  # Find the max and mean percent relative error
  max_err[i] = numpy.max(pin_power_err)
  mean_err[i] = numpy.mean(pin_power_err)

  # Plot the 
  fig = plt.figure("Pin Powers")
  plot = plt.imshow(pin_power_err, interpolation='nearest')
  plt.xlabel('X')
  plt.ylabel('Y')
  plt.title('Batch {0} Pin Powers'.format(batch))
  plt.colorbar()
  plt.savefig('powers/pin-powers-{0}.png'.format(i), bbox_inches='tight')
  plt.close()

# Plot the pin power error convergence rate
fig = plt.figure()
plt.plot(numpy.arange(10, 101), max_err, linewidth=2)
plt.plot(numpy.arange(10, 101), mean_err, linewidth=2)
plt.title('Pin Power Error')
plt.xlabel('Batch')
plt.ylabel('% Error')
plt.legend(['max', 'mean'])
plt.grid()
plt.savefig('powers/error-convergence.png', bbox_inches='tight')
plt.close()
