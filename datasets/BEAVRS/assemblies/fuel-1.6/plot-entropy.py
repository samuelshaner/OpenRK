import matplotlib.pyplot as plt
import numpy as np
import h5py

num_batches = 25

file = h5py.File('statepoint.' + str(num_batches) + '.h5')
entropy = file['entropy'][...]
batches = np.arange(num_batches)

fig = plt.figure()
plt.plot(batches, entropy)
plt.title('Shannon Entropy - 1.6% Enr. Assembly')
plt.xlabel('Batch #')
plt.ylabel('Entropy')
plt.show()

file.close()