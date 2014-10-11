# Arithmetic for arrays of correlated variables stored with both mean
# and standard deviation stored in NumPy arrays.
#
# NOTE: The mathematical routines are based on error propagation theory:
#        https://en.wikipedia.org/wiki/Propagation_of_uncertainty
#
# NOTE: These routines are used by the MultiGroupXS class to compute
#       group cross-sections from Monte Carlo tallies. This module assumes that
#       random variables may or may not be correlated and includes a second
#       term in the error propagation of the standard deviation.

import numpy as np
import itertools


def cov(a, b, micro=True):

  # Find the complete shape of the covariance array we must return
  # Assume that the only difference between a,b is in the final (nuclide) index
  if a.size == b.size:
    full_shape = b.shape
  elif a.size > b.size:
    full_shape = a.shape
    num_tiles = a.size / b.size
    b = np.repeat(b, num_tiles, axis=-1)
  else:
    full_shape = b.shape
    num_tiles = b.size / a.size
    a = np.repeat(a, num_tiles, axis=-1)

  # Allocate an empty NumPy array for the covariance
  covariance = np.zeros(full_shape)

  # If only one region, there is nothing to covary, so return 0
  num_domains = a.shape[0]
  if num_domains == 1:
    return covariance

  # Treat cross-sections for each group, nuclide as different random variables,
  # with different instances for different domains (e.g., distribcells)
  if micro:
    num_nuclides = a.shape[-1]
    num_groups = a.shape[-2]
    rand_var_shape = (num_groups, num_nuclides)
  else:
    num_groups = a.shape[-1]
    rand_var_shape = (num_groups, )

  # Create an iterator over each of the different cross-section random variables
  rand_var_map = tuple()
  for dim in rand_var_shape:
    rand_var_map += (range(dim),)

  rand_var_iterator = itertools.product(*rand_var_map)

  # Compute the covariance for each of the random variables
  for rand_var_index in rand_var_iterator:
    rv1 = a[..., rand_var_index].ravel()
    rv2 = b[..., rand_var_index].ravel()
    covariance[..., rand_var_index] = np.cov(rv1, rv2)[0][1]

  return covariance


def add(a, b, micro=True):

  # Extract the mean and standard deviation for a,b from NumPy arrays
  a_mean = a[0, ...]
  a_std_dev = a[1, ...]
  b_mean = b[0, ...]
  b_std_dev = b[1, ...]

  # Compute the mean and standard deviation of a+b
  c_mean = a_mean + b_mean
  c_std_dev = np.sqrt(a_std_dev**2 + b_std_dev**2 + \
                      2.*cov(a_mean, b_mean, micro))

  # Return a+b as a unumpy array
  return np.array([c_mean, c_std_dev])


def sub(a, b, micro=True):

  # Extract the mean and standard deviation for a,b from NumPy arrays
  a_mean = a[0, ...]
  a_std_dev = a[1, ...]
  b_mean = b[0, ...]
  b_std_dev = b[1, ...]

  # Compute the mean and standard deviation of a-b
  c_mean = a_mean - b_mean
  c_std_dev = np.sqrt(a_std_dev**2 + b_std_dev**2 + \
                      2.*cov(a_mean, b_mean, micro))

  # Return a-b as a unumpy array
  return np.array([c_mean, c_std_dev])


def multiply_by_array(a, b, micro=True):

  # Extract the mean and standard deviation for a,b from NumPy arrays
  a_mean = a[0, ...]
  a_std_dev = a[1, ...]
  b_mean = b[0, ...]
  b_std_dev = b[1, ...]

  # Compute the mean and standard deviation of a*b
  c_mean = a_mean * b_mean
  c_std_dev = np.abs(c_mean) * \
              np.sqrt((a_std_dev/a_mean)**2+(b_std_dev/b_mean)**2 + \
                      (2.*cov(a_mean, b_mean, micro) / c_mean))

  # Return a*b as a unumpy array
  return np.array([c_mean, c_std_dev])


def multiply_by_scalar(a, b, micro=True):

  # Extract the mean and standard deviation for a from NumPy array
  a_mean = a[0, ...]
  a_std_dev = a[1, ...]
  b_mean = b[0, ...]
  b_std_dev = b[1, ...]

  # Compute the mean and standard deviation of a*b
  c_mean = a_mean * b
  c_std_dev = a_std_dev * b

  # Return a*b as a unumpy array
  return np.array([c_mean, c_std_dev])


def divide_by_array(a, b, micro=True):

  # Extract the mean and standard deviation for a,b from NumPy arrays
  a_mean = a[0, ...]
  a_std_dev = a[1, ...]
  b_mean = b[0, ...]
  b_std_dev = b[1, ...]

  # Compute the mean and standard deviation of a/b
  c_mean = a_mean / b_mean
  c_std_dev = np.abs(c_mean) * \
              np.sqrt((a_std_dev/a_mean)**2+(b_std_dev/b_mean)**2 - \
                      (2.*cov(a_mean, b_mean, micro) / (a_mean * b_mean)))

  # Return a/b as a unumpy array
  return np.array([c_mean, c_std_dev])


def divide_by_scalar(a, b, micro=True):

  # Extract the mean and standard deviation for a from NumPy array
  a_mean = a[0, ...]
  a_std_dev = a[1, ...]

  # Compute the mean and standard deviation of a/b
  c_mean = a_mean / b
  c_std_dev = (1./b) * a_std_dev

  # Return a/b as a unumpy array
  return np.array([c_mean, c_std_dev])