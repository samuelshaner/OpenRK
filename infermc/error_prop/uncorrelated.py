# Arithmetic for arrays of uncorrelated variables stored with both mean
# and standard deviation stored in NumPy arrays.
#
# NOTE: The mathematical routines are based on error propagation theory:
#        https://en.wikipedia.org/wiki/Propagation_of_uncertainty
#
# NOTE: These routines are used by the MultiGroupXS class to compute
#       group cross-sections from Monte Carlo tallies. The assumption that
#       the variables are uncorrelated is not correct - in reality
#       a covariance matrix should be used in a higher order term
#       for the standard deviation!!!

import numpy as np


def add(a, b):

  # Extract the mean and standard deviation for a,b from NumPy arrays
  a_mean = a[0, ...]
  a_std_dev = a[1, ...]
  b_mean = b[0, ...]
  b_std_dev = b[1, ...]

  # Compute the mean and standard deviation of a+b
  c_mean = a_mean + b_mean
  c_std_dev = np.sqrt(a_std_dev**2 + b_std_dev**2)

  # Return a+b as a unumpy array
  return np.array([c_mean, c_std_dev])


def sub(a, b):

  # Extract the mean and standard deviation for a,b from NumPy arrays
  a_mean = a[0, ...]
  a_std_dev = a[1, ...]
  b_mean = b[0, ...]
  b_std_dev = b[1, ...]

  # Compute the mean and standard deviation of a-b
  c_mean = a_mean - b_mean
  c_std_dev = np.sqrt(a_std_dev**2 + b_std_dev**2)

  # Return a-b as a unumpy array
  return np.array([c_mean, c_std_dev])


def multiply_by_array(a, b):

  # Extract the mean and standard deviation for a,b from NumPy arrays
  a_mean = a[0, ...]
  a_std_dev = a[1, ...]
  b_mean = b[0, ...]
  b_std_dev = b[1, ...]

  # Compute the mean and standard deviation of a*b
  c_mean = a_mean * b_mean
  c_std_dev = np.abs(c_mean)*np.sqrt((a_std_dev/a_mean)**2+(b_std_dev/b_mean)**2)

  # Return a*b as a unumpy array
  return np.array([c_mean, c_std_dev])


def multiply_by_scalar(a, b):

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


def divide_by_array(a, b):

  # Extract the mean and standard deviation for a,b from NumPy arrays
  a_mean = a[0, ...]
  a_std_dev = a[1, ...]
  b_mean = b[0, ...]
  b_std_dev = b[1, ...]

  # Compute the mean and standard deviation of a/b
  c_mean = a_mean / b_mean
  c_std_dev = np.abs(c_mean)*np.sqrt((a_std_dev/a_mean)**2+(b_std_dev/b_mean)**2)

  # Return a/b as a unumpy array
  return np.array([c_mean, c_std_dev])


def divide_by_scalar(a, b):

  # Extract the mean and standard deviation for a from NumPy array
  a_mean = a[0, ...]
  a_std_dev = a[1, ...]

  # Compute the mean and standard deviation of a/b
  c_mean = a_mean / b
  c_std_dev = (1./b) * a_std_dev

  # Return a/b as a unumpy array
  return np.array([c_mean, c_std_dev])