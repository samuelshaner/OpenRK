# Arithmetic for arrays of uncorrelated variables stored in
# unumpy arrays with the Python "uncertainties" package:
# http://pythonhosted.org/uncertainties/
#
# NOTE: The mathematical routines are based on error propagation theory:
#        https://en.wikipedia.org/wiki/Propagation_of_uncertainty
#
# NOTE: These routines are used by the MultiGroupXS class to compute
#       group cross-sections from Monte Carlo tallies. The assumption that
#       the variables are uncorrelated is not correct - in reality
#       a covariance matrix should be used in a higher order term
#       for the standard deviation!!!

import uncertainties.unumpy as unumpy
import numpy as np


def unumpy_uncorr_add(a, b):

  # Extract the mean and standard deviation for a,b from unumpy arrays
  a_mean = unumpy.nominal_values(a)
  a_std_dev = unumpy.std_devs(a)
  b_mean = unumpy.nominal_values(b)
  b_std_dev = unumpy.std_devs(b)

  # Compute the mean and standard deviation of a+b
  c_mean = a_mean + b_mean
  c_std_dev = np.sqrt(a_std_dev**2 + b_std_dev**2)

  # Return a+b as a unumpy array
  return unumpy.uarray(c_mean, c_std_dev)


def unumpy_uncorr_sub(a, b):

  # Extract the mean and standard deviation for a,b from unumpy arrays
  a_mean = unumpy.nominal_values(a)
  a_std_dev = unumpy.std_devs(a)
  b_mean = unumpy.nominal_values(b)
  b_std_dev = unumpy.std_devs(b)

  # Compute the mean and standard deviation of a-b
  c_mean = a_mean - b_mean
  c_std_dev = np.sqrt(a_std_dev**2 + b_std_dev**2)

  # Return a-b as a unumpy array
  return unumpy.uarray(c_mean, c_std_dev)


def unumpy_uncorr_multiply(a, b):

  # Extract the mean and standard deviation for a,b from unumpy arrays
  a_mean = unumpy.nominal_values(a)
  a_std_dev = unumpy.std_devs(a)
  b_mean = unumpy.nominal_values(b)
  b_std_dev = unumpy.std_devs(b)

  # Compute the mean and standard deviation of a*b
  c_mean = a_mean * b_mean
  c_std_dev = c_mean**2 * ((a_std_dev/a_mean)**2 + (b_std_dev/b_mean)**2)

  # Return a*b as a unumpy array
  return unumpy.uarray(c_mean, c_std_dev)


def unumpy_uncorr_divide(a, b):

  # Extract the mean and standard deviation for a,b from unumpy arrays
  a_mean = unumpy.nominal_values(a)
  a_std_dev = unumpy.std_devs(a)
  b_mean = unumpy.nominal_values(b)
  b_std_dev = unumpy.std_devs(b)

  # Compute the mean and standard deviation of a/b
  c_mean = a_mean / b_mean
  c_std_dev = c_mean**2 * ((a_std_dev/a_mean)**2 + (b_std_dev/b_mean)**2)

  # Return a/b as a unumpy array
  return unumpy.uarray(c_mean, c_std_dev)
