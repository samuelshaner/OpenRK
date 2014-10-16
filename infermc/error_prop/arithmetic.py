import correlated, uncorrelated


def add(a, b, corr=False, micro=True):
  if corr:
    return correlated.add(a, b, micro=micro)
  else:
    return uncorrelated.add(a, b)


def sub(a, b, corr=False, micro=True):
  if corr:
    return correlated.sub(a, b, micro=micro)
  else:
    return uncorrelated.sub(a, b)


def multiply_by_array(a, b, corr=False, micro=True):
  if corr:
    return correlated.multiply_by_array(a, b, micro=micro)
  else:
    return uncorrelated.multiply_by_array(a, b)


def multiply_by_scalar(a, b, corr=False, micro=True):
  if corr:
    return correlated.multiply_by_scalar(a, b, micro=micro)
  else:
    return uncorrelated.multiply_by_scalar(a, b)


def divide_by_array(a, b, corr=False, micro=True):
  if corr:
    return correlated.divide_by_array(a, b, micro=micro)
  else:
    return uncorrelated.divide_by_array(a, b)


def divide_by_scalar(a, b, corr=False, micro=True):
  if corr:
    return correlated.divide_by_scalar(a, b, micro=micro)
  else:
    return uncorrelated.divide_by_scalar(a, b)