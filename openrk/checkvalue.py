__author__ = 'Samuel Shaner'
__email__ = 'shaner@mit.edu'

# Import modules
import numpy as np
import clock as ck


def is_integer(val):
    return isinstance(val, (int, np.int32, np.int64))


def is_float(val):
    return isinstance(val, (float, np.float32, np.float64))


def is_string(val):
    return isinstance(val, (str, np.str))


def is_list(val):
    return isinstance(val, (list, np.ndarray))


def is_bool(val):
    return isinstance(val, bool)


def check_set(val, set_name, val_name):
    if val is None:
        msg = 'Unable to get/set {0} since {1} has not been set.'. \
            format(set_name, val_name)
        raise ValueError(msg)


def check_list_of_floats(val, name):
    check_is_list(val, name)

    for i in val:
        if not is_float(i):
            msg = 'Unable to get/set {0} since {1} it is not a list of all floats.'. \
                format(name, val)
            raise ValueError(msg)


def check_list_of_ints(val, name):
    check_is_list(val, name)

    for i in val:
        if not is_integer(i):
            msg = 'Unable to get/set {0} since {1} it is not a list of all ints.'. \
                format(name, val)
            raise ValueError(msg)


def check_list_of_floats_or_ints(val, name):
    check_is_list(val, name)

    for i in val:
        if not is_integer(i) and not is_float(i):
            msg = 'Unable to get/set {0} since {1} it is not a list of floats or ints.'. \
                format(name, val)
            raise ValueError(msg)


def check_is_list(val, name):
    if not is_list(val):
        msg = 'Unable to set {0} with a ' \
              'non-list value {1}'.format(name, val)
        raise ValueError(msg)


def check_is_float(val, name, param):
    if not is_float(val):
        msg = 'Unable to get/set {0} with a non-float value {1} for' \
              ' parameter {2}'.format(name, val, param)
        raise ValueError(msg)


def check_is_int(val, name, param):
    if not is_integer(val):
        msg = 'Unable to get/set {0} with a non-int value {1} for' \
              ' parameter {2}'.format(name, val, param)
        raise ValueError(msg)


def check_is_float_or_int(val, name, param):
    if not is_integer(val) and not is_float(val):
        msg = 'Unable to get/set {0} with a non-int and not-float value {1} for' \
              ' parameter {2}'.format(name, val, param)
        raise ValueError(msg)

def check_is_bool(val, name, param):
    if not is_bool(val):
        msg = 'Unable to get/set {0} with a non-bool value {1} for' \
              ' parameter {2}'.format(name, val, param)
        raise ValueError(msg)
