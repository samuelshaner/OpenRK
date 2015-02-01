__author__ = 'Samuel Shaner'
__email__ = 'shaner@mit.edu'

# Import modules
import checkvalue as cv


class Clock(object):
    def __init__(self, start_time=0.0, end_time=3.0, dt_outer=1e-2, dt_inner=1e-3):

        self._positions = ['START', 'PREVIOUS_OUT', 'PREVIOUS_IN', 'CURRENT', 'FORWARD_OUT', 'END']

        # Initialize class attributes
        self._times = dict()
        self._times['START'] = start_time
        self._times['PREVIOUS_OUT'] = 0.0
        self._times['PREVIOUS_IN'] = 0.0
        self._times['CURRENT'] = 0.0
        self._times['FORWARD_OUT'] = 0.0
        self._times['END'] = end_time
        self._dt_outer = None
        self._dt_inner = None

        self.set_dt_outer(dt_outer)
        self.set_dt_inner(dt_inner)

    def get_positions(self):

        return self._positions

    def set_dt_outer(self, dt_outer):

        if not cv.is_float(dt_outer):
            msg = 'Unable to set DT outer for non-float value {1}' \
                .format(dt_outer)
            raise ValueError(msg)

        elif dt_outer <= 0.0:
            msg = 'Unable to set DT outer for non positive value {1}' \
                .format(dt_outer)

            print msg

        else:
            self._dt_outer = dt_outer

    def set_dt_inner(self, dt_inner):

        if not cv.is_float(dt_inner):
            msg = 'Unable to set DT inner for non-float value {1}' \
                .format(dt_inner)
            raise ValueError(msg)

        elif dt_inner <= 0.0:
            msg = 'Unable to set DT inner for non positive value {1}' \
                .format(dt_inner)
            raise ValueError(msg)

        else:
            self._dt_inner = dt_inner

    def take_inner_step(self):

        # Check to make sure inner step is multiple of outer time step size
        if self._dt_outer % self._dt_inner > 1e-8:
            msg = 'Unable to take inner step since DT outer is not an integer ' \
                  'multiple of DT inner. DT inner: {0}, DT outer: {1}' \
                .format(self._dt_inner, self._dt_outer)
            print msg

            # Check to make sure that CURRENT is less than FORWARD_OUT
        elif self._times['CURRENT'] >= self._times['FORWARD_OUT']:
            msg = 'Unable to take inner step since CURRENT time is not ' \
                  'less than the FORWARD_OUT time. CURRENT: {0}, FORWARD_OUT: {1}' \
                .format(self._times['CURRENT'], self._times['FORWARD_OUT'])
            print msg

        else:

            self._times['PREVIOUS_IN'] = self._times['CURRENT']
            self._times['CURRENT'] = self._times['CURRENT'] + self._dt_inner

    def take_outer_step(self):

        # Check to make sure inner step is multiple of outer time step size
        if self._dt_outer % self._dt_inner > 1e-8:
            msg = 'Unable to take outer step since DT outer is not an integer ' \
                  'multiple of DT inner. DT inner: {0}, DT outer: {1}' \
                .format(self._dt_inner, self._dt_outer)
            print msg

            # Check to make sure that CURRENT time equals FORWARD_OUT
        elif self._times['CURRENT'] != self._times['FORWARD_OUT']:
            msg = 'Unable to take outer step since CURRENT time is not equal to ' \
                  'FORWARD_OUT time. CURRENT: {0}, FORWARD_OUT: {1}' \
                .format(self._times['CURRENT'], self._times['FORWARD_OUT'])
            print msg

        else:

            self._times['PREVIOUS_OUT'] = self._times['FORWARD_OUT']
            self._times['FORWARD_OUT'] = self._times['FORWARD_OUT'] + self._dt_outer

            if self._times['FORWARD_OUT'] > self._times['END']:
                self._times['FORWARD_OUT'] = self._times['END']

            # set CURRENT and PREVIOUS_IN to PREVIOUS_OUT
            self._times['CURRENT'] = self._times['PREVIOUS_OUT']
            self._times['PREVIOUS_IN'] = self._times['PREVIOUS_OUT']

    def reset_to_previous_outer_step(self):

        # Check to make sure that CURRENT time equals FORWARD_OUT
        if self._times['CURRENT'] != self._times['FORWARD_OUT']:
            msg = 'Unable to reset to previous out since CURRENT time is not equal to ' \
                  'FORWARD_OUT time. CURRENT: {0}, FORWARD_OUT: {1}' \
                .format(self._times['CURRENT'], self._times['FORWARD_OUT'])
            print msg

        else:

            self._times['CURRENT'] = self._times['PREVIOUS_OUT']
            self._times['PREVIOUS_IN'] = self._times['PREVIOUS_OUT']

    def __repr__(self):

        string = 'OpenRK Clock\n'
        string += ' START time \t\t\t = {:6.5f} \n'.format(self._times['START'])
        string += ' PREVIOUS_OUT time \t\t = {:6.5f} \n'.format(self._times['PREVIOUS_OUT'])
        string += ' PREVIOUS_IN time \t\t = {:6.5f} \n'.format(self._times['PREVIOUS_IN'])
        string += ' CURRENT time \t\t\t = {:6.5f} \n'.format(self._times['CURRENT'])
        string += ' FORWARD_OUT time \t\t = {:6.5f} \n'.format(self._times['FORWARD_OUT'])
        string += ' END time \t\t\t = {:6.5f} \n'.format(self._times['END'])
        string += ' DT outer \t\t\t = {:6.5f} \n'.format(self._dt_outer)
        string += ' DT inner \t\t\t = {:6.5f} \n'.format(self._dt_inner)

        return string
