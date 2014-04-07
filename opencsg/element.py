

class Element(object):

    def __init__(self, name='', xs=None):

        # Initialize class attributes
        self._name = ''
        self._xs = None

        # Set the Material class attributes
        self.setName(name)

        if not xs is None:
            self.setXS(xs)


    def getName(self):
        return self._name


    def getXS(self):
        return self._xs


    def setName(self, name):

        if not isinstance(name, str):
            exit('Unable to set name for Element with a non-string value %s',
                 str(name))

        else:
            self._name = name


    def setXS(self, xs):

        if not isinstance(xs, str):
            exit('Unable to set cross-section identifier xs for Element with '
                 'a non-string value %s', str(xs))

        else:
            self._xs = xs


    def toString(self):

        string = 'Element  -  ' + self._name

        xs = '{0: <16}'.format('\tXS') + '=\t' + self._xs
        string += xs + '\n'

        return string


    def printString(self):
        print(self.toString())
