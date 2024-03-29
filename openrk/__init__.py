import sys, os
import random
import datetime
import signal

# For Python 2.X.X
if (sys.version_info[0] == 2):
  from openrk import *
# For Python 3.X.X
else:
  from openrk.openrk import *

# Tell Python to recognize CTRL+C and stop the C++ extension module
# when this is passed in from the keyboard
signal.signal(signal.SIGINT, signal.SIG_DFL)

from clock import *
from checkvalue import *
from material import *
from mesh import *
from checkvalue import *
from cell import *
from solver import *
from plotter import *
from transient import *