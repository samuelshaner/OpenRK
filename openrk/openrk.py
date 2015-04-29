# This file was automatically generated by SWIG (http://www.swig.org).
# Version 2.0.11
#
# Do not make changes to this file unless you know what you are doing--modify
# the SWIG interface file instead.





from sys import version_info
if version_info >= (2,6,0):
    def swig_import_helper():
        from os.path import dirname
        import imp
        fp = None
        try:
            fp, pathname, description = imp.find_module('_openrk', [dirname(__file__)])
        except ImportError:
            import _openrk
            return _openrk
        if fp is not None:
            try:
                _mod = imp.load_module('_openrk', fp, pathname, description)
            finally:
                fp.close()
            return _mod
    _openrk = swig_import_helper()
    del swig_import_helper
else:
    import _openrk
del version_info
try:
    _swig_property = property
except NameError:
    pass # Python < 2.2 doesn't have 'property'.
def _swig_setattr_nondynamic(self,class_type,name,value,static=1):
    if (name == "thisown"): return self.this.own(value)
    if (name == "this"):
        if type(value).__name__ == 'SwigPyObject':
            self.__dict__[name] = value
            return
    method = class_type.__swig_setmethods__.get(name,None)
    if method: return method(self,value)
    if (not static):
        self.__dict__[name] = value
    else:
        raise AttributeError("You cannot add attributes to %s" % self)

def _swig_setattr(self,class_type,name,value):
    return _swig_setattr_nondynamic(self,class_type,name,value,0)

def _swig_getattr(self,class_type,name):
    if (name == "thisown"): return self.this.own()
    method = class_type.__swig_getmethods__.get(name,None)
    if method: return method(self)
    raise AttributeError(name)

def _swig_repr(self):
    try: strthis = "proxy of " + self.this.__repr__()
    except: strthis = ""
    return "<%s.%s; %s >" % (self.__class__.__module__, self.__class__.__name__, strthis,)

try:
    _object = object
    _newclass = 1
except AttributeError:
    class _object : pass
    _newclass = 0



def eigenvalueSolve(*args, **kwargs):
  return _openrk.eigenvalueSolve(*args, **kwargs)
eigenvalueSolve = _openrk.eigenvalueSolve

def eigenvalueSolve2d(*args, **kwargs):
  return _openrk.eigenvalueSolve2d(*args, **kwargs)
eigenvalueSolve2d = _openrk.eigenvalueSolve2d

def linearSolve(*args, **kwargs):
  return _openrk.linearSolve(*args, **kwargs)
linearSolve = _openrk.linearSolve

def linearSolve2d(*args, **kwargs):
  return _openrk.linearSolve2d(*args, **kwargs)
linearSolve2d = _openrk.linearSolve2d

def matrix_multiplication(*args, **kwargs):
  return _openrk.matrix_multiplication(*args, **kwargs)
matrix_multiplication = _openrk.matrix_multiplication

def matrix_multiplication2d(*args, **kwargs):
  return _openrk.matrix_multiplication2d(*args, **kwargs)
matrix_multiplication2d = _openrk.matrix_multiplication2d

def vector_scale(*args, **kwargs):
  return _openrk.vector_scale(*args, **kwargs)
vector_scale = _openrk.vector_scale

def setNumThreads(*args, **kwargs):
  return _openrk.setNumThreads(*args, **kwargs)
setNumThreads = _openrk.setNumThreads

def matMultA(*args, **kwargs):
  return _openrk.matMultA(*args, **kwargs)
matMultA = _openrk.matMultA
DEBUG = _openrk.DEBUG
INFO = _openrk.INFO
NORMAL = _openrk.NORMAL
SEPARATOR = _openrk.SEPARATOR
HEADER = _openrk.HEADER
TITLE = _openrk.TITLE
WARNING = _openrk.WARNING
CRITICAL = _openrk.CRITICAL
RESULT = _openrk.RESULT
UNITTEST = _openrk.UNITTEST
ERROR = _openrk.ERROR

def set_err(*args, **kwargs):
  return _openrk.set_err(*args, **kwargs)
set_err = _openrk.set_err

def set_output_directory(*args, **kwargs):
  return _openrk.set_output_directory(*args, **kwargs)
set_output_directory = _openrk.set_output_directory

def get_output_directory():
  return _openrk.get_output_directory()
get_output_directory = _openrk.get_output_directory

def set_log_filename(*args, **kwargs):
  return _openrk.set_log_filename(*args, **kwargs)
set_log_filename = _openrk.set_log_filename

def get_log_filename():
  return _openrk.get_log_filename()
get_log_filename = _openrk.get_log_filename

def set_separator_character(*args, **kwargs):
  return _openrk.set_separator_character(*args, **kwargs)
set_separator_character = _openrk.set_separator_character

def get_separator_character():
  return _openrk.get_separator_character()
get_separator_character = _openrk.get_separator_character

def set_header_character(*args, **kwargs):
  return _openrk.set_header_character(*args, **kwargs)
set_header_character = _openrk.set_header_character

def get_header_character():
  return _openrk.get_header_character()
get_header_character = _openrk.get_header_character

def set_title_character(*args, **kwargs):
  return _openrk.set_title_character(*args, **kwargs)
set_title_character = _openrk.set_title_character

def get_title_character():
  return _openrk.get_title_character()
get_title_character = _openrk.get_title_character

def set_line_length(*args, **kwargs):
  return _openrk.set_line_length(*args, **kwargs)
set_line_length = _openrk.set_line_length

def set_log_level(*args, **kwargs):
  return _openrk.set_log_level(*args, **kwargs)
set_log_level = _openrk.set_log_level

def get_log_level():
  return _openrk.get_log_level()
get_log_level = _openrk.get_log_level

def log_printf(*args, **kwargs):
  return _openrk.log_printf(*args, **kwargs)
log_printf = _openrk.log_printf

def create_multiline_msg(*args, **kwargs):
  return _openrk.create_multiline_msg(*args, **kwargs)
create_multiline_msg = _openrk.create_multiline_msg

def material_id():
  return _openrk.material_id()
material_id = _openrk.material_id
class Material(_object):
    __swig_setmethods__ = {}
    __setattr__ = lambda self, name, value: _swig_setattr(self, Material, name, value)
    __swig_getmethods__ = {}
    __getattr__ = lambda self, name: _swig_getattr(self, Material, name)
    __repr__ = _swig_repr
    def __init__(self, id=0): 
        this = _openrk.new_Material(id)
        try: self.this.append(this)
        except: self.this = this
    __swig_destroy__ = _openrk.delete_Material
    __del__ = lambda self : None;
    def setNumEnergyGroups(self, *args, **kwargs): return _openrk.Material_setNumEnergyGroups(self, *args, **kwargs)
    def setNumDelayedGroups(self, *args, **kwargs): return _openrk.Material_setNumDelayedGroups(self, *args, **kwargs)
    def setEnergyPerFission(self, *args, **kwargs): return _openrk.Material_setEnergyPerFission(self, *args, **kwargs)
    def setSigmaT(self, *args, **kwargs): return _openrk.Material_setSigmaT(self, *args, **kwargs)
    def setSigmaA(self, *args, **kwargs): return _openrk.Material_setSigmaA(self, *args, **kwargs)
    def setSigmaS(self, *args, **kwargs): return _openrk.Material_setSigmaS(self, *args, **kwargs)
    def setSigmaF(self, *args, **kwargs): return _openrk.Material_setSigmaF(self, *args, **kwargs)
    def setNuSigmaF(self, *args, **kwargs): return _openrk.Material_setNuSigmaF(self, *args, **kwargs)
    def setChi(self, *args, **kwargs): return _openrk.Material_setChi(self, *args, **kwargs)
    def setDifCoef(self, *args, **kwargs): return _openrk.Material_setDifCoef(self, *args, **kwargs)
    def setVelocity(self, *args, **kwargs): return _openrk.Material_setVelocity(self, *args, **kwargs)
    def setPrecursorConc(self, *args, **kwargs): return _openrk.Material_setPrecursorConc(self, *args, **kwargs)
    def setSigmaTByGroup(self, *args, **kwargs): return _openrk.Material_setSigmaTByGroup(self, *args, **kwargs)
    def setSigmaAByGroup(self, *args, **kwargs): return _openrk.Material_setSigmaAByGroup(self, *args, **kwargs)
    def setSigmaFByGroup(self, *args, **kwargs): return _openrk.Material_setSigmaFByGroup(self, *args, **kwargs)
    def setNuSigmaFByGroup(self, *args, **kwargs): return _openrk.Material_setNuSigmaFByGroup(self, *args, **kwargs)
    def setSigmaSByGroup(self, *args, **kwargs): return _openrk.Material_setSigmaSByGroup(self, *args, **kwargs)
    def setChiByGroup(self, *args, **kwargs): return _openrk.Material_setChiByGroup(self, *args, **kwargs)
    def setDifCoefByGroup(self, *args, **kwargs): return _openrk.Material_setDifCoefByGroup(self, *args, **kwargs)
    def setVelocityByGroup(self, *args, **kwargs): return _openrk.Material_setVelocityByGroup(self, *args, **kwargs)
    def setPrecursorConcByGroup(self, *args, **kwargs): return _openrk.Material_setPrecursorConcByGroup(self, *args, **kwargs)
    def setTemperatureConversionFactor(self, *args, **kwargs): return _openrk.Material_setTemperatureConversionFactor(self, *args, **kwargs)
    def setClock(self, *args, **kwargs): return _openrk.Material_setClock(self, *args, **kwargs)
    def getId(self): return _openrk.Material_getId(self)
    def getNumEnergyGroups(self): return _openrk.Material_getNumEnergyGroups(self)
    def getNumDelayedGroups(self): return _openrk.Material_getNumDelayedGroups(self)
    def getEnergyPerFission(self): return _openrk.Material_getEnergyPerFission(self)
    def getSigmaT(self): return _openrk.Material_getSigmaT(self)
    def getSigmaA(self): return _openrk.Material_getSigmaA(self)
    def getSigmaS(self): return _openrk.Material_getSigmaS(self)
    def getSigmaF(self): return _openrk.Material_getSigmaF(self)
    def getNuSigmaF(self): return _openrk.Material_getNuSigmaF(self)
    def getChi(self): return _openrk.Material_getChi(self)
    def getDifCoef(self): return _openrk.Material_getDifCoef(self)
    def getVelocity(self): return _openrk.Material_getVelocity(self)
    def getPrecursorConc(self): return _openrk.Material_getPrecursorConc(self)
    def getSigmaTByGroup(self, *args, **kwargs): return _openrk.Material_getSigmaTByGroup(self, *args, **kwargs)
    def getSigmaAByGroup(self, *args, **kwargs): return _openrk.Material_getSigmaAByGroup(self, *args, **kwargs)
    def getSigmaSByGroup(self, *args, **kwargs): return _openrk.Material_getSigmaSByGroup(self, *args, **kwargs)
    def getSigmaFByGroup(self, *args, **kwargs): return _openrk.Material_getSigmaFByGroup(self, *args, **kwargs)
    def getNuSigmaFByGroup(self, *args, **kwargs): return _openrk.Material_getNuSigmaFByGroup(self, *args, **kwargs)
    def getChiByGroup(self, *args, **kwargs): return _openrk.Material_getChiByGroup(self, *args, **kwargs)
    def getDifCoefByGroup(self, *args, **kwargs): return _openrk.Material_getDifCoefByGroup(self, *args, **kwargs)
    def getVelocityByGroup(self, *args, **kwargs): return _openrk.Material_getVelocityByGroup(self, *args, **kwargs)
    def getPrecursorConcByGroup(self, *args, **kwargs): return _openrk.Material_getPrecursorConcByGroup(self, *args, **kwargs)
    def getTemperatureConversionFactor(self): return _openrk.Material_getTemperatureConversionFactor(self)
    def isFissionable(self): return _openrk.Material_isFissionable(self)
    def toString(self): return _openrk.Material_toString(self)
    def printString(self): return _openrk.Material_printString(self)
    def clone(self): return _openrk.Material_clone(self)
    def copy(self, *args, **kwargs): return _openrk.Material_copy(self, *args, **kwargs)
Material_swigregister = _openrk.Material_swigregister
Material_swigregister(Material)

class FunctionalMaterial(Material):
    __swig_setmethods__ = {}
    for _s in [Material]: __swig_setmethods__.update(getattr(_s,'__swig_setmethods__',{}))
    __setattr__ = lambda self, name, value: _swig_setattr(self, FunctionalMaterial, name, value)
    __swig_getmethods__ = {}
    for _s in [Material]: __swig_getmethods__.update(getattr(_s,'__swig_getmethods__',{}))
    __getattr__ = lambda self, name: _swig_getattr(self, FunctionalMaterial, name)
    __repr__ = _swig_repr
    def __init__(self, id=0): 
        this = _openrk.new_FunctionalMaterial(id)
        try: self.this.append(this)
        except: self.this = this
    __swig_destroy__ = _openrk.delete_FunctionalMaterial
    __del__ = lambda self : None;
    def setNumEnergyGroups(self, *args, **kwargs): return _openrk.FunctionalMaterial_setNumEnergyGroups(self, *args, **kwargs)
    def setTimeSteps(self, *args, **kwargs): return _openrk.FunctionalMaterial_setTimeSteps(self, *args, **kwargs)
    def setSigmaT(self, *args, **kwargs): return _openrk.FunctionalMaterial_setSigmaT(self, *args, **kwargs)
    def setSigmaA(self, *args, **kwargs): return _openrk.FunctionalMaterial_setSigmaA(self, *args, **kwargs)
    def setSigmaS(self, *args, **kwargs): return _openrk.FunctionalMaterial_setSigmaS(self, *args, **kwargs)
    def setSigmaF(self, *args, **kwargs): return _openrk.FunctionalMaterial_setSigmaF(self, *args, **kwargs)
    def setNuSigmaF(self, *args, **kwargs): return _openrk.FunctionalMaterial_setNuSigmaF(self, *args, **kwargs)
    def setChi(self, *args, **kwargs): return _openrk.FunctionalMaterial_setChi(self, *args, **kwargs)
    def setDifCoef(self, *args, **kwargs): return _openrk.FunctionalMaterial_setDifCoef(self, *args, **kwargs)
    def setVelocity(self, *args, **kwargs): return _openrk.FunctionalMaterial_setVelocity(self, *args, **kwargs)
    def setDopplerCoefficients(self, *args, **kwargs): return _openrk.FunctionalMaterial_setDopplerCoefficients(self, *args, **kwargs)
    def setSigmaTByGroup(self, *args, **kwargs): return _openrk.FunctionalMaterial_setSigmaTByGroup(self, *args, **kwargs)
    def setSigmaAByGroup(self, *args, **kwargs): return _openrk.FunctionalMaterial_setSigmaAByGroup(self, *args, **kwargs)
    def setSigmaFByGroup(self, *args, **kwargs): return _openrk.FunctionalMaterial_setSigmaFByGroup(self, *args, **kwargs)
    def setNuSigmaFByGroup(self, *args, **kwargs): return _openrk.FunctionalMaterial_setNuSigmaFByGroup(self, *args, **kwargs)
    def setSigmaSByGroup(self, *args, **kwargs): return _openrk.FunctionalMaterial_setSigmaSByGroup(self, *args, **kwargs)
    def setChiByGroup(self, *args, **kwargs): return _openrk.FunctionalMaterial_setChiByGroup(self, *args, **kwargs)
    def setDifCoefByGroup(self, *args, **kwargs): return _openrk.FunctionalMaterial_setDifCoefByGroup(self, *args, **kwargs)
    def setVelocityByGroup(self, *args, **kwargs): return _openrk.FunctionalMaterial_setVelocityByGroup(self, *args, **kwargs)
    def getSigmaTByGroup(self, *args, **kwargs): return _openrk.FunctionalMaterial_getSigmaTByGroup(self, *args, **kwargs)
    def getSigmaAByGroup(self, *args, **kwargs): return _openrk.FunctionalMaterial_getSigmaAByGroup(self, *args, **kwargs)
    def getSigmaSByGroup(self, *args, **kwargs): return _openrk.FunctionalMaterial_getSigmaSByGroup(self, *args, **kwargs)
    def getSigmaFByGroup(self, *args, **kwargs): return _openrk.FunctionalMaterial_getSigmaFByGroup(self, *args, **kwargs)
    def getNuSigmaFByGroup(self, *args, **kwargs): return _openrk.FunctionalMaterial_getNuSigmaFByGroup(self, *args, **kwargs)
    def getChiByGroup(self, *args, **kwargs): return _openrk.FunctionalMaterial_getChiByGroup(self, *args, **kwargs)
    def getDifCoefByGroup(self, *args, **kwargs): return _openrk.FunctionalMaterial_getDifCoefByGroup(self, *args, **kwargs)
    def getVelocityByGroup(self, *args, **kwargs): return _openrk.FunctionalMaterial_getVelocityByGroup(self, *args, **kwargs)
    def getDopplerCoefficientByGroup(self, *args, **kwargs): return _openrk.FunctionalMaterial_getDopplerCoefficientByGroup(self, *args, **kwargs)
    def getTimeStep(self, *args, **kwargs): return _openrk.FunctionalMaterial_getTimeStep(self, *args, **kwargs)
    def toString(self): return _openrk.FunctionalMaterial_toString(self)
    def clone(self): return _openrk.FunctionalMaterial_clone(self)
    def copy(self, *args, **kwargs): return _openrk.FunctionalMaterial_copy(self, *args, **kwargs)
FunctionalMaterial_swigregister = _openrk.FunctionalMaterial_swigregister
FunctionalMaterial_swigregister(FunctionalMaterial)

START = _openrk.START
PREVIOUS_OUT = _openrk.PREVIOUS_OUT
PREVIOUS_IN = _openrk.PREVIOUS_IN
CURRENT = _openrk.CURRENT
FORWARD_IN_OLD = _openrk.FORWARD_IN_OLD
FORWARD_OUT = _openrk.FORWARD_OUT
FORWARD_OUT_OLD = _openrk.FORWARD_OUT_OLD
END = _openrk.END
class Clock(_object):
    __swig_setmethods__ = {}
    __setattr__ = lambda self, name, value: _swig_setattr(self, Clock, name, value)
    __swig_getmethods__ = {}
    __getattr__ = lambda self, name: _swig_getattr(self, Clock, name)
    __repr__ = _swig_repr
    def __init__(self, start_time=0.0, end_time=3.0, dt_outer=1.e-1, dt_inner=1.e-2): 
        this = _openrk.new_Clock(start_time, end_time, dt_outer, dt_inner)
        try: self.this.append(this)
        except: self.this = this
    __swig_destroy__ = _openrk.delete_Clock
    __del__ = lambda self : None;
    def getTime(self, *args, **kwargs): return _openrk.Clock_getTime(self, *args, **kwargs)
    def getDtInner(self): return _openrk.Clock_getDtInner(self)
    def getDtOuter(self): return _openrk.Clock_getDtOuter(self)
    def getPositionName(self, *args, **kwargs): return _openrk.Clock_getPositionName(self, *args, **kwargs)
    def setDtOuter(self, *args, **kwargs): return _openrk.Clock_setDtOuter(self, *args, **kwargs)
    def setDtInner(self, *args, **kwargs): return _openrk.Clock_setDtInner(self, *args, **kwargs)
    def setTime(self, *args, **kwargs): return _openrk.Clock_setTime(self, *args, **kwargs)
    def setStartTime(self, *args, **kwargs): return _openrk.Clock_setStartTime(self, *args, **kwargs)
    def setEndTime(self, *args, **kwargs): return _openrk.Clock_setEndTime(self, *args, **kwargs)
    def takeInnerStep(self): return _openrk.Clock_takeInnerStep(self)
    def takeOuterStep(self): return _openrk.Clock_takeOuterStep(self)
    def resetToPreviousOuterStep(self): return _openrk.Clock_resetToPreviousOuterStep(self)
    def toString(self): return _openrk.Clock_toString(self)
    def printString(self): return _openrk.Clock_printString(self)
Clock_swigregister = _openrk.Clock_swigregister
Clock_swigregister(Clock)

REFLECTIVE = _openrk.REFLECTIVE
VACUUM = _openrk.VACUUM
class Mesh(_object):
    __swig_setmethods__ = {}
    __setattr__ = lambda self, name, value: _swig_setattr(self, Mesh, name, value)
    __swig_getmethods__ = {}
    __getattr__ = lambda self, name: _swig_getattr(self, Mesh, name)
    __repr__ = _swig_repr
    def __init__(self, width=1.0, height=1.0, depth=1.0): 
        this = _openrk.new_Mesh(width, height, depth)
        try: self.this.append(this)
        except: self.this = this
    __swig_destroy__ = _openrk.delete_Mesh
    __del__ = lambda self : None;
    def setNumShapeEnergyGroups(self, *args, **kwargs): return _openrk.Mesh_setNumShapeEnergyGroups(self, *args, **kwargs)
    def setNumAmpEnergyGroups(self, *args, **kwargs): return _openrk.Mesh_setNumAmpEnergyGroups(self, *args, **kwargs)
    def setNumDelayedGroups(self, *args, **kwargs): return _openrk.Mesh_setNumDelayedGroups(self, *args, **kwargs)
    def setXMin(self, *args, **kwargs): return _openrk.Mesh_setXMin(self, *args, **kwargs)
    def setXMax(self, *args, **kwargs): return _openrk.Mesh_setXMax(self, *args, **kwargs)
    def setYMin(self, *args, **kwargs): return _openrk.Mesh_setYMin(self, *args, **kwargs)
    def setYMax(self, *args, **kwargs): return _openrk.Mesh_setYMax(self, *args, **kwargs)
    def setZMin(self, *args, **kwargs): return _openrk.Mesh_setZMin(self, *args, **kwargs)
    def setZMax(self, *args, **kwargs): return _openrk.Mesh_setZMax(self, *args, **kwargs)
    def setKeff0(self, *args, **kwargs): return _openrk.Mesh_setKeff0(self, *args, **kwargs)
    def setBoundary(self, *args, **kwargs): return _openrk.Mesh_setBoundary(self, *args, **kwargs)
    def setMaterial(self, *args, **kwargs): return _openrk.Mesh_setMaterial(self, *args, **kwargs)
    def setBuckling(self, *args, **kwargs): return _openrk.Mesh_setBuckling(self, *args, **kwargs)
    def setDecayConstants(self, *args, **kwargs): return _openrk.Mesh_setDecayConstants(self, *args, **kwargs)
    def setDelayedFractions(self, *args, **kwargs): return _openrk.Mesh_setDelayedFractions(self, *args, **kwargs)
    def getKeff0(self): return _openrk.Mesh_getKeff0(self)
    def getClock(self): return _openrk.Mesh_getClock(self)
    def getWidth(self): return _openrk.Mesh_getWidth(self)
    def getHeight(self): return _openrk.Mesh_getHeight(self)
    def getDepth(self): return _openrk.Mesh_getDepth(self)
    def getBoundary(self, *args, **kwargs): return _openrk.Mesh_getBoundary(self, *args, **kwargs)
    def getNumShapeEnergyGroups(self): return _openrk.Mesh_getNumShapeEnergyGroups(self)
    def getNumAmpEnergyGroups(self): return _openrk.Mesh_getNumAmpEnergyGroups(self)
    def getNumDelayedGroups(self): return _openrk.Mesh_getNumDelayedGroups(self)
    def getBuckling(self): return _openrk.Mesh_getBuckling(self)
    def getDecayConstants(self): return _openrk.Mesh_getDecayConstants(self)
    def getDelayedFractions(self): return _openrk.Mesh_getDelayedFractions(self)
    def getDecayConstantByGroup(self, *args, **kwargs): return _openrk.Mesh_getDecayConstantByGroup(self, *args, **kwargs)
    def getDelayedFractionByGroup(self, *args, **kwargs): return _openrk.Mesh_getDelayedFractionByGroup(self, *args, **kwargs)
    def getXMin(self): return _openrk.Mesh_getXMin(self)
    def getXMax(self): return _openrk.Mesh_getXMax(self)
    def getYMin(self): return _openrk.Mesh_getYMin(self)
    def getYMax(self): return _openrk.Mesh_getYMax(self)
    def getZMin(self): return _openrk.Mesh_getZMin(self)
    def getZMax(self): return _openrk.Mesh_getZMax(self)
    def getFlux(self, *args, **kwargs): return _openrk.Mesh_getFlux(self, *args, **kwargs)
    def getPower(self, *args, **kwargs): return _openrk.Mesh_getPower(self, *args, **kwargs)
    def getTemperature(self, *args, **kwargs): return _openrk.Mesh_getTemperature(self, *args, **kwargs)
    def getPowerByValue(self, *args, **kwargs): return _openrk.Mesh_getPowerByValue(self, *args, **kwargs)
    def getTemperatureByValue(self, *args, **kwargs): return _openrk.Mesh_getTemperatureByValue(self, *args, **kwargs)
    def getMaterial(self, *args, **kwargs): return _openrk.Mesh_getMaterial(self, *args, **kwargs)
Mesh_swigregister = _openrk.Mesh_swigregister
Mesh_swigregister(Mesh)

class StructuredMesh(Mesh):
    __swig_setmethods__ = {}
    for _s in [Mesh]: __swig_setmethods__.update(getattr(_s,'__swig_setmethods__',{}))
    __setattr__ = lambda self, name, value: _swig_setattr(self, StructuredMesh, name, value)
    __swig_getmethods__ = {}
    for _s in [Mesh]: __swig_getmethods__.update(getattr(_s,'__swig_getmethods__',{}))
    __getattr__ = lambda self, name: _swig_getattr(self, StructuredMesh, name)
    __repr__ = _swig_repr
    def __init__(self, width=1.0, height=1.0, depth=1.0, num_x=1, num_y=1, num_z=1): 
        this = _openrk.new_StructuredMesh(width, height, depth, num_x, num_y, num_z)
        try: self.this.append(this)
        except: self.this = this
    __swig_destroy__ = _openrk.delete_StructuredMesh
    __del__ = lambda self : None;
    def setNumX(self, *args, **kwargs): return _openrk.StructuredMesh_setNumX(self, *args, **kwargs)
    def setNumY(self, *args, **kwargs): return _openrk.StructuredMesh_setNumY(self, *args, **kwargs)
    def setNumZ(self, *args, **kwargs): return _openrk.StructuredMesh_setNumZ(self, *args, **kwargs)
    def setClock(self, *args, **kwargs): return _openrk.StructuredMesh_setClock(self, *args, **kwargs)
    def getCurrent(self, *args, **kwargs): return _openrk.StructuredMesh_getCurrent(self, *args, **kwargs)
    def getDifLinear(self, *args, **kwargs): return _openrk.StructuredMesh_getDifLinear(self, *args, **kwargs)
    def getDifNonlinear(self, *args, **kwargs): return _openrk.StructuredMesh_getDifNonlinear(self, *args, **kwargs)
    def getNeighborCell(self, *args, **kwargs): return _openrk.StructuredMesh_getNeighborCell(self, *args, **kwargs)
    def getNeighborMaterial(self, *args, **kwargs): return _openrk.StructuredMesh_getNeighborMaterial(self, *args, **kwargs)
    def getNumX(self): return _openrk.StructuredMesh_getNumX(self)
    def getNumY(self): return _openrk.StructuredMesh_getNumY(self)
    def getNumZ(self): return _openrk.StructuredMesh_getNumZ(self)
    def getCellWidth(self): return _openrk.StructuredMesh_getCellWidth(self)
    def getCellHeight(self): return _openrk.StructuredMesh_getCellHeight(self)
    def getCellDepth(self): return _openrk.StructuredMesh_getCellDepth(self)
    def getCellVolume(self): return _openrk.StructuredMesh_getCellVolume(self)
    def findCell(self, *args, **kwargs): return _openrk.StructuredMesh_findCell(self, *args, **kwargs)
    def computeFuelVolume(self): return _openrk.StructuredMesh_computeFuelVolume(self)
    def uniquifyMaterials(self): return _openrk.StructuredMesh_uniquifyMaterials(self)
    def getMaxTemperature(self, *args, **kwargs): return _openrk.StructuredMesh_getMaxTemperature(self, *args, **kwargs)
    def copyPower(self, *args, **kwargs): return _openrk.StructuredMesh_copyPower(self, *args, **kwargs)
    def copyTemperature(self, *args, **kwargs): return _openrk.StructuredMesh_copyTemperature(self, *args, **kwargs)
    def setTemperature(self, *args, **kwargs): return _openrk.StructuredMesh_setTemperature(self, *args, **kwargs)
StructuredMesh_swigregister = _openrk.StructuredMesh_swigregister
StructuredMesh_swigregister(StructuredMesh)

class StructuredShapeMesh(StructuredMesh):
    __swig_setmethods__ = {}
    for _s in [StructuredMesh]: __swig_setmethods__.update(getattr(_s,'__swig_setmethods__',{}))
    __setattr__ = lambda self, name, value: _swig_setattr(self, StructuredShapeMesh, name, value)
    __swig_getmethods__ = {}
    for _s in [StructuredMesh]: __swig_getmethods__.update(getattr(_s,'__swig_getmethods__',{}))
    __getattr__ = lambda self, name: _swig_getattr(self, StructuredShapeMesh, name)
    __repr__ = _swig_repr
    def __init__(self, width=1.0, height=1.0, depth=1.0, num_x=1, num_y=1, num_z=1): 
        this = _openrk.new_StructuredShapeMesh(width, height, depth, num_x, num_y, num_z)
        try: self.this.append(this)
        except: self.this = this
    __swig_destroy__ = _openrk.delete_StructuredShapeMesh
    __del__ = lambda self : None;
    def setAmpMesh(self, *args, **kwargs): return _openrk.StructuredShapeMesh_setAmpMesh(self, *args, **kwargs)
    def setFluxByValue(self, *args, **kwargs): return _openrk.StructuredShapeMesh_setFluxByValue(self, *args, **kwargs)
    def setCurrentByValue(self, *args, **kwargs): return _openrk.StructuredShapeMesh_setCurrentByValue(self, *args, **kwargs)
    def setDifLinearByValue(self, *args, **kwargs): return _openrk.StructuredShapeMesh_setDifLinearByValue(self, *args, **kwargs)
    def setDifNonlinearByValue(self, *args, **kwargs): return _openrk.StructuredShapeMesh_setDifNonlinearByValue(self, *args, **kwargs)
    def setGroupStructure(self, *args, **kwargs): return _openrk.StructuredShapeMesh_setGroupStructure(self, *args, **kwargs)
    def getFluxByValue(self, *args, **kwargs): return _openrk.StructuredShapeMesh_getFluxByValue(self, *args, **kwargs)
    def getCurrentByValue(self, *args, **kwargs): return _openrk.StructuredShapeMesh_getCurrentByValue(self, *args, **kwargs)
    def getDifLinearByValue(self, *args, **kwargs): return _openrk.StructuredShapeMesh_getDifLinearByValue(self, *args, **kwargs)
    def getDifNonlinearByValue(self, *args, **kwargs): return _openrk.StructuredShapeMesh_getDifNonlinearByValue(self, *args, **kwargs)
    def getAmpGroup(self, *args, **kwargs): return _openrk.StructuredShapeMesh_getAmpGroup(self, *args, **kwargs)
    def clone(self): return _openrk.StructuredShapeMesh_clone(self)
    def uniformRefine(self, refine_x=1, refine_y=1, refine_z=1): return _openrk.StructuredShapeMesh_uniformRefine(self, refine_x, refine_y, refine_z)
    def initialize(self): return _openrk.StructuredShapeMesh_initialize(self)
    def synthesizeFlux(self, *args, **kwargs): return _openrk.StructuredShapeMesh_synthesizeFlux(self, *args, **kwargs)
    def reconstructFlux(self, *args, **kwargs): return _openrk.StructuredShapeMesh_reconstructFlux(self, *args, **kwargs)
    def computePower(self, *args, **kwargs): return _openrk.StructuredShapeMesh_computePower(self, *args, **kwargs)
    def computeInitialPrecursorConc(self, *args, **kwargs): return _openrk.StructuredShapeMesh_computeInitialPrecursorConc(self, *args, **kwargs)
    def integratePrecursorConc(self, *args, **kwargs): return _openrk.StructuredShapeMesh_integratePrecursorConc(self, *args, **kwargs)
    def integrateTemperature(self, *args, **kwargs): return _openrk.StructuredShapeMesh_integrateTemperature(self, *args, **kwargs)
    def computeDifCoefs(self, *args, **kwargs): return _openrk.StructuredShapeMesh_computeDifCoefs(self, *args, **kwargs)
    def copyFlux(self, *args, **kwargs): return _openrk.StructuredShapeMesh_copyFlux(self, *args, **kwargs)
    def copyCurrent(self, *args, **kwargs): return _openrk.StructuredShapeMesh_copyCurrent(self, *args, **kwargs)
    def copyDifLinear(self, *args, **kwargs): return _openrk.StructuredShapeMesh_copyDifLinear(self, *args, **kwargs)
    def copyDifNonlinear(self, *args, **kwargs): return _openrk.StructuredShapeMesh_copyDifNonlinear(self, *args, **kwargs)
    def scaleFlux(self, *args, **kwargs): return _openrk.StructuredShapeMesh_scaleFlux(self, *args, **kwargs)
    def computeAveragePower(self, *args, **kwargs): return _openrk.StructuredShapeMesh_computeAveragePower(self, *args, **kwargs)
    def computePowerL2Norm(self, *args, **kwargs): return _openrk.StructuredShapeMesh_computePowerL2Norm(self, *args, **kwargs)
    def findAmpCell(self, *args, **kwargs): return _openrk.StructuredShapeMesh_findAmpCell(self, *args, **kwargs)
    def saveShape(self): return _openrk.StructuredShapeMesh_saveShape(self)
StructuredShapeMesh_swigregister = _openrk.StructuredShapeMesh_swigregister
StructuredShapeMesh_swigregister(StructuredShapeMesh)

class AmpMesh(StructuredMesh):
    __swig_setmethods__ = {}
    for _s in [StructuredMesh]: __swig_setmethods__.update(getattr(_s,'__swig_setmethods__',{}))
    __setattr__ = lambda self, name, value: _swig_setattr(self, AmpMesh, name, value)
    __swig_getmethods__ = {}
    for _s in [StructuredMesh]: __swig_getmethods__.update(getattr(_s,'__swig_getmethods__',{}))
    __getattr__ = lambda self, name: _swig_getattr(self, AmpMesh, name)
    __repr__ = _swig_repr
    def __init__(self, width=1.0, height=1.0, depth=1.0, num_x=1, num_y=1, num_z=1): 
        this = _openrk.new_AmpMesh(width, height, depth, num_x, num_y, num_z)
        try: self.this.append(this)
        except: self.this = this
    __swig_destroy__ = _openrk.delete_AmpMesh
    __del__ = lambda self : None;
    def setOpticallyThick(self, *args, **kwargs): return _openrk.AmpMesh_setOpticallyThick(self, *args, **kwargs)
    def setShapeMesh(self, *args, **kwargs): return _openrk.AmpMesh_setShapeMesh(self, *args, **kwargs)
    def setFluxByValue(self, *args, **kwargs): return _openrk.AmpMesh_setFluxByValue(self, *args, **kwargs)
    def setCurrentByValue(self, *args, **kwargs): return _openrk.AmpMesh_setCurrentByValue(self, *args, **kwargs)
    def setDifLinearByValue(self, *args, **kwargs): return _openrk.AmpMesh_setDifLinearByValue(self, *args, **kwargs)
    def setDifNonlinearByValue(self, *args, **kwargs): return _openrk.AmpMesh_setDifNonlinearByValue(self, *args, **kwargs)
    def setGroupStructure(self, *args, **kwargs): return _openrk.AmpMesh_setGroupStructure(self, *args, **kwargs)
    def getFluxByValue(self, *args, **kwargs): return _openrk.AmpMesh_getFluxByValue(self, *args, **kwargs)
    def getCurrentByValue(self, *args, **kwargs): return _openrk.AmpMesh_getCurrentByValue(self, *args, **kwargs)
    def getDifLinearByValue(self, *args, **kwargs): return _openrk.AmpMesh_getDifLinearByValue(self, *args, **kwargs)
    def getDifNonlinearByValue(self, *args, **kwargs): return _openrk.AmpMesh_getDifNonlinearByValue(self, *args, **kwargs)
    def clone(self): return _openrk.AmpMesh_clone(self)
    def initialize(self): return _openrk.AmpMesh_initialize(self)
    def condenseMaterials(self, *args, **kwargs): return _openrk.AmpMesh_condenseMaterials(self, *args, **kwargs)
    def computePower(self, *args, **kwargs): return _openrk.AmpMesh_computePower(self, *args, **kwargs)
    def computeCurrent(self, *args, **kwargs): return _openrk.AmpMesh_computeCurrent(self, *args, **kwargs)
    def computeDifCorrect(self, *args, **kwargs): return _openrk.AmpMesh_computeDifCorrect(self, *args, **kwargs)
    def computeDifCoefs(self, *args, **kwargs): return _openrk.AmpMesh_computeDifCoefs(self, *args, **kwargs)
    def copyFlux(self, *args, **kwargs): return _openrk.AmpMesh_copyFlux(self, *args, **kwargs)
    def copyCurrent(self, *args, **kwargs): return _openrk.AmpMesh_copyCurrent(self, *args, **kwargs)
    def copyDifLinear(self, *args, **kwargs): return _openrk.AmpMesh_copyDifLinear(self, *args, **kwargs)
    def copyDifNonlinear(self, *args, **kwargs): return _openrk.AmpMesh_copyDifNonlinear(self, *args, **kwargs)
    def computeAveragePower(self, *args, **kwargs): return _openrk.AmpMesh_computeAveragePower(self, *args, **kwargs)
    def computePowerL2Norm(self, *args, **kwargs): return _openrk.AmpMesh_computePowerL2Norm(self, *args, **kwargs)
    def interpolateDifNonlinear(self, *args, **kwargs): return _openrk.AmpMesh_interpolateDifNonlinear(self, *args, **kwargs)
AmpMesh_swigregister = _openrk.AmpMesh_swigregister
AmpMesh_swigregister(AmpMesh)

class Solver(_object):
    __swig_setmethods__ = {}
    __setattr__ = lambda self, name, value: _swig_setattr(self, Solver, name, value)
    __swig_getmethods__ = {}
    __getattr__ = lambda self, name: _swig_getattr(self, Solver, name)
    __repr__ = _swig_repr
    def __init__(self, *args, **kwargs): 
        this = _openrk.new_Solver(*args, **kwargs)
        try: self.this.append(this)
        except: self.this = this
    __swig_destroy__ = _openrk.delete_Solver
    __del__ = lambda self : None;
    def getAMShape(self): return _openrk.Solver_getAMShape(self)
    def getAShape(self): return _openrk.Solver_getAShape(self)
    def getMShape(self): return _openrk.Solver_getMShape(self)
    def getBShape(self): return _openrk.Solver_getBShape(self)
    def getAMAmp(self): return _openrk.Solver_getAMAmp(self)
    def getAAmp(self): return _openrk.Solver_getAAmp(self)
    def getMAmp(self): return _openrk.Solver_getMAmp(self)
    def getBAmp(self): return _openrk.Solver_getBAmp(self)
    def makeAMShapeInitial(self, *args, **kwargs): return _openrk.Solver_makeAMShapeInitial(self, *args, **kwargs)
    def computeInitialShape(self, *args, **kwargs): return _openrk.Solver_computeInitialShape(self, *args, **kwargs)
    def makeAMAmp(self, *args, **kwargs): return _openrk.Solver_makeAMAmp(self, *args, **kwargs)
    def makeAMShape(self, *args, **kwargs): return _openrk.Solver_makeAMShape(self, *args, **kwargs)
    def computeAmpFrequency(self): return _openrk.Solver_computeAmpFrequency(self)
Solver_swigregister = _openrk.Solver_swigregister
Solver_swigregister(Solver)

FORWARD_EULER = _openrk.FORWARD_EULER
BACKWARD_EULER = _openrk.BACKWARD_EULER
CRANK_NICOLSON = _openrk.CRANK_NICOLSON
CUSTOM = _openrk.CUSTOM
class Transient(_object):
    __swig_setmethods__ = {}
    __setattr__ = lambda self, name, value: _swig_setattr(self, Transient, name, value)
    __swig_getmethods__ = {}
    __getattr__ = lambda self, name: _swig_getattr(self, Transient, name)
    __repr__ = _swig_repr
    def __init__(self, *args, **kwargs): 
        this = _openrk.new_Transient(*args, **kwargs)
        try: self.this.append(this)
        except: self.this = this
    __swig_destroy__ = _openrk.delete_Transient
    __del__ = lambda self : None;
    def setInnerMethod(self, *args, **kwargs): return _openrk.Transient_setInnerMethod(self, *args, **kwargs)
    def setOuterMethod(self, *args, **kwargs): return _openrk.Transient_setOuterMethod(self, *args, **kwargs)
    def setInitialPower(self, *args, **kwargs): return _openrk.Transient_setInitialPower(self, *args, **kwargs)
    def setClock(self, *args, **kwargs): return _openrk.Transient_setClock(self, *args, **kwargs)
    def setSolver(self, *args, **kwargs): return _openrk.Transient_setSolver(self, *args, **kwargs)
    def setShapeMesh(self, *args, **kwargs): return _openrk.Transient_setShapeMesh(self, *args, **kwargs)
    def setAmpMesh(self, *args, **kwargs): return _openrk.Transient_setAmpMesh(self, *args, **kwargs)
    def computeInitialShape(self): return _openrk.Transient_computeInitialShape(self)
    def takeInnerStep(self): return _openrk.Transient_takeInnerStep(self)
    def takeOuterStep(self, tol=1.e-6): return _openrk.Transient_takeOuterStep(self, tol)
    def takeOuterStepOnly(self, tol=1.e-6): return _openrk.Transient_takeOuterStepOnly(self, tol)
    def broadcastToActive(self, *args, **kwargs): return _openrk.Transient_broadcastToActive(self, *args, **kwargs)
    def broadcastToAll(self, *args, **kwargs): return _openrk.Transient_broadcastToAll(self, *args, **kwargs)
    def broadcastToOne(self, *args, **kwargs): return _openrk.Transient_broadcastToOne(self, *args, **kwargs)
Transient_swigregister = _openrk.Transient_swigregister
Transient_swigregister(Transient)

# This file is compatible with both classic and new-style classes.


