from openmoc import *

# Keys are Surface string names and values are OpenMOC Surfaces.
surfaces = {}


################################################################################
#############################  Fuel Pin Cell  ##################################
################################################################################

surfaces['Fuel Radius'] = Circle(x=0.0, y=0.0, radius=0.39218)
surfaces['Gap Radius'] = Circle(x=0.0, y=0.0, radius=0.40005)


################################################################################
##########################  Guide Tube Pin Cell  ###############################
################################################################################

surfaces['Fuel Clad Radius'] = Circle(x=0.0, y=0.0, radius=0.45720)
surfaces['GT Radius'] = Circle(x=0.0, y=0.0, radius=0.50419)
surfaces['GT Clad Radius'] = Circle(x=0.0, y=0.0, radius=0.54610)


################################################################################
#######################  Burnable Absorber Pin Cell  ###########################
################################################################################

surfaces['BA Radius 1'] = Circle(x=0.0, y=0.0, radius=0.21400)
surfaces['BA Radius 2'] = Circle(x=0.0, y=0.0, radius=0.23051)
surfaces['BA Radius 3'] = Circle(x=0.0, y=0.0, radius=0.24130)
surfaces['BA Radius 4'] = Circle(x=0.0, y=0.0, radius=0.42672)
surfaces['BA Radius 5'] = Circle(x=0.0, y=0.0, radius=0.43688)
surfaces['BA Radius 6'] = Circle(x=0.0, y=0.0, radius=0.48387)
surfaces['BA Radius 7'] = Circle(x=0.0, y=0.0, radius=0.56134)
surfaces['BA Radius 8'] = Circle(x=0.0, y=0.0, radius=0.60198)