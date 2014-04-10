from opencsg import *

# Keys are Surface string names and values are OpenCSG Surfaces.
surfaces = dict()


################################################################################
#############################  Fuel Pin Cell  ##################################
################################################################################

surfaces['Fuel Radius-1'] = ZCylinder(x0=0.0, y0=0.0, R=0.39218)
surfaces['Fuel Radius-2'] = ZCylinder(x0=0.0, y0=0.0, R=0.40005)
surfaces['Fuel Radius-3'] = ZCylinder(x0=0.0, y0=0.0, R=0.45720)


################################################################################
##########################  Guide Tube Pin Cell  ###############################
################################################################################

surfaces['GT Radius-1'] = ZCylinder(x0=0.0, y0=0.0, R=0.56134)
surfaces['GT Radius-2'] = ZCylinder(x0=0.0, y0=0.0, R=0.60198)


################################################################################
#######################  Instrument Tube Pin Cell  #############################
################################################################################

surfaces['IT Radius-1'] = ZCylinder(x0=0.0, y0=0.0, R=0.43688)
surfaces['IT Radius-2'] = ZCylinder(x0=0.0, y0=0.0, R=0.48387)


################################################################################
#######################  Burnable Absorber Pin Cell  ###########################
################################################################################

surfaces['BA Radius-1'] = ZCylinder(x0=0.0, y0=0.0, R=0.21400)
surfaces['BA Radius-2'] = ZCylinder(x0=0.0, y0=0.0, R=0.23051)
surfaces['BA Radius-3'] = ZCylinder(x0=0.0, y0=0.0, R=0.24130)
surfaces['BA Radius-4'] = ZCylinder(x0=0.0, y0=0.0, R=0.42672)
surfaces['BA Radius-5'] = ZCylinder(x0=0.0, y0=0.0, R=0.43688)
surfaces['BA Radius-6'] = ZCylinder(x0=0.0, y0=0.0, R=0.48387)
surfaces['BA Radius-7'] = ZCylinder(x0=0.0, y0=0.0, R=0.56134)
surfaces['BA Radius-8'] = ZCylinder(x0=0.0, y0=0.0, R=0.60198)