from openmoc import *

# Keys are Isotope string names and values are OpenMOC Isotopes.
isotopes = {}

################################################################################
###########################  Initialize Isotopes  ##############################
################################################################################

isotopes['H-1'] = Isotope(isotope_id(), 'H-1')
isotopes['H-2'] = Isotope(isotope_id(), 'H-2')
isotopes['He-4'] = Isotope(isotope_id(), 'He-4')
isotopes['B-10'] = Isotope(isotope_id(), 'B-10')
isotopes['B-11'] = Isotope(isotope_id(), 'B-11')
isotopes['C-12'] = Isotope(isotope_id(), 'C-12')
isotopes['C-13'] = Isotope(isotope_id(), 'C-13')
isotopes['N-14'] = Isotope(isotope_id(), 'N-14')
isotopes['N-15'] = Isotope(isotope_id(), 'N-15')
isotopes['O-16'] = Isotope(isotope_id(), 'O-16')
isotopes['Al-27'] = Isotope(isotope_id(), 'Al-27')
isotopes['Si-28'] = Isotope(isotope_id(), 'Si-28')
isotopes['Si-29'] = Isotope(isotope_id(), 'Si-29')
isotopes['Si-30'] = Isotope(isotope_id(), 'Si-30')
isotopes['Ar-36'] = Isotope(isotope_id(), 'Ar-36')
isotopes['Ar-38'] = Isotope(isotope_id(), 'Ar-38')
isotopes['Ar-40'] = Isotope(isotope_id(), 'Ar-40')
isotopes['Cr-50'] = Isotope(isotope_id(), 'Cr-50')
isotopes['Cr-52'] = Isotope(isotope_id(), 'Cr-52')
isotopes['Cr-53'] = Isotope(isotope_id(), 'Cr-53')
isotopes['Cr-54'] = Isotope(isotope_id(), 'Cr-54')
isotopes['Fe-54'] = Isotope(isotope_id(), 'Fe-54')
isotopes['Fe-56'] = Isotope(isotope_id(), 'Fe-56')
isotopes['Fe-57'] = Isotope(isotope_id(), 'Fe-57')
isotopes['Fe-58'] = Isotope(isotope_id(), 'Fe-58')
isotopes['Mn-55'] = Isotope(isotope_id(), 'Mn-55')
isotopes['Ni-58'] = Isotope(isotope_id(), 'Ni-58')
isotopes['Ni-60'] = Isotope(isotope_id(), 'Ni-60')
isotopes['Ni-61'] = Isotope(isotope_id(), 'Ni-61')
isotopes['Ni-62'] = Isotope(isotope_id(), 'Ni-62')
isotopes['Ni-64'] = Isotope(isotope_id(), 'Ni-64')
isotopes['Zr-90'] = Isotope(isotope_id(), 'Zr-90')
isotopes['Zr-91'] = Isotope(isotope_id(), 'Zr-91')
isotopes['Zr-92'] = Isotope(isotope_id(), 'Zr-92')
isotopes['Zr-94'] = Isotope(isotope_id(), 'Zr-94')
isotopes['Zr-96'] = Isotope(isotope_id(), 'Zr-96')
isotopes['Sn-112'] = Isotope(isotope_id(), 'Sn-112')
isotopes['Sn-114'] = Isotope(isotope_id(), 'Sn-114')
isotopes['Sn-115'] = Isotope(isotope_id(), 'Sn-115')
isotopes['Sn-116'] = Isotope(isotope_id(), 'Sn-116')
isotopes['Sn-117'] = Isotope(isotope_id(), 'Sn-117')
isotopes['Sn-118'] = Isotope(isotope_id(), 'Sn-118')
isotopes['Sn-119'] = Isotope(isotope_id(), 'Sn-119')
isotopes['Sn-120'] = Isotope(isotope_id(), 'Sn-120')
isotopes['Sn-122'] = Isotope(isotope_id(), 'Sn-122')
isotopes['Sn-124'] = Isotope(isotope_id(), 'Sn-124')
isotopes['U-234'] = Isotope(isotope_id(), 'U-234')
isotopes['U-235'] = Isotope(isotope_id(), 'U-235')
isotopes['U-238'] = Isotope(isotope_id(), 'U-238')


################################################################################
########################  Set Dummy Energy Groups  #############################
################################################################################

# Set the number of energy groups for each isotope to a dummy value to
# satisfy OpenMOC's error checking requirements
for key in isotopes.keys():
  isotopes[key].setNumEnergyGroups(1)