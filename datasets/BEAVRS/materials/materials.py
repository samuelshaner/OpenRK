from isotopes import *


# Keys are string material names and values are OpenMOC IsoMaterials filled
# with Isotopes
materials = {}


################################################################################
############################  1.6% Enriched Fuel  ##############################
################################################################################

materials['1.6% Fuel'] = IsoMaterial(material_id(), '1.6% Fuel')
materials['1.6% Fuel'].setNumEnergyGroups(1)
materials['1.6% Fuel'].setDensity(10.31341, 'g/cc')
materials['1.6% Fuel'].addIsotope(isotopes['U-234'], 3.0131e-6)
materials['1.6% Fuel'].addIsotope(isotopes['U-235'], 3.7503e-4)
materials['1.6% Fuel'].addIsotope(isotopes['U-238'], 2.2625e-2)

# Includes O-17 and O-18 number densities
materials['1.6% Fuel'].addIsotope(isotopes['O-16'], 4.6007e-2)


################################################################################
############################  2.4% Enriched Fuel  ##############################
################################################################################

materials['2.4% Fuel'] = IsoMaterial(material_id(), '2.4% Fuel')
materials['2.4% Fuel'].setNumEnergyGroups(1)
materials['2.4% Fuel'].setDensity(10.29748, 'g/cc')
materials['2.4% Fuel'].addIsotope(isotopes['U-234'], 4.4842E-6)
materials['2.4% Fuel'].addIsotope(isotopes['U-235'], 5.5814e-4)
materials['2.4% Fuel'].addIsotope(isotopes['U-238'], 2.2407e-2)

# Includes O-17 and O-18 number densities
materials['2.4% Fuel'].addIsotope(isotopes['O-16'], 4.5940e-2)


################################################################################
############################  3.1% Enriched Fuel  ##############################
################################################################################

materials['3.1% Fuel'] = IsoMaterial(material_id(), '3.1% Fuel')
materials['3.1% Fuel'].setNumEnergyGroups(1)
materials['3.1% Fuel'].setDensity(10.30166, 'g/cc')
materials['3.1% Fuel'].addIsotope(isotopes['U-234'], 5.7987e-6)
materials['3.1% Fuel'].addIsotope(isotopes['U-235'], 7.2175e-4)
materials['3.1% Fuel'].addIsotope(isotopes['U-238'], 2.2253e-2)

# Includes O-17 and O-18 number densities
materials['3.1% Fuel'].addIsotope(isotopes['O-16'], 4.5940e-2)


################################################################################
###############################  Borated Water  ################################
################################################################################

materials['Borated Water'] = IsoMaterial(material_id(), 'Borated Water')
materials['Borated Water'].setNumEnergyGroups(1)
materials['Borated Water'].setDensity(0.740582, 'g/cc')
materials['Borated Water'].addIsotope(isotopes['B-10'], 8.0042e-6)
materials['Borated Water'].addIsotope(isotopes['B-11'], 3.2218e-5)
materials['Borated Water'].addIsotope(isotopes['H-1'], 4.9457e-2)
materials['Borated Water'].addIsotope(isotopes['H-2'], 7.4196e-6)

# Includes O-17 and O-18 number densities
materials['Borated Water'].addIsotope(isotopes['O-16'], 2.4732e-2)


################################################################################
####################################  Gap  #####################################
################################################################################

materials['Gap'] = IsoMaterial(material_id(), 'Gap')
materials['Gap'].setNumEnergyGroups(1)
materials['Gap'].setDensity(0.001598, 'g/cc')
materials['Gap'].addIsotope(isotopes['He-4'], 2.4044e-4)


################################################################################
##################################  Zircaloy  ##################################
################################################################################

materials['Zircaloy'] = IsoMaterial(material_id(), 'Zircaloy')
materials['Zircaloy'].setNumEnergyGroups(1)
materials['Zircaloy'].setDensity(6.55, 'g/cc')
materials['Zircaloy'].addIsotope(isotopes['Cr-50'], 3.2962e-6)
materials['Zircaloy'].addIsotope(isotopes['Cr-52'], 6.3564e-5)
materials['Zircaloy'].addIsotope(isotopes['Cr-53'], 7.2076e-6)
materials['Zircaloy'].addIsotope(isotopes['Cr-54'], 1.7941e-6)
materials['Zircaloy'].addIsotope(isotopes['Fe-54'], 8.6699e-6)
materials['Zircaloy'].addIsotope(isotopes['Fe-56'], 1.3610e-4)
materials['Zircaloy'].addIsotope(isotopes['Fe-57'], 3.1431e-6)
materials['Zircaloy'].addIsotope(isotopes['Fe-58'], 4.1829e-7)
materials['Zircaloy'].addIsotope(isotopes['Zr-90'], 2.1827e-2)
materials['Zircaloy'].addIsotope(isotopes['Zr-91'], 4.7600e-3)
materials['Zircaloy'].addIsotope(isotopes['Zr-92'], 7.2758e-3)
materials['Zircaloy'].addIsotope(isotopes['Zr-94'], 7.3734e-3)
materials['Zircaloy'].addIsotope(isotopes['Zr-96'], 1.1879e-3)
materials['Zircaloy'].addIsotope(isotopes['Sn-112'], 4.6735e-6)
materials['Zircaloy'].addIsotope(isotopes['Sn-114'], 3.1799e-6)
materials['Zircaloy'].addIsotope(isotopes['Sn-115'], 1.6381e-6)
materials['Zircaloy'].addIsotope(isotopes['Sn-116'], 7.0055e-5)
materials['Zircaloy'].addIsotope(isotopes['Sn-117'], 3.7003e-5)
materials['Zircaloy'].addIsotope(isotopes['Sn-118'], 1.1669e-4)
materials['Zircaloy'].addIsotope(isotopes['Sn-119'], 4.1387e-5)
materials['Zircaloy'].addIsotope(isotopes['Sn-120'], 1.5697e-4)
materials['Zircaloy'].addIsotope(isotopes['Sn-122'], 2.2308e-5)
materials['Zircaloy'].addIsotope(isotopes['Sn-124'], 2.7897e-5)

# Includes O-17 and O-18 number densities
materials['Zircaloy'].addIsotope(isotopes['O-16'], 3.0818e-4)


################################################################################
####################################  Air  #####################################
################################################################################

materials['Air'] = IsoMaterial(material_id(), 'Air')
materials['Air'].setNumEnergyGroups(1)
materials['Air'].setDensity(0.000616, 'g/cc')
materials['Air'].addIsotope(isotopes['N-14'], 1.9681e-5)
materials['Air'].addIsotope(isotopes['N-15'], 7.1900e-8)
materials['Air'].addIsotope(isotopes['Ar-36'], 7.9414e-10)
materials['Air'].addIsotope(isotopes['Ar-38'], 1.4915e-10)
materials['Air'].addIsotope(isotopes['Ar-40'], 2.3506e-7)

# Includes O-17 and O-18 number densities
materials['Air'].addIsotope(isotopes['O-16'], 5.2993e-6)

# Unavailable in cross-section library?
#materials['Air'].addIsotope(isotopes['C-12'], 6.7565e-9)
#materials['Air'].addIsotope(isotopes['C-13'], 7.3706e-11)


################################################################################
#############################  Borosilicate Glass  #############################
################################################################################

materials['Boro. Glass'] = IsoMaterial(material_id(), 'Borosilicate Glass')
materials['Boro. Glass'].setNumEnergyGroups(1)
materials['Boro. Glass'].setDensity(2.26000, 'g/cc')
materials['Boro. Glass'].addIsotope(isotopes['B-10'], 9.6506e-4)
materials['Boro. Glass'].addIsotope(isotopes['B-11'], 3.9189e-3)
materials['Boro. Glass'].addIsotope(isotopes['Al-27'], 1.7352e-3)
materials['Boro. Glass'].addIsotope(isotopes['Si-28'], 1.6924e-2)
materials['Boro. Glass'].addIsotope(isotopes['Si-29'], 8.5977e-4)
materials['Boro. Glass'].addIsotope(isotopes['Si-30'], 5.6743e-4)

# Includes O-17 and O-18 number densities
materials['Boro. Glass'].addIsotope(isotopes['O-16'], 4.6624e-2)


################################################################################
##############################  Stainless Steel  ###############################
################################################################################

materials['Steel'] = IsoMaterial(material_id(), 'Stainless Steel')
materials['Steel'].setNumEnergyGroups(1)
materials['Steel'].setDensity(8.0300, 'g/cc')
materials['Steel'].addIsotope(isotopes['Si-28'], 9.5274e-4)
materials['Steel'].addIsotope(isotopes['Si-29'], 4.8400e-5)
materials['Steel'].addIsotope(isotopes['Si-30'], 3.1943e-5)
materials['Steel'].addIsotope(isotopes['Cr-50'], 7.6778e-4)
materials['Steel'].addIsotope(isotopes['Cr-52'], 1.4806e-2)
materials['Steel'].addIsotope(isotopes['Cr-53'], 1.6789e-3)
materials['Steel'].addIsotope(isotopes['Cr-54'], 4.1791e-4)
materials['Steel'].addIsotope(isotopes['Mn-55'], 1.7604e-3)
materials['Steel'].addIsotope(isotopes['Fe-54'], 3.4620e-3)
materials['Steel'].addIsotope(isotopes['Fe-56'], 5.4345e-2)
materials['Steel'].addIsotope(isotopes['Fe-57'], 1.2551e-3)
materials['Steel'].addIsotope(isotopes['Fe-58'], 1.6703e-4)
materials['Steel'].addIsotope(isotopes['Ni-58'], 5.6089e-3)
materials['Steel'].addIsotope(isotopes['Ni-60'], 2.1605e-3)
materials['Steel'].addIsotope(isotopes['Ni-61'], 9.3917e-5)
materials['Steel'].addIsotope(isotopes['Ni-62'], 2.9945e-4)
materials['Steel'].addIsotope(isotopes['Ni-64'], 7.6261e-5)