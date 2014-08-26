from openmc.material import Material
from openmc.opencsg_compatible import get_opencsg_material
from datasets.BEAVRS import nuclides


# Keys are string material names and values are OpenMC Materials filled
# with nuclides
openmc_materials = dict()

# Keys are string material names and values are OpenCSG Materials
opencsg_materials = dict()


################################################################################
############################  1.6% Enriched Fuel  ##############################
################################################################################

openmc_materials['1.6% Fuel'] = Material(name='1.6% Fuel')
openmc_materials['1.6% Fuel'].setDensity('g/cm3', 10.31341)
openmc_materials['1.6% Fuel'].addNuclide(nuclides['U-234'], 3.0131e-6)
openmc_materials['1.6% Fuel'].addNuclide(nuclides['U-235'], 3.7503e-4)
openmc_materials['1.6% Fuel'].addNuclide(nuclides['U-238'], 2.2625e-2)

# Includes O-17 and O-18 number densities
openmc_materials['1.6% Fuel'].addNuclide(nuclides['O-16'], 4.6007e-2)


################################################################################
############################  2.4% Enriched Fuel  ##############################
################################################################################

openmc_materials['2.4% Fuel'] = Material(name='2.4% Fuel')
openmc_materials['2.4% Fuel'].setDensity('g/cm3', 10.29748)
openmc_materials['2.4% Fuel'].addNuclide(nuclides['U-234'], 4.4842E-6)
openmc_materials['2.4% Fuel'].addNuclide(nuclides['U-235'], 5.5814e-4)
openmc_materials['2.4% Fuel'].addNuclide(nuclides['U-238'], 2.2407e-2)

# Includes O-17 and O-18 number densities
openmc_materials['2.4% Fuel'].addNuclide(nuclides['O-16'], 4.5940e-2)


################################################################################
############################  3.1% Enriched Fuel  ##############################
################################################################################

openmc_materials['3.1% Fuel'] = Material(name='3.1% Fuel')
openmc_materials['3.1% Fuel'].setDensity('g/cm3', 10.30166)
openmc_materials['3.1% Fuel'].addNuclide(nuclides['U-234'], 5.7987e-6)
openmc_materials['3.1% Fuel'].addNuclide(nuclides['U-235'], 7.2175e-4)
openmc_materials['3.1% Fuel'].addNuclide(nuclides['U-238'], 2.2253e-2)

# Includes O-17 and O-18 number densities
openmc_materials['3.1% Fuel'].addNuclide(nuclides['O-16'], 4.5940e-2)


################################################################################
###############################  Borated Water  ################################
################################################################################

openmc_materials['Borated Water'] = Material(name='Borated Water')
openmc_materials['Borated Water'].setDensity('g/cm3', 0.740582)
openmc_materials['Borated Water'].addNuclide(nuclides['B-10'], 8.0042e-6)
openmc_materials['Borated Water'].addNuclide(nuclides['B-11'], 3.2218e-5)
openmc_materials['Borated Water'].addNuclide(nuclides['H-1'], 4.9457e-2)
openmc_materials['Borated Water'].addNuclide(nuclides['H-2'], 7.4196e-6)

# Includes O-17 and O-18 number densities
openmc_materials['Borated Water'].addNuclide(nuclides['O-16'], 2.4732e-2)


################################################################################
####################################  Gap  #####################################
################################################################################

openmc_materials['Helium'] = Material(name='Gap')
openmc_materials['Helium'].setDensity('g/cm3', 0.001598)
openmc_materials['Helium'].addNuclide(nuclides['He-4'], 2.4044e-4)


################################################################################
##################################  Zircaloy  ##################################
################################################################################

openmc_materials['Zircaloy'] = Material(name='Zircaloy')
openmc_materials['Zircaloy'].setDensity('g/cm3', 6.55)
openmc_materials['Zircaloy'].addNuclide(nuclides['Cr-50'], 3.2962e-6)
openmc_materials['Zircaloy'].addNuclide(nuclides['Cr-52'], 6.3564e-5)
openmc_materials['Zircaloy'].addNuclide(nuclides['Cr-53'], 7.2076e-6)
openmc_materials['Zircaloy'].addNuclide(nuclides['Cr-54'], 1.7941e-6)
openmc_materials['Zircaloy'].addNuclide(nuclides['Fe-54'], 8.6699e-6)
openmc_materials['Zircaloy'].addNuclide(nuclides['Fe-56'], 1.3610e-4)
openmc_materials['Zircaloy'].addNuclide(nuclides['Fe-57'], 3.1431e-6)
openmc_materials['Zircaloy'].addNuclide(nuclides['Fe-58'], 4.1829e-7)
openmc_materials['Zircaloy'].addNuclide(nuclides['Zr-90'], 2.1827e-2)
openmc_materials['Zircaloy'].addNuclide(nuclides['Zr-91'], 4.7600e-3)
openmc_materials['Zircaloy'].addNuclide(nuclides['Zr-92'], 7.2758e-3)
openmc_materials['Zircaloy'].addNuclide(nuclides['Zr-94'], 7.3734e-3)
openmc_materials['Zircaloy'].addNuclide(nuclides['Zr-96'], 1.1879e-3)
openmc_materials['Zircaloy'].addNuclide(nuclides['Sn-112'], 4.6735e-6)
openmc_materials['Zircaloy'].addNuclide(nuclides['Sn-114'], 3.1799e-6)
openmc_materials['Zircaloy'].addNuclide(nuclides['Sn-115'], 1.6381e-6)
openmc_materials['Zircaloy'].addNuclide(nuclides['Sn-116'], 7.0055e-5)
openmc_materials['Zircaloy'].addNuclide(nuclides['Sn-117'], 3.7003e-5)
openmc_materials['Zircaloy'].addNuclide(nuclides['Sn-118'], 1.1669e-4)
openmc_materials['Zircaloy'].addNuclide(nuclides['Sn-119'], 4.1387e-5)
openmc_materials['Zircaloy'].addNuclide(nuclides['Sn-120'], 1.5697e-4)
openmc_materials['Zircaloy'].addNuclide(nuclides['Sn-122'], 2.2308e-5)
openmc_materials['Zircaloy'].addNuclide(nuclides['Sn-124'], 2.7897e-5)

# Includes O-17 and O-18 number densities
openmc_materials['Zircaloy'].addNuclide(nuclides['O-16'], 3.0818e-4)


################################################################################
####################################  Air  #####################################
################################################################################

openmc_materials['Air'] = Material(name='Air')
openmc_materials['Air'].setDensity('g/cm3', 0.000616)
openmc_materials['Air'].addNuclide(nuclides['N-14'], 1.9681e-5)
openmc_materials['Air'].addNuclide(nuclides['N-15'], 7.1900e-8)
openmc_materials['Air'].addNuclide(nuclides['Ar-36'], 7.9414e-10)
openmc_materials['Air'].addNuclide(nuclides['Ar-38'], 1.4915e-10)
openmc_materials['Air'].addNuclide(nuclides['Ar-40'], 2.3506e-7)

# Includes O-17 and O-18 number densities
openmc_materials['Air'].addNuclide(nuclides['O-16'], 5.2993e-6)

# Unavailable in cross-section library?
#openmc_materials['Air'].addNuclide(nuclides['C-12'], 6.7565e-9)
#openmc_materials['Air'].addNuclide(nuclides['C-13'], 7.3706e-11)


################################################################################
#############################  Borosilicate Glass  #############################
################################################################################

openmc_materials['Boro. Glass'] = Material(name='Borosilicate Glass')
openmc_materials['Boro. Glass'].setDensity('g/cm3', 2.26000)
openmc_materials['Boro. Glass'].addNuclide(nuclides['B-10'], 9.6506e-4)
openmc_materials['Boro. Glass'].addNuclide(nuclides['B-11'], 3.9189e-3)
openmc_materials['Boro. Glass'].addNuclide(nuclides['Al-27'], 1.7352e-3)
openmc_materials['Boro. Glass'].addNuclide(nuclides['Si-28'], 1.6924e-2)
openmc_materials['Boro. Glass'].addNuclide(nuclides['Si-29'], 8.5977e-4)
openmc_materials['Boro. Glass'].addNuclide(nuclides['Si-30'], 5.6743e-4)

# Includes O-17 and O-18 number densities
openmc_materials['Boro. Glass'].addNuclide(nuclides['O-16'], 4.6624e-2)


################################################################################
##############################  Stainless Steel  ###############################
################################################################################

openmc_materials['Steel'] = Material(name='Stainless Steel')
openmc_materials['Steel'].setDensity('g/cm3', 8.0300)
openmc_materials['Steel'].addNuclide(nuclides['Si-28'], 9.5274e-4)
openmc_materials['Steel'].addNuclide(nuclides['Si-29'], 4.8400e-5)
openmc_materials['Steel'].addNuclide(nuclides['Si-30'], 3.1943e-5)
openmc_materials['Steel'].addNuclide(nuclides['Cr-50'], 7.6778e-4)
openmc_materials['Steel'].addNuclide(nuclides['Cr-52'], 1.4806e-2)
openmc_materials['Steel'].addNuclide(nuclides['Cr-53'], 1.6789e-3)
openmc_materials['Steel'].addNuclide(nuclides['Cr-54'], 4.1791e-4)
openmc_materials['Steel'].addNuclide(nuclides['Mn-55'], 1.7604e-3)
openmc_materials['Steel'].addNuclide(nuclides['Fe-54'], 3.4620e-3)
openmc_materials['Steel'].addNuclide(nuclides['Fe-56'], 5.4345e-2)
openmc_materials['Steel'].addNuclide(nuclides['Fe-57'], 1.2551e-3)
openmc_materials['Steel'].addNuclide(nuclides['Fe-58'], 1.6703e-4)
openmc_materials['Steel'].addNuclide(nuclides['Ni-58'], 5.6089e-3)
openmc_materials['Steel'].addNuclide(nuclides['Ni-60'], 2.1605e-3)
openmc_materials['Steel'].addNuclide(nuclides['Ni-61'], 9.3917e-5)
openmc_materials['Steel'].addNuclide(nuclides['Ni-62'], 2.9945e-4)
openmc_materials['Steel'].addNuclide(nuclides['Ni-64'], 7.6261e-5)


################################################################################
#############################  OpenCSG Materials  ##############################
################################################################################

# Store each OpenMC Material object as an OpenCSG Material object
for key in openmc_materials.keys():
  opencsg_materials[key] = get_opencsg_material(openmc_materials[key])
