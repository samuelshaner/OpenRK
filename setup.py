__author__ = 'Will Boyd'
__email__ = 'wboyd@mit.edu'


from distutils.core import setup

setup(name='InferMC',
      version='0.1',
      description='A package for inference of Monte Carlo tally data.',
      author='Will Boyd',
      author_email='wboyd@mit.edu',
      packages=['datasets', 'datasets.BEAVRS', 'datasets.BEAVRS.materials',
                'datasets.BEAVRS.geometry', 'datasets.BEAVRS.geometry.openmoc',
                'datasets.BEAVRS.geometry.openmoc.assemblies',
                'datasets.BEAVRS.geometry.openmoc.colorsets']
     )