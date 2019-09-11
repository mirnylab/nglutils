from setuptools import setup, find_packages

with open('README.md') as f:
  readme = f.read()

setup(
    name='nglutils',
    version='0.0.1',
    description='NGL utils for trajectory visualization',
    long_description=readme,
    author='Simon Grosse-Holz',
    url='https://github.com/mirnylab/nglutils',
    packages=find_packages(exclude=('tests', 'docs')),
    install_requires=[
        'matplotlib', 'numpy', 'mdtraj', 'openmmlib', 'openmm', 'nglview>=2.1.0,<=2.1.1'
        # Jupyter and widgets are currently not included in the dependencies
    ])

# For installation of the nglview widget for ipython one might need: 
# $jupyter-nbextension enable --py --user widgetsnbextension
# $jupyter-nbextension enable --py --user nglview