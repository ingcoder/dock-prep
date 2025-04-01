from setuptools import setup, Extension, find_packages

c_extension = Extension("dock_prep.c_tutorial", sources=['dock_prep/c_tutorial.c'])

setup(
    name="dock_prep",
    version='0.1.0-alpha',
    description='A tool for preparing PDB files for molecular docking',
    author='Ingrid Barbosa-Farias',
    author_email='ingrid@simatomic.com',
    packages=find_packages(),
    ext_modules=[c_extension],
    install_requires=[
        'numpy',
        'biopython',
        'pdbfixer',
        'openmm',
        'pdb2pqr',
        'openbabel'
    ],
    entry_points={
        'console_scripts': [
            'dock-prep=dock_prep.run:main',
        ],
    },
    python_requires='>=3.6',
    classifiers=[
        'Development Status :: 3 - Alpha',
        'Intended Audience :: Science/Research',
        'Topic :: Scientific/Engineering :: Bio-Informatics',
        'Programming Language :: Python :: 3',
        'Programming Language :: C',
    ],
)