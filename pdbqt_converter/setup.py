from setuptools import setup, Extension, find_packages

c_extension = Extension("pdbqt_converter.c_tutorial", sources=['pdbqt_converter/c_tutorial.c'])

setup(
    name="pdbqt_converter",
    version='0.1.0-alpha',
    description='A tool for converting PDB files to PDBQT format for molecular docking',
    author='Ingrid Barbosa-Farias',
    author_email='ingrid@simatomic.com',
    packages=find_packages(),
    ext_modules=[c_extension],
    install_requires=[
        'numpy',
        'biopython',
        'pdbfixer',
        'openmm',
    ],
    entry_points={
        'console_scripts': [
            'pdbqt-converter=pdbqt_converter.run:main',
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