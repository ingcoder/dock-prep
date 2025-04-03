from setuptools import setup, Extension, find_packages

c_extension = Extension("dock_prep.c_extension_sample", sources=['dock_prep/c_extension_sample.c'])

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
        'pdb2pqr'
    ],
    extras_require={
        'dev': [
            'pytest>=6.0.0',
            'pytest-cov>=2.10.0',
            'black>=21.5b2',
            'flake8>=3.9.0',
            'mypy>=0.812',
            'sphinx>=4.0.0',
            'sphinx-rtd-theme>=0.5.0',
        ],
    },
    entry_points={
        'console_scripts': [
            'dock-prep=dock_prep.run:main',
            'dock-prep-check=dock_prep.check_dependencies:check_all_dependencies',
        ],
    },
    python_requires='>=3.6',
    classifiers=[
        'Development Status :: 3 - Alpha',
        'Intended Audience :: Science/Research',
        'Topic :: Scientific/Engineering :: Bio-Informatics',
        'Programming Language :: Python :: 3',
        'Programming Language :: C',
        'License :: OSI Approved :: MIT License',
        'Operating System :: OS Independent',
    ],
    project_urls={
        'Bug Reports': 'https://github.com/ingcoder/pdbqt-converter/issues',
        'Source': 'https://github.com/ingcoder/pdbqt-converter',
    },
)