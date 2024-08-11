#!/usr/bin/env python

"""The setup script."""

import os
from setuptools import setup, find_packages

with open('README.md', encoding='utf8') as readme_file:
    readme = readme_file.read()

install_requires = ['pip', 'numpy', 'pandas','netCDF4', 'scipy',
                    'xarray', 'julian', 'scipy', 'geopandas',
                    'matplotlib', 'rasterio', 'shapely', 'tqdm','cartopy', 
                    'ipykernel', 'jupyterlab', 'rioxarray',],
setup(
    author="Amin Shakya",
    author_email='aminshk50@gmail.com',
    python_requires='>=3.9',
    classifiers=[
        'Development Status :: 3',
        'Intended Audience :: Scientific/Engineering',
        'License :: OSI Approved :: GNU General Public License v3 (GPLv3)',
        'Natural Language :: English',
        'Programming Language :: Python :: 3.9',
        'Programming Language :: Python :: 3.12',
    ],
    description="PySHbundle: A Python implementation of GRACE Spherical Harmonics Synthesis MATLAB codes SHbundle",
    package_data={"pyshbundle": ["data/*"]},
    install_requires=install_requires,
    license="GNU General Public License v3",
    long_description=readme,
    long_description_content_type='text/markdown',
    include_package_data=True,
    keywords='pyshbundle',
    name=['pyshbundle', 'GRACE', 'Spherical Harmonics', 'Synthesis Harmonics Synthesis', 'Spherical Harmonics Analysis'],
    packages=find_packages(include=['pyshbundle', 'pyshbundle.*']),
    setup_requires=[],
    test_suite='tests',
    tests_require=install_requires,
    url='https://github.com/mn5hk/pyshbundle',
    version='0.3.0',
    zip_safe=False,
)