#!/usr/bin/env python

"""The setup script."""

from setuptools import setup, find_packages

with open('README.md', encoding='utf8') as readme_file:
    readme = readme_file.read()

install_requires = ['pip', 'numpy', 'pandas','netCDF4', 'scipy',
                    'xarray', 'julian', 'geopandas',
                    'matplotlib', 'rasterio', 'shapely', 'tqdm','cartopy', 
                    'ipykernel', 'jupyterlab', 'rioxarray',],
setup(
    name='pyshbundle',
    version='0.3.0',
    python_requires='>=3.9',
    packages=find_packages(include=['pyshbundle', 'pyshbundle.*']),
    description="PySHbundle: A Python implementation of GRACE Spherical Harmonics Synthesis MATLAB codes SHbundle",
    license="GNU General Public License v3",
    author="Vivek Kumar Yadav",
    author_email='viveky@iisc.ac.in',
    url='https://github.com/mn5hk/pyshbundle',
    install_requires=install_requires,
    long_description=readme,
    long_description_content_type='text/markdown',
    package_data={'pyshbundle': ['data/*/*',],},
    include_package_data=True,
    keywords=['pyshbundle', 'GRACE', 'Spherical Harmonics', 'Synthesis Harmonics Synthesis', 'Spherical Harmonics Analysis'],
    test_suite='tests',
    tests_require=install_requires,
    zip_safe=False,
    classifiers=[
        'Development Status :: 3',
        'Intended Audience :: Scientific/Engineering',
        'License :: OSI Approved :: GNU General Public License v3 (GPLv3)',
        'Natural Language :: English',
        'Programming Language :: Python :: 3.9',
        'Programming Language :: Python :: 3.12',
    ],
)