# Welcome to PySHbundle

[![image](https://img.shields.io/pypi/v/pyshbundle.svg)](https://pypi.python.org/pypi/pyshbundle)


**PySHbundle: A Python implementation of MATLAB codes SHbundle**

Our package, `PySHbundle` has been developed in a modularized manner. The package provides tools to process GRACE data, such as, the computation of anomalies, substitution of poor quality low degree coefficients, reducing noise in GRACE data using filtering approaches, signal leakage correction using `GDDC`, etc. In addition, the package provides a flexibility for future development and addition of further processing choices for handling GRACE data for hydrological application.<br>

It is hoped the contribution will make GRACE L2 data processing more accessible to a wider audience of young researchers. Our python package is titled `PySHbundle` and the working code can be accessed at [GitHub](https://github.com/mn5hk/pyshbundle). <br>


## How to install <br>
```shell
# creating a new virtual environment
$ python3 -m venv <name-env>
# activate the virtual environment environment
$ source </location-of-virt-env/name-env/bin/activate>
# install package into virtual environment
$ pip install pyshbundle

# clone the repository in order to access the notebooks and data
$ git clone git@github.com:mn5hk/pyshbundle.git
```
For more details refer to the installation section.


## Group

[Geodesy for Earth system science (GESS) research Group at ICWaR, IISc](https://ultra-pluto-7f6d1.netlify.app/)





### Credits for the theme

This package was created with [Cookiecutter](https://github.com/cookiecutter/cookiecutter) and the [giswqs/pypackage](https://github.com/giswqs/pypackage) project template.
