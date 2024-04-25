# Installation

## Stable release

To install pyshbundle, run this command in your terminal:

This is the preferred method to install pyshbundle, as it will always install the most recent stable release.

If you don't have [pip](https://pip.pypa.io) installed, this [Python installation guide](http://docs.python-guide.org/en/latest/starting/installation/) can guide you through the process.

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

## From the source - for Devs/Contributors

Developers can access the latest development branch and 

```shell
# clone the repo and fetch the dev branch
$ git clone git@github.com:mn5hk/pyshbundle.git

# creating a new virtual environment
$ python3 -m venv <name-env>

# install the dependencies from the requirements-dev file
$ pip install -r ../pyshbundle/requirements-dev.txt

# activate the virtual environment environment
$ source </location-of-virt-env/name-env/bin/activate>

# install package into virtual environment
$ pip install ../pyshbundle/dist/<required-version>.tar.gz

# you also have the option to build the module using, be careful of 
$ python setup.py sdist
```
