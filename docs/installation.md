# Installation

## For Users

The simplest way to install is via pip:

```shell
$ pip install pyshbundle
```

To also access the example notebooks and data, clone the repository:

```shell
# Clone the repository
$ git clone https://github.com/GESS-research-group/pyshbundle.git
$ cd pyshbundle

# Create and activate a virtual environment
$ python3 -m venv <name-env>
$ source <name-env>/bin/activate  # On Windows: <name-env>\Scripts\activate

# Install the package
$ pip install .
```

## From Source — For Devs/Contributors

```shell
# Clone the repository
$ git clone https://github.com/GESS-research-group/pyshbundle.git
$ cd pyshbundle

# Create and activate a virtual environment
$ python3 -m venv <name-env>
$ source <name-env>/bin/activate  # On Windows: <name-env>\Scripts\activate

# Install the package in editable mode with dev dependencies
$ pip install -r requirements-dev.txt
$ pip install -e .

# To build a source distribution
$ python -m build
```
