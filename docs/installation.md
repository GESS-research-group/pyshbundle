# Installation

## For Users

The recommended installation method is to clone the repository and install locally. This also gives you access to the example notebooks and data.

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

> **Note:** The package is available on PyPI but is currently broken.
> Please avoid installing via `pip install pyshbundle` until this is resolved.

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
