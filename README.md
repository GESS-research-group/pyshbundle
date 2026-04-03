# PySHbundle: A Python implementation of GRACE Spherical Harmonics Synthesis MATLAB codes SHbundle <br>

![](https://visitor-badge.glitch.me/badge?page_id=GESS-research-group.pyshbundle)
[![Build and Test](https://github.com/GESS-research-group/pyshbundle/actions/workflows/python-package-conda.yml/badge.svg)](https://github.com/GESS-research-group/pyshbundle/actions/workflows/python-package-conda.yml)
[![License: GPL v3](https://img.shields.io/badge/License-GPLv3-blue.svg)](https://www.gnu.org/licenses/gpl-3.0)
[![Python 3.9 | 3.12](https://img.shields.io/badge/python-3.9%20%7C%203.12-blue)](https://www.python.org)
[![Documentation](https://img.shields.io/badge/docs-github.io-green)](https://gess-research-group.github.io/pyshbundle/) <br>

This package, `PySHbundle` provides tools to process GRACE data, such as, the computation of anomalies, substitution of poor quality low degree coefficients, reducing noise in GRACE data using filtering approaches, signal leakage correction using `GDDC`, etc. In addition, the package provides a flexibility for future development and addition of further processing choices for handling GRACE data for hydrological application.

PySHBundle is a tool to process GRACE L2 data and re-implements the popular [SHBundle](https://www.gis.uni-stuttgart.de/en/research/downloads/shbundle/) and [DataDrivenCorrection Bundle](https://www.gis.uni-stuttgart.de/en/research/downloads/datadrivencorrectionbundle/) tools originally written using MATLAB. 


## Usage
1. Read and Load level-2 spherical harmonic data
2. Create basin time series for TWS
3. Perform grace data driven correction
4. Plot spherical harmonic related plots

## 1. How to install <br>
### 1.1 For Users
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

### 1.2 For Devs/Contributors
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

## Trying it out

Data for trying out this new tool is included in the repo. After installing and cloning the repo, go to the notebooks directory in order to find explainatory ipython jupyter notebooks. Simply activate the virtual environment and fire up these jupyter notebooks. Available notebooks:

1. Introduction to Spherical Harmonics
2. Loading the data
3. Visualizations
4. Terrestrial Water Storage (TWS) Time Series
5. Tests and Validation notebook

## Docs

Please find the docs here - [PySHBundle](https://gess-research-group.github.io/pyshbundle/)

## Testing

The test suite validates the accuracy of the TWS computation against a MATLAB reference solution.

### Framework
Tests are written in `pytest` and live in the `tests/` directory.

### Running the tests
```shell
# From the project root
pytest tests/ -v
```

By default, the tests look for GRACE input data in `data/JPL_input/`. You can override this with an environment variable:
```shell
PYSHBUNDLE_DATA_DIR=/path/to/your/data pytest tests/ -v
```

### What is tested
The suite runs 6 tests using 60 months of JPL GRACE RL06 data compared against a MATLAB-generated reference TWS field (`tws_sh.mat`):

| Test | Description | Threshold |
|---|---|---|
| `test_tws_output_shape` | Computed and reference arrays have identical shape | — |
| `test_tws_output_dtype` | Output is `float32` | — |
| `test_gridwise_rmse` | Gridwise RMSE (mm) between computed and reference TWS | < 1e-3 |
| `test_gridwise_nrmse` | Gridwise NRMSE (normalised by reference std) | < 1e-5 |
| `test_no_nan_in_output` | No NaN values in computed TWS | — |
| `test_no_nan_in_reference` | No NaN values in MATLAB reference (sanity check) | — |

### Continuous Integration
Tests run automatically on every push via GitHub Actions across 6 combinations:

- **OS:** Ubuntu, macOS, Windows
- **Python:** 3.9, 3.12

## Contributing

Contributions are welcome, and they are greatly appreciated! Every
little bit helps, and credit will always be given.

You can contribute in many ways:

## Types of Contributions

### Report Bugs

Report bugs at [GitHub Issues](https://github.com/GESS-research-group/pyshbundle/issues)


If you are reporting a bug, please include:

-   Your operating system name and version.
-   Any details about your local setup that might be helpful in troubleshooting.
-   Detailed steps to reproduce the bug.

### Fix Bugs

Look through the GitHub issues for bugs. Anything tagged with `bug` and
`help wanted` is open to whoever wants to implement it.

### Implement Features

Look through the GitHub issues for features. Anything tagged with
`enhancement` and `help wanted` is open to whoever wants to implement it.

### Write Documentation

pyshbundle could always use more documentation,
whether as part of the official pyshbundle docs,
in docstrings, or even on the web in blog posts, articles, and such.

### Submit Feedback

The best way to send feedback is to file an issue at
[GitHub Issues](https://github.com/GESS-research-group/pyshbundle/issues)

If you are proposing a feature:

-   Explain in detail how it would work.
-   Keep the scope as narrow as possible, to make it easier to implement.
-   Remember that this is a volunteer-driven project, and that contributions are welcome :).

## License Statement

This file is part of PySHbundle. <br>
    PySHbundle is free software: you can redistribute it and/or modify<br>
    it under the terms of the GNU General Public License as published by<br>
    the Free Software Foundation, either version 3 of the License, or<br>
    (at your option) any later version.<br>
<br>
    This program is distributed in the hope that it will be useful,<br>
    but WITHOUT ANY WARRANTY; without even the implied warranty of<br>
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the<br>
    GNU General Public License for more details.<br>
<br>
    You should have received a copy of the GNU General Public License<br>
    along with this program.  If not, see <http://www.gnu.org/licenses/>.<br>
    

## Acknowledgement:
Please note that PySHbundle has adapted the following code packages,both licensed under GNU General Public License

  1. SHbundle: https://www.gis.uni-stuttgart.de/en/research/downloads/shbundle/ 
  2. Downscaling GRACE Total Water Storage Change using Partial Least Squares Regression: https://springernature.figshare.com/collections/downscaling_GRACE_Total_Water_Storage_Change_using_Partial_Least_Squares_Regression/5054564 


## Key Papers Referred:
 1. Vishwakarma, B. D., Horwath, M., Devaraju, B., Groh, A., & Sneeuw, N. (2017). 
    A data‐driven approach for repairing the hydrological catchment signal damage 
    due to filtering of GRACE products. Water Resources Research, 
    53(11), 9824-9844. https://doi.org/10.1002/2017WR021150 

 2. Vishwakarma, B. D., Zhang, J., & Sneeuw, N. (2021). 
    Downscaling GRACE total water storage change using 
    partial least squares regression. Scientific data, 8(1), 95.
    https://doi.org/10.1038/s41597-021-00862-6 


## How to Cite?
Coming soon!


## Follow the Research Group
<a href="https://ultra-pluto-7f6d1.netlify.app" rel="Geodesy for Earth system science (GESS) research Group at ICWaR, IISc">
    <img src="./examples/imgs/logoGESS.jpg" alt="Geodesy for Earth system science (GESS) research Group at ICWaR, IISc" width="300" height="300">
</a>
