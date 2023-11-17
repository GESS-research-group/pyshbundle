# PySHbundle: A Python implementation of MATLAB codes SHbundle <br>

![](https://visitor-badge.glitch.me/badge?page_id=mn5hk.mat2py) <br>


PySHBundle is the python implementation of the popular SHBundle toolbox originally written using MATLAB. 

TODO: Badges and and banner for the project


## Usage

1. Read and Load level-2 spherical harmonic data
2. Create basin time series for TWS
3. Perform grace data driven correction
4. Plot spherical harmonic related plots

## 1. How to install <br>
### 1.1 For Users
Currently the package is not yet finalized hence the version on PyPI is outdated and should not be used (as of now). Please follow the steps mentioned in the following section till things are finalized. <br>

### 1.2 Till things get finalized

1.  Fork the pyshbundle repo on GitHub.

2.  Clone your fork locally:

    ```shell
    $ git clone git@github.com:your_name_here/pyshbundle.git
    ```

3.  Create a new virtual environment or conda environment to install all the necessary dependencies. using `conda` is recommended along with `jupyter lab`. Use of python 3.x is recommended

    ```shell
    $ conda create -n pyshbundle-env
    $ conda activate pyshbundle-env
    $ conda install -c conda-forge --file requirements_dev.txt -y
    ```
4. Note that the base path to the entire repo is important while importing (this is tempoary only, after PyPi module gets updated this approach will not be required) 
   ```
   Example base repo path -> ../open_source/pyshbundle

   pyshbundle (base container repo)
    | - pyshbundle (codes and functions reside here)
        | - all the codes
    | - notebooks
    | - docs
    | - and rest of the other folders
   ```

5. In order to use any of the functions change the current directory to the pyshbundle (base repo) then use import as usual
   ```python
   import os
   os.chdir(../open_source/pyshbundle)

   import pyshbundle
   # or import individual functions
   from pyshbundle import read_jpl
   ```
   after importing the fucntions this way you are all good to go.


## Read the Docs

Please find the docs here - [PySHBundle](https://abhimhamane.github.io/pyshbundle/)


## Contributing

Contributions are welcome, and they are greatly appreciated! Every
little bit helps, and credit will always be given.

You can contribute in many ways:

## Types of Contributions

### Report Bugs

Report bugs at [GitHub Issues](https://github.com/mn5hk/pyshbundle/issues)


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
[GitHub Issues](https://github.com/mn5hk/pyshbundle/issues)

If you are proposing a feature:

-   Explain in detail how it would work.
-   Keep the scope as narrow as possible, to make it easier to implement.
-   Remember that this is a volunteer-driven project, and that contributions are welcome :)

## Known Issues

1. New implementation of reading data has some bug. It works for CSR but not properly properly for JPL and ITSG sources. For more information check the current issues.
2. There is some inconsistency related to importing functions it is recommended to always use `import pyshbundle` along with any of the targeted imports you want.


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
TODO: Add citing info







## Follow the Research Group

<a href="https://ultra-pluto-7f6d1.netlify.app" rel="Geodesy for Earth system science (GESS) research Group at ICWaR, IISc">![Geodesy for Earth system science (GESS) research Group at ICWaR, IISc](../notebooks/imgs/logoGESS.jpg)</a>
