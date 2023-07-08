# PySHbundle: A Python implementation of MATLAB codes SHbundle <br>

PySHBundle is the python implementation of the popular SHBundle toolbox originally written using MATLAB. 

TODO: Badges and and banner for the project


## Usage

1. Read and Load level-2 spherical harmonic data
2. create basin time series for TWS
3. perform grace data driven correction

## 1. How to install <br>
### 1.1 For Users
We recommend using Mamba to install required packages <br>
`conda create pysh` <br>
`conda activate pysh` <br>
`conda install -c conda-forge mamba` <br>
`mamba install -c conda-forge numpy pandas netCDF4 scipy xarray julian scipy geopandas matplotlib rasterio salem shapely` <br>

### 1.2 For Contributors
1. Fork the repository and clone it onto your system
2. In order to use the codebase use 
   
   ```
   Base repo path - ../open-source/pyshbundle

   pyshbundle (base repo)
    | - pyshbundle (actual module)
        | - all the codes
    | - notebooks
    | - docs
    | - and rest other folders
   ```



## Read the Docs

Please find the docs here - [link goes here]


## Contribution Guidelines
1. Found a bug?
 - Please check the know issues section before reporting
 - If the bug not mentioned in known issue section, report the bug in the given format
   
2. Want a new feature to be added?
 - Try creating a feature request in the GitHub issues section

3. Want to start the open source journey?
 - Try helping out in improving the documentation for the project
 - Try the module and come up with innovative tutorial ideas and implement them


### Reporting bugs/issues
> Expected format
[Update]

## Known Issues

[To Be Updated...]


# License Statement

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
    Please note that PySHbundle has adapted the following code packages, 
    both licensed under GNU General Public License
  1. SHbundle: https://www.gis.uni-stuttgart.de/en/research/downloads/shbundle/ 
  2. Downscaling GRACE Total Water Storage Change using Partial Least Squares Regression
    https://springernature.figshare.com/collections/Downscaling_GRACE_Total_Water_Storage_Change_using_Partial_Least_Squares_Regression/5054564 


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



![](https://visitor-badge.glitch.me/badge?page_id=mn5hk.mat2py) <br>
[Geodesy for Earth system science (GESS) research Group at ICWaR, IISc](https://ultra-pluto-7f6d1.netlify.app/)
