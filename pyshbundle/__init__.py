"""Top-level package for pyshbundle.

This file is part of PySHbundle. 
    PySHbundle is free software: you can redistribute it and/or modify
    it under the terms of the GNU General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.

    This program is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
    GNU General Public License for more details.

    You should have received a copy of the GNU General Public License
    along with this program.  If not, see <http://www.gnu.org/licenses/>.
    
- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 

Acknowledgement Statement:
    Please note that PySHbundle has adapted the following code packages, 
    both licensed under GNU General Public License
    1. SHbundle: https://www.gis.uni-stuttgart.de/en/research/downloads/shbundle/

    2. Downscaling GRACE Total Water Storage Change using 
    Partial Least Squares Regression
    https://springernature.figshare.com/collections/Downscaling_GRACE_Total_Water_Storage_Change_using_Partial_Least_Squares_Regression/5054564 
    
Key Papers Referred:
    1. Vishwakarma, B. D., Horwath, M., Devaraju, B., Groh, A., & Sneeuw, N. (2017). 
    A data‐driven approach for repairing the hydrological catchment signal damage 
    due to filtering of GRACE products. Water Resources Research, 
    53(11), 9824-9844. https://doi.org/10.1002/2017WR021150

    2. Vishwakarma, B. D., Zhang, J., & Sneeuw, N. (2021). 
    Downscaling GRACE total water storage change using 
    partial least squares regression. Scientific data, 8(1), 95.
    https://doi.org/10.1038/s41597-021-00862-6 
    """
# 2024-06-10, updated: Vivek Kumar Yadav, IISc Bengaluru

__author__ = """Vivek Kumar Yadav"""
__email__ = 'viveky@iisc.ac.in'
__version__ = '0.3.0'

# __init__.py with initialization code
print("Initializing PySHbundle v0.3.0")

# __init__.py with __all__
__all__ = ['GRACEpy',
           'io',
           'viz_utils',
           'reshape_SH_coefficients',
           'shutils',
           'hydro',
           'pysh_core'
           ]


from .GRACEpy import upwcon, lovenr, lovenrPREM
from .io import extract_SH_data, extract_deg1_coeff_tn13, extract_deg2_3_coeff_tn14, \
    parse_lines, \
    parse_jpl_file, parse_csr_file, parse_itsg_file, parse_jpl_header, parse_csr_header, parse_itsg_header, \
    parse_tn13_header, parse_tn14_header, \
    find_date_in_replacemnt_file, extract_C10_11_replcmnt_coeff, extract_C20_replcmnt_coeff, \
    extract_C30_replcmnt_coeff, \
    read_GRACE_SH_paths, load_longterm_mean
from .sc2cs import sc2cs
from .reshape_SH_coefficients import sc2cs, clm2cs, clm2sc, cs2sc, klm2sc, cklm2sc_new
from .hydro import TWSCalc, area_weighting, Basinaverage
from .shutils import plm, iplm, ispec, eigengrav, grule, Gaussian, neumann, naninterp, normalklm
from .pysh_core import gshs, gsha, GRACE_Data_Driven_Correction_Vishwakarma, PhaseCalc
from .viz_utils import sc_triplot, cs_sqplot, polar_plot, mapfield, ylm, ylm_plot, gshs_prepare