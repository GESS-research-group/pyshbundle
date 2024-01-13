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

__author__ = """Amin Shakya"""
__email__ = 'aminshk50@gmail.com'
__version__ = '0.1.0'

# __init__.py with initialization code
print("Initializing PySHbundle v0.0.1")

# __init__.py with __all__
__all__ = ['basin_avg', 
           'clm2cs',
           'clm2sc',
           'eigengrav',
           'gaussian',
           'GRACE_Data_Driven_Correction_Vishwakarma',
           'GRACEpy',
           'grule',
           'gsha',
           'gshs',
           'iplm',
           'ispec',
           'klm2sc',
           'load_longterm_mean',
           'naninterp',
           'neumann',
           'new_io',
           'normalklm',
           'Phase_calc',
           'plm',
           'read_GRACE_SH_paths',
           'reader_replacer',
           'reader_replacer_csr',
           'reader_replacer_itsg',
           'reader_replacer_jpl',
           'sc2cs',
           'tws_cal',
           'visualisation_utils']

#Import individual modules
from .basin_avg import BasinAvg
from .clm2cs import clm2cs, clm2cs_new
from .clm2sc import clm2sc, clm2sc_new
# from .delta_sc import eigengrav
from .eigengrav import eigengrav
from .gaussian import Gaussian
from .GRACE_Data_Driven_Correction_Vishwakarma import deg_to_rad, GRACE_Data_Driven_Correction_Vishwakarma
# from .GRACEconstants import 
from .GRACEpy import upwcon, lovenr, lovenrPREM
from .grule import grule
from .gsha import gsha
from .gshs import GSHS
from .iplm import iplm
from .ispec import ispec
from .klm2sc import klm2sc #, klm2sc_new
from .load_longterm_mean import load_longterm_mean
from .naninterp import naninterp
from .neumann import neumann
from .new_io import clm2cs_new, read_jpl, parse_jpl_header, parse_jpl_data, parse_lines, read_csr, find_word, parse_csr_header, parse_csr_data, read_itsg, parse_itsg_header, parse_itsg_data, read_tn13, parse_tn13_header, parse_tn13_data, read_tn14, parse_tn14_header, parse_tn14_data, find_date_in_replacemnt_file, extract_C10_11_replcmnt_coeff, extract_C20_replcmnt_coeff, extract_C30_replcmnt_coeff, replace_zonal_coeff, klm2sc_new, sub2ind, cklm2sc_new, check_format
from .normalklm import normalklm
from .Phase_calc import PhaseCalc
from .plm import PLM, secrecur, lrecur, derivALF
# from .pyshbundle import 
from .read_GRACE_SH_paths import read_GRACE_SH_paths
from .reader_replacer import reader, TIME, last_4chars, reader_replacer
from .reader_replacer_csr import reader_replacer_csr 
from .reader_replacer_itsg import reader_replacer_itsg 
from .reader_replacer_jpl import reader_replacer_jpl
from .sc2cs import sc2cs
from .tws_cal import TWSCalc
# from .tws_py import 
from .visualisation_utils import sc_triplot, cs_sqplot, polar_plot, mapfield, ylm, ylm_plot, gshs_prepare