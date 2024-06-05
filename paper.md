---
title: 'PySHbundle: A Python software to convert GRACE Spherical Harmonic Coefficients to gridded mass change fields'
tags:
  - Python
  - GRACE
  - Spherical Harmonic Analysis
  - Spherical Harmonic Synthesis
  - GRACE Data Driven Correction
authors:
  - name: Amin Shakya
    orcid: 0000-0002-4706-826X
    equal-contrib: true # (This is how you can denote equal contributions between multiple authors)
    corresponding: true # (This is how to denote the corresponding author)
    affiliation: "1,2"
  - name: Vivek Kumar Yadav
    orcid: 0009-0000-7156-4450
    equal-contrib: true # (This is how you can denote equal contributions between multiple authors)
    affiliation: 2
  - name: Abhishek Mhamane
    orcid: 0000-0001-9788-0371
    equal-contrib: true # (This is how you can denote equal contributions between multiple authors)
    affiliation: 3
  - name: Tsungrojungla Walling
    orcid: 0009-0006-9323-1191
    affiliation: 4
  - name: Shard Chander
    affiliation: 5
  - name: Bhaskar R. Nikam
    affiliation: 6
  - name: Nagesh Kumar Dasika
    orcid: 0009-0006-9323-1191
    affiliation: 7
  - name: Bramha Dutt Vishwakarma
    orcid: 0000-0003-4787-8470
    affiliation: "2,8" # (Multiple affiliations must be quoted)
affiliations:
 - name: Faculty of Geo-Information Science and Earth Observation, University of Twente, the Netherlands
   index: 1
 - name: Interdisciplinary Centre for Water Research, Indian Institute of Science, India
   index: 2
 - name: National Centre for Geodesy, Indian Institute of Technology Kanpur, India
   index: 3
 - name: Undergraduate Programme, Indian Institute of Science, India
   index: 4
 - name: Land Hydrology Division, Space Applications Centre, Indian Space Research Organisation, India 
   index: 5
 - name: Earth Observation Applications & Disaster Management Support Programme Office (EDPO), Indian Space Research Organisation, India
   index: 6
 - name: Department of Civil Engineering, Indian Institute of Science, India
   index: 7
 - name: Centre of Earth Science, Indian Institute of Science, India
   index: 8

date: 15 November 2023
bibliography: paper.bib

# Optional fields if submitting to a AAS journal too, see this blog post:
# https://blog.joss.theoj.org/2018/12/a-new-collaboration-with-aas-publishing
# aas-doi: 
# aas-journal: 
# IMPORTANA NOTE: do not enclose the /begin{equation} tag withing $$, this leads to 'latex math ennvironment error' causing issues with compilation of manuscript. Either use $$ or /begin{equation} syntax for a equation block, and $ for inline equation.
---

# Summary

`GRACE` (Gravity Recovery and Climate Experiment) satellite mission has been mapping mass changes near the surface of the Earth since 2002. Since mass redistribution at shorter temporal scales is dominated by hydrology, GRACE has transformed our understanding of changes in the hydrosphere. GRACE data has been used for monitoring and studying groundwater depletion, floods, droughts, etc. GRACE satellite products are typically released at various levels of complexity, often referred to as processing levels. Level 1 is the satellite instrument data that is processed to obtain Level 2 data: the spherical harmonic coefficients representing the mean monthly gravity field of the Earth. Level 2 are further processed to obtain level 3 products; gridded mass change estimates expressed as terrestrial water storage anomalies (`TWSA`). Level 2 data are unconstrained gravity field solutions and are noisy, which are filtered and corrected for known artifacts and signals from solid Earth processes to obtain level 3 products that are useful for hydrology. Processing choices, such as filter properties and type, have a significant impact on the accuracy and the resolution of final gridded output. Therefore, level 3 users must be cautious when using GRACE data for specific applications. Since, majority of the GRACE data user community is not well versed with level 2 data processing, they often use off the shelf product with doubts on the efficacy of GRACE mission. Here we developed an open-source GRACE level 2 processing toolbox to provide users with more control over processing choices. A python module, called PySHbundle, is developed that converts GRACE level 2 (`L2`) Spherical Harmonics data products to Level 3 (`L3`) `TWSA` products while applying the data-driven correction algorithm for reducing the impact of filtering on signal. With this contribution, we hope to enable further usage of GRACE data for Earth system science.

# Introduction

GRACE stands for the Gravity Recovery and Climate Experiment, a joint satellite mission by the National Aeronautics and Space Administration (`NASA`) in the USA, and the Geoforschung Zentrum (`GFZ`) in Germany. GRACE mission launced on 17 March, 2002, and ended on 27 October, 2017. Some details of the GRACE mission are provided in `Table 1`. GRACE has a successor, GRACE-FO, which was successfully launched on 22 May 2018 and is currently operational. GRACE consists of two identical satellites in the same orbit separated by approx. 220 km. The mission measures changes in the intersatellite distance with a microwave ranging system that gives an accuracy in the range of micrometers [@wahr1998time]. When the satellite system comes in the vicinity of a temporal mass anomaly, the relative intersatellite distance changes and it can be inverted to estimate the mass change near the surface of the Earth. Over the continental land surface, the hydrological processes are the major driver of the variation in mass anomaly at monthly to decadal scales. However various other signals such as oceanic and atmospheric variations, high frequency tidal mass changes, systemic correlated errors, etc. are also part of the obtained GRACE signals [@humphrey2023using]. Some of the unwanted signals, such as the high frequency tidal mass changes in the ocean and the atmosphere, are modelled and removed at level 1 processing [@flechtner2007aod1b], while noise is still present at level 2 and it requires filtering [@wahr1998time; @vishwakarma2017understanding]. The choice of filter and/or subsequent steps to counter the signal loss due to filtering have an impact on the quality of GRACE products that are of interest to hydrologists [@humphrey2023using; @vishwakarma2020monitoring]. The estimated hydrological signal is represented in terms of  `total water storage anomaly` (`TWSA`), which is the change in the water mass over a vertical column. Conventionally, it is represented in terms of `equivalent water height` (`m`). <br>

<i>Table 1: Summary of GRACE satellite mission [[source]](https://www2.csr.utexas.edu/grace/mission/mdetail.html)</i>

| Parameter        |    Details     | 
| --------------   |:--------------:| 
| Start of Mission | 17 March 2002  | 
| End of Mission   | 27 October 2017| 
| Altitude         | 485 km         |
| Intersatellite distance | ~ 220 km 
| Inclination      | 89.0°          | 
| Period           | 94.5 minutes   |  

Various research centers provide GRACE data, such as the University of Texas Center for Space Research (`CSR`), Jet Propulsion Laboratory (`JPL`), and the German Research Center for Geosciences (`GFZ`) and so on. `Level 2` data products are the spherical harmonic coefficients representing the mean monthly gravity field of the Earth. `Level 3` consists of gridded mass anomalies or other standardized products, such as the Monthly Ocean/Land Water Equivalent Thickness Surface-Mass Anomaly. Obtaining `level 3` products from `level 2` requires filtering that affects the signal quality and resolution, which is why another suite of products called mass concentration blocks or `mascons` are also available that are designed to take care of signal degradation due to filtering inherently. More details on the mascon approaches for studying gravity fields and the approaches used by the different data centers for generating mascon products may be referred to in [@antoni2022review]. Mascon products from various centres differ due to the difference in the post-processing strategy used by these centres, which emphasizes that the processing choices have an impact of the gridded GRACE data. Processing `Level 2` data gives the user the freedom and the flexibility to explore GRACE data for a specific application over a specific region as per their convenience and belief. <br>

Converting spherical harmonic coefficients (`level 2`) to gridded field (`level 3`) is called spherical harmonic synthesis and vice-versa is spherical harmonic analysis. `Level 3` products may further be processed to obtain region-averaged timeseries data, labelled as `Level 4` products. Various tools exist to process GRACE data and to analyze it. Some of these are developed in the `MATLAB` programming language: [`SHbundle`](https://www.gis.uni-stuttgart.de/en/research/downloads/shbundle) [@SHbundle], [`GRACE Data Driven Correction`](https://www.gis.uni-stuttgart.de/en/research/downloads/datadrivencorrectionbundle) [@vishwakarma2017understanding],  [`LUH-GRACE2018`](https://www.ife.uni-hannover.de/en/services/luh-grace) [@koch2020luh], [`GRAMAT`](https://link.springer.com/article/10.1007/s12145-018-0368-0) [@feng2019gramat], [`SHADE`](https://www.sciencedirect.com/science/article/pii/S0098300418302760) [@piretzidis2018shade], [`GRACETOOLS`](https://www.mdpi.com/2076-3263/8/9/350) [@darbeheshti2018gracetools], [`SSAS GRACE filter`](https://github.com/shuang-yi/SSAS-GRACE-filter)[@yi2022novel], etc. Similarly, some GRACE data processing tools are also available based on the python programming language. These include [`gravity-toolkit`](https://gravity-toolkit.readthedocs.io/en/latest/) [@gravity-toolkit], [`ggtools`](https://pypi.org/project/ggtools/1.1.0/) [@ggtools] and [`GRACE-filter`](https://github.com/strawpants) [@GRACEfilter]. General tools for spheric harmonic analysis are also available, such as [`SHTools`](https://agupubs.onlinelibrary.wiley.com/doi/full/10.1029/2018GC007529) [@wieczorek2018shtools]. [`SHbundle`](https://www.gis.uni-stuttgart.de/en/research/downloads/shbundle) provide MATLAB scripts for `spheric harmonic synthesis` and `spheric harmonic analysis`. The first version of the code was developed in 1994 while the latest version with upgrades can be found dated 2018. [`GRAMAT`](https://link.springer.com/article/10.1007/s12145-018-0368-0) provides a similar MATLAB-based scripts for processing GRACE spherical harmonics data to obtain spatiotemporal global mass variations. The [`GRAMAT`](https://link.springer.com/article/10.1007/s12145-018-0368-0) toolbox includes Gaussian smoothening filter to reduce noise that appears strongly as North-South stripes, spherical harmonic analysis and synthesis routines, signal leakage reduction routines, harmonic analysis of times series over regions, and uncertainty analysis of GRACE estimates [@feng2019gramat]. [`SHADE`](https://www.sciencedirect.com/science/article/pii/S0098300418302760) provides a MATLAB-based toolbox for the empirical de-correlation of GRACE monthly spherical harmonics [@piretzidis2018shade]. [`Gravity-toolkit`](https://gravity-toolkit.readthedocs.io/en/latest/) is a python-based package meant to handle GRACE L2 data products. Its functionalities include visualization of GRACE and GRACE-FO L2 data products, and the estimation of GRACE and GRACE-FO L2 data product errors. [`Gg-tools`](https://pypi.org/project/ggtools/1.1.0/) too contain similar tools for signal correction and for conversion of GRACE L2 products to L3. `GRACE-filter` provides tool for filtering of GRACE L2 product using DDK filter based on [@kusche2009decorrelated].
 
# Statement of need

As discussed in the previous section, GRACE level 2 products can be complicated to process to obtain level 3 products, such as the Total Water Storage Anomaly (TWSA). While level 3 mascon products exist, this software intends to lower the barrier for researchers and end users to use the level 2 products, with the ease of selecting different processing choices.

Furthermore, various research groups provide the level 2 spherical harmonics products. This software has been developed to be able to process, three of the widely used level 2 products, from:  the University of Texas Center for Space Research (CSR), Jet Propulsion Laboratory (JPL), and the German Research Center for Geosciences (GFZ). The software is further developed in modular form. As such, it can be extended to handle Level 2 products from other providers, if needed.

Thirdly, the software has been developed, closely in line with a popular GRACE spherical harmonics processing Matlab code bundle, called [`SHbundle`](https://www.gis.uni-stuttgart.de/en/research/downloads/shbundle) [@SHbundle] and the `GRACE Data Driven Correction (GDDC)` [@vishwakarma2017understanding]. Both the code bundles are written in Matlab programming language, and distributed under the GNU license. While `PySHbundle` makes the spherical harmonics processing codes available in Python, the similarity in its architecture to `SHbundle` and `GDDC` makes it easy for researchers from different programming language preferences to understand each other's processing steps.<br>

Furthermore, our package, `PySHbundle` has been developed in a modularized manner. The package provides tools to process GRACE data, such as, the computation of anomalies, substitution of poor quality low degree coefficients, reducing noise in GRACE data using filtering approaches, signal leakage correction using `GDDC`, etc. In addition, the package provides a flexibility for future development and addition of further processing choices for handling GRACE data for hydrological application.

Lastly, with modularization, use of open source programming language (Python) and distribution under the GNU license, we hope to make GRACE level 2 data processing more accessible to a wider audience, both by reducing the technical barrier due to processing complexities, as well as, accessible globally alleviating financial barriers, adhering to the [FAIR principle](https://www.go-fair.org/fair-principles/). The software can be used by researchers as well as young students. We anticipate the software's use in short courses, such as the [GRACE Hackweek](https://www.quantumfrontiers.de/de/aktuelles/veranstaltungen/details/news/grace-hackweek-3) organized by IIT Kanpur, India, where many young researchers from low-income countries participate.

Our python package is titled `PySHbundle` and the working code can be accessed at [GitHub](https://github.com/mn5hk/pyshbundle).

# Implementation

Obtaining gridded fields from GRACE spherical harmonic coefficients consists of several steps, including obtainnig the spherical harmonic coefficiets from the data providers, replacement of poor coefficients, reducing noise using filtering approaches, etc. Mathematical details of the steps involved can be referred in [@vishwakarma2017understanding].

A schematic diagram of the code workflow is presented in Fig 01. The module codes can be categorized into four categories: load data, convert data formats, core functionality, and auxiliary codes. The <i>load data</i> codes can read data from either of the `JPL`, `ITSG`, or `CSR` GRACE data centers. The codes further perform the necessary preprocessing, including the conversion of data formats, replacement of some Legendre coefficients as well as removing the long-term mean. The <i>convert data format</i> codes can convert `L2` GRACE Spherical Harmonics data from one format to another. These codes are invoked by the  <i>load data</i> codes for data format conversion. Further, these codes can be independently invoked as well by the user for their needs.

The `core functionalities` of the module are the `gsha`, `gshs`, `gddc`, `tws_cal`, and `basin_avg` codes. `GSHS` module inputs the GRACE L2 spherical harmonic coefficients and performs the `GRACE Spherical Harmonics Synthesis (GSHS)` algorithm. The algorithm converts the input L2 spherical harmonic coefficients into gridded values at the user-desired grid resolution. An inverse module is also provided, called the `gsha.py` module. This module performs the `GRACE Spherical Harmonics Analysis (GSHA)` algorithm. The algorithm converts the gridded `TWSA` values into the GRACE L2 spherical harmonics coefficients. In addition to the translation of the `SHbundle` MATLAB package, this contribution further includes the GRACE Data Driven Correction function, detailed and first coded in MATLAB by [@vishwakarma2017understanding]. The implementation is done via the `gddc` module. More details on the `gddc` implementation can be referred to in the paper cited above. `tws_cal ` computes the total water storage values at each grid from the GRACE L2 spherical harmonics by first applying the `Gaussian filter` on the L2 data products, and then calling the `gshs` module. `basin_avg` module computes the `L4` GRACE basin average time-series product, for any basin input by the user as a GIS shapefile.
<br>

![Schematic diagram of code workflow. \label{fig:code_workflow}](./pic/01_flowchart_20240327.jpg)<br>
<i>Fig 01: Schematic Diagram of the Code Workflow</i>

<br>

The rest of the codes are bundled together as <i>auxiliary codes</i> in Fig 01. An important part of the `GSHS` algorithm implementation is the implementation of the `PLM` algorithm. The `PLM` algorithm inputs degree, order, and co-latitude values and computes the output from normalized Legendre function of the first kind. The `plm.py` module can also provide the first and second derivatives of the Legendre functions. The implementation of the integrals of the Legendre functions is also done; this is available through the `iplm.py` module. `IPLM` inputs the degree, order and co-latitude, and returns the integrated Legendre functions.<br>

Some important modules for the spherical harmonic synthesis step are `normalklm`, `eigengrav`, and `ispec`. `normalklm.py` returns the hydrostatic equilibrium ellipsoid for the earth's surface based on [@lambeck1988geophysical]. `eigengrav.py` provides the isotropic spectral transfer to obtain the equivalent water thickness (m). Lastly, `ispec.py` inputs the sine and cosine coefficients and returns the field function `F`. <br>

The `Global Spherical Harmonic Analysis` code depends upon `neumann` along with the `IPLM` and `sc2cs`. `neumann.py` returns the weights and nodes for Neumann's numerical integration scheme on the sphere. The `gshs.py` code provides options for spherical harmonic synthesis to compute the sine and cosine components of the Legendre function. These include `least squares`, `weighted least squares`, `approximate quadrature`, `first neumann method`, `second neumann method`, and `block mean values`. The `neumann.py` code is required for the implementation of the `first neumann method` and `second neumann method`.<br>

Further, our code has been tested using the `SHbundle` outputs as the benchmark. The validation results are provided as part of the tutorials in the repository.
<br>

# Concluding Remarks

In this paper, we have introduced a new Python software package named `PySHbundle`. The software can process Stokes coefficients for Earth's gravity field to provide gridded products representing changes in mass, geoid height anomalies, equivalent water height anomalies, etc. This package has been specially developed to process the level 2  spherical harmonic data from the GRACE satellite mission, which finds application in numerous disciplines of Earth system science. The software has potential to increase GRACE level-2 data processing accessibility to early-career researchers and researchers from low-income regions, adhering to the principles of open science.

# Acknowledgements

The authors would like to thank Dr.-Ing. Markus Antoni and Clara Buetzler, Institute of Geodesy, University of Stuttgart, Germany, for early feedback to the `PySHbundle` software during its development phase. We are extremely grateful to the financial support from IISc-ISRO Space Technology Cell for funding the project titled "Improving the spatial resolution of GRACE TWS for India using remote sensing datasets and modeling approach" under grant number STC0437.

# References

</p>
