---
title: 'PySHbundle: A Python tool for processing GRACE Gravimetry data into Global Surface Mass Change Datasets'
tags:
  - Python
  - GRACE
  - Geodesy
  - Gravimetry
  - Terrestrial Water Storage
  - Spherical Harmonic Analysis
  - Spherical Harmonic Synthesis
  - GRACE Data Driven Correction
authors:
  - name: Vivek Kumar Yadav
    orcid: 0009-0000-7156-4450
    equal-contrib: true # (This is how you can denote equal contributions between multiple authors)
    corresponding: true
    affiliation: 1
  - name: Amin Shakya
    orcid: 0000-0002-4706-826X
    equal-contrib: true # (This is how you can denote equal contributions between multiple authors)
    affiliation: "1,2"
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
 - name: Interdisciplinary Centre for Water Research, Indian Institute of Science, India
   index: 1
 - name: Faculty of Geo-Information Science and Earth Observation, University of Twente, the Netherlands
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

`GRACE` (Gravity Recovery and Climate Experiment) satellite mission has been mapping mass changes near the surface of the Earth since 2002. One of the major mechanisms of short term mass transport is the redistribution of water, GRACE has significantly influenced Geosciences. GRACE satellite products are typically released at various levels of complexity, often referred to as processing levels. Level 1 is the satellite instrument data that is processed to obtain Level 2 (`L2`) the Spherical harmonic coefficients and standard deviations of the static gravity field, for a particular time period. `L2` are further processed to obtain Level 3 products; global gridded mass change estimates (`L3`) expressed as terrestrial water storage anomalies (`TWSA`). The `L2` spherical harmonic data are typically noisy, which necessitates the use of spectral filtering. The data also have to be corrected for known artifacts and contaminating geophysical signals, such as solid Earth processes in the case of isolating TWSA. Processing choices, such as filter properties and type, have a significant impact on the accuracy and the resolution of final gridded output. Therefore, most `L3` users must be cautious when using GRACE data for specific applications. The majority of the GRACE data user community is not well versed with `L2` data processing, and most often use the off-the-shelf `L3` products. Here we developed an open-source processing toolbox to provide users with more control over processing choices. A python module, called PySHbundle, was developed to ease the conversion of GRACE `L2` Spherical Harmonics data products to `L3` `TWSA` products. With this contribution, we hope to enable further usage of GRACE data for Earth system science.


# Introduction

The NASA/DLR GRACE and NASA/GFZ GRACE-FO twin satellite missions measure changes in the Earth's gravitational field by measuring their inter-satellite distance. Changes in the local gravity field affect the orbit of each satellite, which is recorded with the onboard ranging system [@wahr1998time]. When the satellite pair comes in the vicinity of a temporal mass anomaly, the relative inter-satellite distance changes and it can be inverted to estimate the mass change near the surface of the Earth. Over the continental land surface, the hydrological processes are the major driver of the variation in mass anomaly at monthly to decadal scales. However various other signals such as oceanic and atmospheric variations, high frequency tidal mass changes, systemic correlated errors, etc. are also part of the obtained GRACE signals [@humphrey2023using]. 

Several researchers in Geosciences use level three GRACE data, which is obtained from `L2` Spherical harmonic coefficients, except JPL MASCONS which are derived from Level-1B satellite ranges [@watkins2015mascons]. The procedure to convert `L2` to `L3` is called spherical harmonic synthesis. However, there are several pre-processing steps; such as anomaly calculation, replacing poor quality low degree coefficients, filtering, and correcting for signal damage due to filtering.

A few GRACE data processing tools are available based on the Python programming language. These include [`gravity-toolkit`](https://gravity-toolkit.readthedocs.io/en/latest/) [@gravity-toolkit], [`ggtools`](https://pypi.org/project/ggtools/1.1.0/) [@ggtools] and [`shxarray`](https://github.com/ITC-Water-Resources/shxarray) [@shxarray]. General tools for spherical harmonic analysis are also available, such as [`SHTools`](https://agupubs.onlinelibrary.wiley.com/doi/full/10.1029/2018GC007529) [@wieczorek2018shtools]. [`SHbundle`](https://www.gis.uni-stuttgart.de/en/research/downloads/shbundle) provide MATLAB scripts for Spherical Harmonic Synthesis and Spherical Harmonic Analysis. The first version of the code was developed in 1994 while the latest version was released in 2021.
 
# Statement of need

Processing choices introduce subtle differences in the final product, potentially affecting results. Processing `L2` data offers flexibility for users to explore GRACE data for specific applications. This software aims to simplify access to `L2` products, allowing users to select different processing options.

The software processes widely used `L2` products from CSR, JPL, and GFZ. It closely follows the structure of the Matlab-based [`SHbundle`](https://www.gis.uni-stuttgart.de/en/research/downloads/shbundle) and [`GRACE Data Driven Correction (GDDC)`](https://doi.org/10.1002/2017WR021150)[@vishwakarma2017data] codes, enabling cross-compatibility between Python and Matlab users.

`PySHbundle` is modular, offering tools to process GRACE data, including anomaly computation, low-degree coefficient substitution, noise reduction, handling gaps and signal leakage correction. It supports future development for hydrological applications. While excellent tools exist (above) for Spherical Harmonic operation, PySHbundle provides familiarity of existing and beginner level users by focusing on GRACE applications, and translating the legacy software, SHbundle.

By using Python and the GNU license, the package is accessible globally and aligns with the [FAIR principles](https://www.go-fair.org/fair-principles/). We aim to reduce technical and financial barriers, making it useful for researchers, students, and educational programs like the [GRACE Hackweek](https://www.quantumfrontiers.de/de/aktuelles/veranstaltungen/details/news/grace-hackweek-3) at IIT Kanpur.

# Implementation

 Mathematical details of the steps involved can be referred in [@vishwakarma2017understanding].
Accordingly, the package consists of four main modules, `io`, `vizutils`, `pysh_core` and `shutils`.

1. `io`: extract the `L2` coefficients from any of `JPL`, `CSR` and `GFZ` solutions. Followed by replacing the poorly measured degree 1, 2 and 3 spherical harmonics coefficients with recommended datasets. Note that degree 1 coefficients, which represents the center-of-mass of the Earth, are inherently zero since mass of Earth is constant for practical purposes. 

2. `vizutils`: plots the `L2` data to visually understand the coefficients, their uncertainties, mathematical functions used for further processing. 

3. `pysh_core`: Scripts for the global spherical harmonics synthesis `gshs` to convert the `L2` data to global gridded `TWSA` data (`L3`). Calculating signal leakage (`gddc`), and basin-scale average (`Basinaverage`).

4. `shutils`:  Helper scripts for applying `pysh_core`.
Based on the main modules, we provide examples as jupyter notebooks for understanding and using spherical harmonics data and the package.

The accuracy of the core capabilities of `pyshbundle` was validated against the MATLAB `SHbundle` software using 60 months of JPL GRACE RL06 Level-2 data. The gridwise Root Mean Square Error (RMSE) between the computed and reference Terrestrial Water Storage fields was below $10^{-3}$ mm equivalent water height, confirming numerical equivalence between the two implementations. A reproducible validation notebook is available at `examples/validation_pyshbundle`, and the automated test suite in `tests/` enforces this criterion on every commit via continuous integration.


# Acknowledgements

The authors would like to thank Dr.-Ing. Markus Antoni and Clara Buetzler, Institute of Geodesy, University of Stuttgart, Germany, for early feedback. We are grateful for the financial support from IISc-ISRO Space Technology Cell for funding the project titled "Improving the spatial resolution of GRACE TWS for India using remote sensing datasets and modeling approach" under grant number STC0437. BDV would like to acknowledge the financial support from Science and Engineering Research Board, Government of India, under the grant agreement number SRG/2022/000625 for the MATRA project.

# Future Plan

The package will be under continuous development to process data from more research centres, add more filtering and processing algorithms.

# References