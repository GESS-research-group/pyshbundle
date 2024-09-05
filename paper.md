---
title: 'PySHbundle: A Python software to convert GRACE Spherical Harmonic Coefficients to gridded mass change fields'
tags:
  - Python
  - GRACE
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

`GRACE` (Gravity Recovery and Climate Experiment) satellite mission has been mapping mass changes near the surface of the Earth since 2002. Since mass redistribution at shorter temporal scales is dominated by hydrology, GRACE has transformed our understanding of changes in the hydrosphere. GRACE data has been used for monitoring and studying groundwater depletion, floods, droughts, etc. GRACE satellite products are typically released at various levels of complexity, often referred to as processing levels. Level 1 is the satellite instrument data that is processed to obtain GRACE Spherical Harmonics data level 2( `L2`). `L2` data: represents the mean monthly gravity field of the Earth. `L2` are further processed to obtain level 3 products; global gridded mass change estimates (`L3`) expressed as terrestrial water storage anomalies (`TWSA`). `L2` data are unconstrained gravity field solutions and are noisy, which are filtered and corrected for known artifacts and signals from solid Earth processes to obtain `L3` products that are useful for hydrology. Processing choices, such as filter properties and type, have a significant impact on the accuracy and the resolution of final gridded output. Therefore, `L3` users must be cautious when using GRACE data for specific applications. Since, majority of the GRACE data user community is not well versed with `L2` data processing, they often use off the shelf product with doubts on the efficacy of GRACE mission. Here we developed an open-source processing toolbox to provide users with more control over processing choices. A python module, called PySHbundle, is developed that converts GRACE `L2` Spherical Harmonics data products to `L3` `TWSA` products while applying the data-driven correction algorithm for reducing the impact of filtering on signal. With this contribution, we hope to enable further usage of GRACE data for Earth system science.

# Introduction

The mission measures changes in the inter-satellite distance with a microwave ranging system that gives an accuracy in the range of micrometers [@wahr1998time]. When the satellite system comes in the vicinity of a temporal mass anomaly, the relative inter-satellite distance changes and it can be inverted to estimate the mass change near the surface of the Earth. Over the continental land surface, the hydrological processes are the major driver of the variation in mass anomaly at monthly to decadal scales. However various other signals such as oceanic and atmospheric variations, high frequency tidal mass changes, systemic correlated errors, etc. are also part of the obtained GRACE signals [@humphrey2023using]. 

Obtaining `L3` products from `L1` requires 
	a) removing Oceanic, Atmospheric, Tidal and other signals
		Isolating the hydrological signal results in host of different atmosphere, ocean and tidal models to be applied.
	b) filtering to reduce noise
		Several methods have been proposed for filtering noise in `L2` data, ranging from simple Gaussian averaging to dynamic filters which use hydrological models.
	c) processing data to global gridded TWSA
  
		Finally, the method applied to convert `L2` data to `L3` `TWSA` data introduces more subtle differences.
Various tools exist to process GRACE data and to analyze it. Some of these are developed in the `MATLAB` programming language: [`SHbundle`](https://www.gis.uni-stuttgart.de/en/research/downloads/shbundle) [@SHbundle], [`GRACE Data Driven Correction`](https://www.gis.uni-stuttgart.de/en/research/downloads/datadrivencorrectionbundle) [@vishwakarma2017understanding],  [`LUH-GRACE2018`](https://www.ife.uni-hannover.de/en/services/luh-grace) [@koch2020luh], [`GRAMAT`](https://link.springer.com/article/10.1007/s12145-018-0368-0) [@feng2019gramat], [`SHADE`](https://www.sciencedirect.com/science/article/pii/S0098300418302760) [@piretzidis2018shade], [`GRACETOOLS`](https://www.mdpi.com/2076-3263/8/9/350) [@darbeheshti2018gracetools], [`SSAS GRACE filter`](https://github.com/shuang-yi/SSAS-GRACE-filter)[@yi2022novel], etc. Similarly, some GRACE data processing tools are also available based on the python programming language. These include [`gravity-toolkit`](https://gravity-toolkit.readthedocs.io/en/latest/) [@gravity-toolkit], [`ggtools`](https://pypi.org/project/ggtools/1.1.0/) [@ggtools] and [`GRACE-filter`](https://github.com/strawpants/GRACE-filter) [@GRACEfilter]. General tools for spheric harmonic analysis are also available, such as [`SHTools`](https://agupubs.onlinelibrary.wiley.com/doi/full/10.1029/2018GC007529) [@wieczorek2018shtools]. [`SHbundle`](https://www.gis.uni-stuttgart.de/en/research/downloads/shbundle) provide MATLAB scripts for `spheric harmonic synthesis` and `spheric harmonic analysis`. The first version of the code was developed in 1994 while the latest version with upgrades can be found dated 2018.
 
# Statement of need

Processing choices (a, b & c) introduce subtle differences in the final `L3` product, potentially affecting results. Processing `L2` data offers flexibility for users to explore GRACE data for specific applications. This software aims to simplify access to `L2` products, allowing users to select different processing options.

The software processes widely used `L2` products from CSR, JPL, and GFZ. It closely follows the structure of the Matlab-based [`SHbundle`](https://www.gis.uni-stuttgart.de/en/research/downloads/shbundle) and [`GRACE Data Driven Correction (GDDC)`](https://doi.org/10.1002/2017WR021150)[@vishwakarma2017data] codes, enabling cross-compatibility between Python and Matlab users.

`PySHbundle` is modular, offering tools to process GRACE data, including anomaly computation, low-degree coefficient substitution, noise reduction, and signal leakage correction. It supports future development for hydrological applications.

By using Python and the GNU license, the package is accessible globally and aligns with the [FAIR principles](https://www.go-fair.org/fair-principles/). We aim to reduce technical and financial barriers, making it useful for researchers, students, and educational programs like the [GRACE Hackweek](https://www.quantumfrontiers.de/de/aktuelles/veranstaltungen/details/news/grace-hackweek-3) at IIT Kanpur.

# Implementation

Obtaining gridded fields from GRACE spherical harmonic coefficients consists of several steps, including obtaining the spherical harmonic coefficients from the data providers, replacement of poor coefficients, reducing noise using filtering approaches, etc. Mathematical details of the steps involved can be referred in [@vishwakarma2017understanding].
Accordingly, the package consists of four main modules, `io`, ``, `` and `shutils`.
1. `io`: extract the `L2` coefficients from any of `JPl`, `CSR` and `ITGZ` solutions. Followed by replacing the poorly measured degree 1, 2 and 3 spherical harmonics coefficients with data from Satellite Laser Ranging missions. And other pre-processing operations.
2. `vizutils`: plots the `L2` data to visually understand the coefficients, their uncertainties, mathematical functions used for further processing. 
3. `pysh_core`: Scripts for the global spherical harmonics synthesis `gshs` to convert the `L2` data to global gridded `TWSA` data (`L3`). Others, calculating signal leakage (`gddc`) and basin-scale average (`Basinaverage`).
4. `shutils`:  Helper scripts for applying `pysh_core`.
Based on the main modules, we provide examples as jupyter notebooks for understanding and using spherical harmonics data and the package.

# Concluding Remarks

In this paper, we have introduced a new Python software package named `PySHbundle`. The software can process Stokes coefficients for Earth's gravity field to provide gridded products representing changes in mass, geoid height anomalies, equivalent water height anomalies, etc. This package has been specially developed to process the`L2`  spherical harmonic data from the GRACE satellite mission, which finds application in numerous disciplines of Earth system science. The software has potential to increase GRACE level-2 data processing accessibility to early-career researchers and researchers from low-income regions, adhering to the principles of open science.

# Acknowledgements

The authors would like to thank Dr.-Ing. Markus Antoni and Clara Buetzler, Institute of Geodesy, University of Stuttgart, Germany, for early feedback. We are grateful to the financial support from IISc-ISRO Space Technology Cell for funding the project titled "Improving the spatial resolution of GRACE TWS for India using remote sensing datasets and modeling approach" under grant number STC0437.

# Future Plan

The package will be under continuous development to process data from more research centres, add more filtering and processing algorithms.
# References

</p>