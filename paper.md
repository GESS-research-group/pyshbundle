---
title: 'PySHbundle: A Python implementation of MATLAB codes SHbundle'
tags:
  - Python
  - GRACE
authors:
  - name: ABCD
    orcid: 0000-0000-0000-0000
    equal-contrib: true
    affiliation: "1, 2" # (Multiple affiliations must be quoted)
  - name: Author Without ORCID
    equal-contrib: true # (This is how you can denote equal contributions between multiple authors)
    affiliation: 2
  - name: Author with no affiliation
    corresponding: true # (This is how to denote the corresponding author)
    affiliation: 3
affiliations:
 - name: Lyman Spitzer, Jr. Fellow, Princeton University, USA
   index: 1
 - name: Institution Name, Country
   index: 2
 - name: Independent Researcher, Country
   index: 3
date: 13 August 2017
bibliography: paper.bib

# Optional fields if submitting to a AAS journal too, see this blog post:
# https://blog.joss.theoj.org/2018/12/a-new-collaboration-with-aas-publishing
aas-doi: 
aas-journal: 
---

<p style='text-align: justify;'>
# Summary

This is the summary - to be included later.

# Introduction

GRACE stands for the Gravity Recovery and Climate Experiment, a joint satellite mission by NASA, the National Aeronautics and Space Administration and DLR, the German Aerospace Centre. 

| Parameter        |    Details      | 
| -------------    |:--------------:| 
| Start of Mission | 17 March 2002  | 
| End of Mission   | 27 October 2017| 
| Inclination      | 89.0°          | 
| Period           | 94.5 minutes   |  

GRACE consists of two identical satellites orbiting around the earth on the same orbital path. The basic principal of the GRACE satellite operation consists of the monitoring of the intersatellite distance between the twin satellites using microwave ranging system. When the satellite system comes across a mass anomaly, each satellite accelerates or decelerates with a phase lag and the intersatellite distance changes. This change in intersatellite distance is later processed to obtain the magnitude of the mass anamoly. When it comes to the continental land surface, the hydrological processes consist of a major component of the mass anamoly over it. Thus, after making necessary post-processing and substracting different mass anomaly components, such as atmospheric changes, etc, the `total water storage anomaly` (`TWSA`) component can be estimated. The `TWSA` is the sum of the total water components over a vertical extention of the grid area through the earth. Conventionally, it is represented in terms of the `equivalent water height` (`m`). GRACE has a successor, GRACE-FO, which was successfully launched on 22 May 2018.<br>

Three different research centres provide GRACE data. These are the University of Texas Center for Space Research (`CSR`), Jet Propulsion Laboratory (`JPL`), and the German Research Center for Geosciences (GFZ). Further, GRACE data is available at different levels of processing. `Level 1` data refers to the raw satellite data. These are further available as `Level 1A` and `Level 1B`, based on the level of processsing done to the raw data. `Level 2` are the spherical harmonic coefficients for the geospatial potential estimates. These may be accessed through the JPL's Physical Oceanography Distributed Active Archive Center (`PO.DAAC`) or through the Information System and Data Center (`ISDC`). `Level 3` consists of mass anomalies or other standardized products, such as the Monthly Ocean/Land Water Equivalent Thickness Surface-Mass Anomaly. Similarly, mass concentration blocks or `mascons` are also availble. These directly provide the total water storage anomaly over gridded regions, and are available through the three GRACE data centers. The mascon products from the various data centers have some differences, attributed to the difference in post-processing steps and corrections applied by the different data centers. While the mascon results make application of GRACE data easier to a wider audience, use of Level 2 data gives the user the freedom and the flexibility to choose their own post-processing algorithms. The choice of application of mascon data product or Level 2 data product may depend upon the purpose of the exercise and the expertise level of the user on the GRACE data post-processing. In this contribtion, we enable the user to obtain the gridded `TWSA` data over a shapefile from Level 2 data.<br>

`Level 2` GRACE data products may be stored in various data formats. These include `|C\S|`, `/S|C\`, `clm`, vector, and `Colombo` format. Our contribution in its current version can only handle the `|C\S\` and `/S|C\` data formats. In `|C\S|` format, the `Clm` and `Slm` coefficient are stored as lower triangle and upper triangle, respectively, in a matrix of dimention <i>(l + 1) x (l + 1)</i>. In `/S|C\` format, the coefficients are stored in amatrix of dimension <i>(l + 1) x (2 l + 1)</i> with horizontally flipped triangular matrix of `Slm` coefficients on the left half, triangular matrix of `Clm` on the right half, and zeros on the rest of the matrix elements. In our contribution, conversion between the two data formats is made possible with the modules `cs2sc` and `sc2cs`.<br>

# Statement of need
A MATLAB code bundle already exists called `SHbundle` developed by `Sneew et al. (2021)` and distributed under the GNU license. The code bundle can be freely used and modified by anyone giving proper credit to the original developers. However, MATLAB being a proprietary software may have some limitations in terms of accessibility. <i>`Brief description of impact of SHBundle package here`</i><br>

On the other hand, a strong community of programmers also exists for Python, an open-source programming language. In this contribution, we have translated the MATLAB codes from the SHbundle into the Python programming language.<br>
It is hoped the contribution will make GRACE L2 data processing more accessible to a wider audience of programmers. Our python package is titled `PySHbundle` and the working code can be accessed in GitHub : [https://github.com/mn5hk/pyshbundle](https://github.com/mn5hk/pyshbundle)

# Methodology

We have implemented the matlab codes `SHbundle` into the python programming language. More details on the `SHbundle` package may be refered to at `Sneew et al. (2021)`. The naming of the modules and the workflow between the modules has been preserved as much as possible in the `PySHbundle` python implementation. This is to ensure smooth communication between user communitities of the two packages and/or the two different programming language communities. Further, our code has been tested using the `SHbundle` implementation results for validation.

# Features
<i>`[Maybe not necessary; can be covered in the implementation section]`</i><br>
# Implementation
A schematic diagram of the code workflow is presented in the figure below. <br>
![Schematic diagram of code workflow. \label{fig:code_workflow}](https://github.com/mn5hk/pyshbundle/blob/main/pic/flowchart_draft_20221227.png)<br>
<i>Fig 01: Schematic Diagram of the Code Workflow</i><br>

The key module for the package is the `gsha.py` module. This module inputs the GRACE L2 spherical harmonic coefficients and performs the `GRACE Spherical Harmonics Analysis (GSHA)` algorithm. The algorithm converts the input L2 spherical harmonic coefficients into gridded values at the user-desired grid resolution. An inverse module is also provided, called the `gshs.py` module. This module performs the `GRACE Spherical Harmonics Synthesis (GSHS)` algorithm. The algorithm converts the gridded `TWSA` values into the GRACE L2 spherical harmonics coefficients.<br>

An important part of the `GSHS` algorithm implementation is the implementation of the `PLM` algorithm, as shown in the following figure \autoref{fig:code_workflow}. The `PLM` algorithm inputs degree, order and co-latitude values and computes the Legendre functions. The `plm.py` module can also provide the first and second derivatives of the Legendre functions. The implementation of the integrals of the Legendre functions is also done; this is available through the `iplm.py` module. `IPLM` inputs the degree, order and co-latitude, and returns the integrated Legendre functions.<br>

Some important modules for the spherical harmonic synthesis step are `normalklm`, `eigengrav`, and `ispec`. `normalklm.py` returns the hydrostatic equilibrium ellipsoid for the earth surface based on `Lambeck (1988) "Geophysical Geodesy", p.18`. `eigengrav.py` provides the isotrophic spectral transfer to obtain the equivalent water thickness (m). Lastly, `ispec.py` inputs the sine and cosine coefficients and returns the field function `F`. <br>

The `Global Spherical Harmonic Analysis` code depends upon `neumann` along with the `IPLM` and `sc2cs`. `neumann.py` returns the weights and nodes for Neumann's numerical integration scheme on the sphere. The `gshs.py` code provides options for spehrical harmonic synthesis to compute the sine and cosine components of the Legendre function. These include `least squares`, `weighted least squares`, `approximate quadrature`, `first neumann method`, `second neumann method` and `block mean values`. The `neumann.py` code is required for the implementation of the `first neumann method` and `second neumann method`.<br>

In addition to the translation of the `SHbundle` matlab package, this contribution further includes the GRACE Data Driven Correction function, detailed and first coded in matlab in `Vishwakarma et al.  (2017)`. The implementation is done via the `gddc` module. More details on the `gddc` implementation can be refered to in the paper cited above.<br>

# Validation
<i>`[Comparative study over some basins @vivek]`</i>
# Acknowledgements

# References
# End of draft

References: <br>
[^1]. https://www.nasa.gov/mission_pages/Grace/ <br>
2. 
3. 

</p>
# Mathematics

Single dollars ($) are required for inline mathematics e.g. $f(x) = e^{\pi/x}$

Double dollars make self-standing equations:

$$\Theta(x) = \left\{\begin{array}{l}
0\textrm{ if } x < 0\cr
1\textrm{ else}
\end{array}\right.$$

You can also use plain \LaTeX for equations
\begin{equation}\label{eq:fourier}
\hat f(\omega) = \int_{-\infty}^{\infty} f(x) e^{i\omega x} dx
\end{equation}
and refer to \autoref{eq:fourier} from text.

# Citations

- Lambeck, K. (1988). Geophysical geodesy (p. 718). Oxford: Clarendon.
- Nico Sneeuw, Matthias Weigelt, Markus Antoni, Matthias Roth, Balaji Devaraju, et. al. (2021). SHBUNDLE 2021. http://www.gis.uni-stuttgart.de/research/projects/Bundles.
- Vishwakarma, B. D., Horwath, M., Devaraju, B., Groh, A., & Sneeuw, N. (2017). A data‐driven approach for repairing the hydrological catchment signal damage due to filtering of GRACE products. Water Resources Research, 53(11), 9824-9844.

Citations to entries in paper.bib should be in
[rMarkdown](http://rmarkdown.rstudio.com/authoring_bibliographies_and_citations.html)
format.

If you want to cite a software repository URL (e.g. something on GitHub without a preferred
citation) then you can do it with the example BibTeX entry below for @fidgit.

For a quick reference, the following citation commands can be used:
- `@author:2001`  ->  "Author et al. (2001)"
- `[@author:2001]` -> "(Author et al., 2001)"
- `[@author1:2001; @author2:2001]` -> "(Author1 et al., 2001; Author2 et al., 2002)"

# Figures

Figures can be included like this:
![Caption for example figure.\label{fig:example}](figure.png)
and referenced from text using \autoref{fig:example}.

Figure sizes can be customized by adding an optional second parameter:
![Caption for example figure.](figure.png){ width=20% }





