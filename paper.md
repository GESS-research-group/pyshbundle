---
title: 'PySHbundle: A Python implementation of MATLAB codes SHbundle'
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
    affiliation: 1
  - name: Vivek Yadav
    equal-contrib: true # (This is how you can denote equal contributions between multiple authors)
    affiliation: 1
  - name: Tsungrojungla Walling
    affiliation: 2
  - name: Maya Suryawanshi
    affiliation: 1
  - name: Bramha Dutt Vishwakarma
    orcid: 0000-0003-4787-8470
    corresponding: true # (This is how to denote the corresponding author)
    affiliation: "1,3" # (Multiple affiliations must be quoted)
affiliations:
 - name: Interdisciplinary Centre for Water Research, Indian Institute of Science, India
   index: 1
 - name: Undergraduate Programme, Indian Institute of Science, India
   index: 2
 - name: Centre of Earth Science, Indian Institute of Science, India
   index: 3
date: 13 February 2023
bibliography: paper.bib

# Optional fields if submitting to a AAS journal too, see this blog post:
# https://blog.joss.theoj.org/2018/12/a-new-collaboration-with-aas-publishing
# aas-doi: 
# aas-journal: 
---

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

Three different research centres provide GRACE data. These are the University of Texas Center for Space Research (`CSR`), Jet Propulsion Laboratory (`JPL`), and the German Research Center for Geosciences (`GFZ`). Further, GRACE data is available at different levels of processing. `Level 1` data refers to the raw satellite data. These are further available as `Level 1A` and `Level 1B`, based on the level of processsing done to the raw data. `Level 2` are the spherical harmonic coefficients for the geospatial potential estimates. These may be accessed through the JPL's Physical Oceanography Distributed Active Archive Center (`PO.DAAC`) or through the Information System and Data Center (`ISDC`). `Level 3` consists of mass anomalies or other standardized products, such as the Monthly Ocean/Land Water Equivalent Thickness Surface-Mass Anomaly. Similarly, mass concentration blocks or `mascons` are also availble. These directly provide the `TWSA` over gridded regions, and are available through the three GRACE data centers. More details on the mascon approaches for studying gravity fields and the approaches used by the different data centers for generating mascon products may be referred to in `Antoni (2022)`. The mascon products from the various data centers have some differences, attributed to the difference in post-processing steps and corrections applied by the different data centers. An online tool exists developed by the `Colorado Center for Astrodynamics Research` (https://ccar.colorado.edu/grace). This tool can be used for a quick visualization of the `CSR` and `GSFC` mascon products over all regions of the globe. While the mascon results make application of GRACE data easier to a wider audience, use of `Level 2` data gives the user the freedom and the flexibility to choose their own post-processing algorithms. The choice of application of mascon data product or Level 2 data product may depend upon the purpose of the exercise and the expertise level of the user on the GRACE data post-processing. In this contribtion, we enable the user to obtain the gridded `TWSA` data over a shapefile from Level 2 data.<br>

`Level 2` GRACE data products may be stored in various data formats. These include `|C\S|`, `/S|C\`, `clm`, vector, and `Colombo` format. Our contribution in its current version can only handle the `|C\S\` and `/S|C\` data formats. In `|C\S|` format, the `Clm` and `Slm` coefficient are stored as lower triangle and upper triangle, respectively, in a matrix of dimention <i>(l + 1) x (l + 1)</i>. In `/S|C\` format, the coefficients are stored in amatrix of dimension <i>(l + 1) x (2 l + 1)</i> with horizontally flipped triangular matrix of `Slm` coefficients on the left half, triangular matrix of `Clm` on the right half, and zeros on the rest of the matrix elements. In our contribution, conversion between the two data formats is made possible with the modules `cs2sc` and `sc2cs`.<br>

`Level 3` products are the catchment average hydrological estimates of `TWSA`. These are obtained through the further processing of `Level 2` products. `Level 3` products may further be processed to obtain catchment average timeseries data, labelled as `Level 4` products. Various tools exist in the literature to process GRACE data and to analyze it. Some of these available in the `Matlab` programming language are: `SHbundle` (`Sneew et al., 2021`), `GRAMAT` (`Feng, 2019`), `SHADE` (`Piretzidis, D., & Sideris, M. G., 2018`), etc. Similarly, some GRACE data processing tools are also available based on the python programming language. These include `gravity-toolkit` (https://github.com/tsutterley/gravity-toolkit), `ggtools` (https://pypi.org/project/ggtools/) and `GRACE-filter` (https://github.com/strawpants/GRACE-filter). Genral tools for spheric harmonic analysis are also available, such as SHTools (`Wieczorek, M. A., & Meschede, M., 2018`). `SHBundle` provides MATLAB-tools for `spheric harmonic synthesis` and `spheric harmonic analysis`. The earliest version of the code were developed in 1994 while the latest version with upgrades can be found dated 2018. `GRAMAT` provides a similar MATLAB-based tools for processing GRACE spherical harmonics data to obtain spatiotemporal global mass variations. The GRAMAT toolbox includes Gaussian smoothening filter to remove North-South stripes, spherical harmonic analysis and synthesis routines, leakage effect reduction routines, harmonic analysis of times series over regions, and uncertainty analysis of GRACE estimates (`Feng, 2019`). `SHADE` provides a matlab-based toolbox for the empirical de-correlation of GRACE monthly spherical harmonics (`Piretzidis, D., & Sideris, M. G., 2018`). `gravity-toolkit` is a python-based package meant to handle GRACE L2 data products. Its functionalities include visualization of GRACE and GRACE-FO L2 data products, and the estimation of GRACE and GRACE-FO L2 data product errors. `gg-tools` too contain similar tools for signal correction and for conversion of GRACE L2 products to L3. `GRACE-filter` provides tool for filtering of GRACE L2 product using DDK filter based on `Kusche et al. (2009)`.
 
# Statement of need
A MATLAB code bundle already exists called `SHbundle` developed by `Sneew et al. (2021)` and distributed under the GNU license. The code bundle can be freely used and modified by anyone giving proper credit to the original developers. However, MATLAB being a proprietary software may have some limitations in terms of accessibility. <i>`Brief description of impact of SHBundle package here`</i><br>

On the other hand, a strong community of programmers also exists for Python, an open-source programming language. In this contribution, we have translated the MATLAB codes from the SHbundle into the Python programming language. In addition to the SHBundle codes, we have further translated the `GRACE Data Driven Correction (GDDC)` codes from Matlab to Python. `GDDC` allows the correction of filtered GRACE `Level 2` products and restore the signal loss, independent of the catchment size `(Vishwakarma et al., 2017)`.<br>
It is hoped the contribution will make GRACE L2 data processing more accessible to a wider audience of programmers. Our python package is titled `PySHbundle` and the working code can be accessed in GitHub : [https://github.com/mn5hk/pyshbundle](https://github.com/mn5hk/pyshbundle)

# Methodology

We have implemented the matlab codes `SHbundle` into the python programming language. More details on the `SHbundle` package may be refered to at `Sneew et al. (2021)`. The naming of the modules and the workflow between the modules has been preserved as much as possible in the `PySHbundle` python implementation. This is to ensure smooth communication between user communitities of the two packages and/or the two different programming language communities. Further, our code has been tested using the `SHbundle` implementation results for validation.

# Mathematics

In this section, we present a mathematical representation of the spherical harmonics analysis. According to potential theory, the gravitational field of a body fulfils the Laplace equation $\nabla^2\phi = 0$. Laplace's equation in spherical coordinates can be written as follows: 

$$
\begin{equation}
\frac{1}{r^2}\frac{\partial}{\partial r}\bigg( r^2\frac{\partial \phi}{\partial r}\bigg)  
+
\frac{1}{r^2\sin\vartheta}\frac{\partial}{\partial \vartheta}\bigg(\sin\vartheta\frac{\partial \phi}{\partial \vartheta}\bigg) 
+
\frac{1}{r^2\sin^2\vartheta}\frac{\partial^2 \phi}{\partial \lambda^2}
 = 0 ,
\end{equation}\\
$$

where 
$\phi$ is the potential, 
$r$ is the radius, 
$\vartheta$ is the co-latitude and 
$\lambda$ is the longitude. 

We perform a separation of variables and insert $\phi(r,\vartheta,\lambda) =f(r)g(\vartheta)h(\lambda)$ into the Laplace equation to get three independent equations:

$$
\begin{equation}
r^2\frac{d^2f}{dr^2}+2r\frac{df}{dr} - n(n+1)f = 0,
\end{equation}
$$
$$
\begin{equation}
\frac{d^2g}{d\vartheta^2}
+
\frac{dg}{d\vartheta}\cot\vartheta
+
\bigg(  n(n+1) - \frac{m^2}{\sin^2\vartheta}   \bigg) g = 0 ,
\end{equation}
$$
$$
\begin{equation}
\frac{d^2h}{d\lambda^2} + m^2h = 0,
\end{equation}
$$

where $m$ and $n$ are the degree and order respectively. Solving $(2), (3)$ and $(4)$, we obtain: 

$$
\begin{equation}
f(r) \in \{r^n, r^{-(n+1)}\},
\end{equation}
$$
$$
\begin{equation}
g(\vartheta) \in \{P_{n,m}(\cos \vartheta), Q_{n,m}(\cos \vartheta)\} ,
\end{equation}
$$
$$
\begin{equation}
h(\lambda) \in \{\cos m\lambda, \sin m\lambda\}.
\end{equation}\\
$$

Thus, the Laplace equation's solution takes the following form: 

$$
\begin{equation}
\phi(r, \vartheta, \lambda) = \sum_{n=0}^{\infty} \sum_{m=0}^{n} 
\alpha_{n,m}
\begin{Bmatrix}
P_{n,m}(\cos\vartheta)\\
Q_{n,m}(\cos\vartheta)\\
\end{Bmatrix}
\dot{•}
\begin{Bmatrix}
\cos m\lambda\\
\sin m\lambda\\
\end{Bmatrix}
\dot{•}
\begin{Bmatrix}
r^n\\
r^{(n+1)}\\
\end{Bmatrix}
.
\end{equation}
$$

Solutions for $f(r)$ and $h(\lambda)$ are fairly straightforward. Eq - (3) for $g(\vartheta)$ is in the form of a Legendre differential equation and its solutions are $P_{n,m}(\cos \vartheta)$ and $Q_{n,m}(\cos \vartheta)$, the associated Legendre functions of the first and second kind. We now apply two constraints to the solution:

* $\phi \rightarrow 0$ when $r \rightarrow \infty$,
* $\phi$ is limited on the sphere,

which leads us to eliminate $Q_{n,m}(\cos \vartheta)$ and $r^n$.The $4\pi$ normalization of the Associated Legendre functions [8] is utilized in our package and is given by: 

$$
\begin{equation}
\bar{P}_{n,m}(\cos\vartheta) = P_{n,m}(\cos\vartheta)\sqrt{(2-\delta_{m0})(2n+1)\frac{(n-m)!}{(n+m)!}},
\end{equation}
$$
where $\delta_{m0}$ is the Kronecker delta function,
$$
\begin{equation}
P_{n,m}(t) = (1-t^2)^{\frac{m}{2}}\frac{d^mP_n(t)}{dt^m},
\end{equation}
$$
and 
$$
\begin{equation}
nP_n(t)=-(n-1)P_{n-2}(t) + (2n-1)tP_{n-1}(t).
\end{equation}
$$

Spherical harmonics are the angular portion of a set of solutions to Laplace's equation. They take into account $\vartheta$ and $\lambda$. They are functions modelled on the surface of a sphere, denoted by $Y_{n,m}(\vartheta,\lambda)$. They are of three kinds: 

* Zonal harmonics: $m=0$ - they are only latitude dependent,
* Tesseral harmonics: $0 < m < n$, and 
* Sectorial harmonics: $m=n$.

Quantities like the gravitational potential, height of water column, gravity anomaly and so on are the functionals of the gravity field which are obtained by differentiating the potential $\phi$ with respect to the spherical coordinates. 

The gravitational potential anomaly $V$ is given by:
$$
\begin{equation}
    V(r, \vartheta, \lambda) = 
    \frac{GM}{r} \sum_{n=0}^{N_{max}} \sum_{m=0}^{n} 
    \left(\frac{R}{r}\right)^{n+1}
    \bar{P}_{n,m}(\cos \vartheta)[C_{n,m}\cos m\lambda+S_{n,m}\sin m\lambda].
\end{equation}
$$

Here, $R$ refers to the radius of the Earth,
$\bar{P}_{n,m}$ refers to the Associated Legendre functions with $4\pi$ normalization,
$C_{lm}$ and  $S_{lm}$ refer to the spherical harmonic coefficients. Similarly, another functional, the change in surface mass density, is represented by:

$$
\begin{equation}
    \Delta\sigma(\vartheta,\lambda) = 
    \frac{a\rho_{ave}}{3} 
    \sum_{n=0}^{N_{max}} \sum_{m=0}^{n} 
    \left(\frac{R}{r}\right)^{n+1} 
    \bar{P}_{n,m}(\cos \vartheta)
    \frac{2n+1}{1+k_l}
    [C_{n,m}\cos m\lambda+S_{n,m}\sin m\lambda],
\end{equation}
$$
where
$\rho_{ave}$ refers to the average density of the Earth in $g/cm^3$ and
$k_n$ refers to the load Love number of degree $n$.

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

1. https://www.nasa.gov/mission_pages/Grace/ <br>

# Citations

- Antoni, M. (2022). A review of different mascon approaches for regional gravity field modelling since 1968. History of Geo-and Space Sciences, 13(2), 205-217.
- Feng, W. GRAMAT: a comprehensive Matlab toolbox for estimating global mass variations from GRACE satellite data. Earth Sci Inform 12, 389–404 (2019). https://doi.org/10.1007/s12145-018-0368-0
- Kusche, J., Schmidt, R., Petrovic, S., & Rietbroek, R. (2009). Decorrelated GRACE time-variable gravity solutions by GFZ, and their validation using a hydrological model. Journal of geodesy, 83, 903-913.
- Lambeck, K. (1988). Geophysical geodesy (p. 718). Oxford: Clarendon.
- Nico Sneeuw, Matthias Weigelt, Markus Antoni, Matthias Roth, Balaji Devaraju, et. al. (2021). SHBUNDLE 2021. http://www.gis.uni-stuttgart.de/research/projects/Bundles.
- Piretzidis, D., & Sideris, M. G. (2018). SHADE: A MATLAB toolbox and graphical user interface for the empirical de-correlation of GRACE monthly solutions. Computers & Geosciences, 119, 137-150.
- Vishwakarma, B. D., Horwath, M., Devaraju, B., Groh, A., & Sneeuw, N. (2017). A data‐driven approach for repairing the hydrological catchment signal damage due to filtering of GRACE products. Water Resources Research, 53(11), 9824-9844.
- Wieczorek, M. A., & Meschede, M. (2018). SHTools: Tools for working with spherical harmonics. Geochemistry, Geophysics, Geosystems, 19(8), 2574-2592.
</p>




