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

# Summary

This is the summary

# Introduction

GRACE stands for the Gravity Recovery and Climate Experiment, a joint satellite mission by NASA, the National Aeronautics and Space Administration and DLR, the German Aerospace Center. Basic stats: (do I need to add reference links?)

| Parameter        |    Answer      | 
| -------------    |:--------------:| 
| Start of Mission | 17 March 2002  | 
| End of Mission   | 27 October 2017| 
| Inclination      | 89.0°          | 
| Period           | 94.5 minutes   |  

The operating principle is as follows: 
GRACE has a successor, GRACE-FO, which was successfully launched on 22 May 2018. 

References: <br>
1. 
2. 
3. 

# Statement of need

This is the statement of need. Python enables wrapping low-level languages (e.g., C) for speed without losing
flexibility or ease-of-use in the user-interface. The API for `Gala` was
designed to provide a class-based and user-friendly interface to fast (C or
Cython-optimized) implementations of common operations such as gravitational
potential and force evaluation, orbit integration, dynamical transformations,
and chaos indicators for nonlinear dynamics. `Gala` also relies heavily on and
interfaces well with the implementations of physical units and astronomical
coordinate systems in the `Astropy` package [@astropy] (`astropy.units` and
`astropy.coordinates`).

`Gala` was designed to be used by both astronomical researchers and by
students in courses on gravitational dynamics or astronomy. It has already been
used in a number of scientific publications [@Pearson:2017] and has also been
used in graduate courses on Galactic dynamics to, e.g., provide interactive
visualizations of textbook material [@Binney:2008]. The combination of speed,
design, and support for Astropy functionality in `Gala` will enable exciting
scientific explorations of forthcoming data releases from the *Gaia* mission
[@gaia] by students and experts alike.

# Methodology

We have implemented the matlab codes `SHbundle` into the python programming language. More details on the `SHbundle` package may be refered to at `@Sneew:2018`. The naming of the modules and the workflow between the modules has been preserved as much as possible in the `PySHbundle` python implementation. This is to ensure smooth communication between user communitities of the two packages and/or the two different programming language communities. Further, our results have been tested using the `SHbundle` implementation results as the validation data.

# Implementation
A schematic diagram of the code workflow is presented in the figure below. <br>

The key module for the package is the `gsha.py` module. This module inputs the GRACE L2 spherical harmonic coefficients and performs the `GRACE Spherical Harmonics Analysis (GSHA)` algorithm. The algorithm converts the input L2 spherical harmonic coefficients into gridded values at the user-desired grid resolution. An inverse module is also provided, called the `gshs.py` module. This module performs the `GRACE Spherical Harmonics Synthesis (GSHS)` algorithm. The algorithm converts the gridded total water storage anomaly (TWSA) values into the GRACE L2 spherical harmonics coefficients.<br>

An important part of the `GSHS` algorithm implementation is the implementation of the `PLM` algorithm, as shown in the following figure. The `PLM` algorithm inputs degree, order and co-latitude values and computes the Legendre functions. The `plm.py` module can also provide the first and second derivatives of the Legendre functions. The implementation of the integrals of the Legendre functions is also done; this is available through the `iplm.py` module. `IPLM` inputs the degree, order and colatitude, and returns the integrated Legendre functions.<br>

Briefly also describe `eigengrav.py` here.<br>

Another important module available in the package is the `gddc` module. This module implements the GRACE data driven corrections.

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




# Features

# Acknowledgements



# References
