# Load Data

## GRACE Processing Centers

There are 3 major research centers which disseminate GRACE data. These are:

+   [University of Texas at Austin, Center for Space Research (CSR)](https://www.csr.utexas.edu/)
+   [Jet Propulsion Laboratory (JPL)](https://grace.jpl.nasa.gov/)
+   [German Research Center for Geosciences (GFZ)](https://www.gfz-potsdam.de/en/)

## GRACE Data Levels

GRACE data is available at different processing levels:

+ Level 1: Raw satellite data
    * Level 1A: 
    * Level 1B:
+ Level 2: Spherical Harmonic Coefficients for the geopotential estimates
+ Level 3: consists of mass anomalies or other standardized products such as Monthly Ocean/Land Water Equivalent Thickness, Surface-Mass Anomaly. Similarly mass concentration blocks or `mascons` are also available.
+ Level 4: Time-series of catchment level hydrological estimates of TWSA

`PySHBundle` provides the capability to obtain grided Total Water Storage Anomaly(TWSA) from Level 2 data.

## Pre-Processing of Data

::: pyshbundle.read_GRACE_SH_paths
::: pyshbundle.reader_replacer_jpl.reader_replacer_jpl
::: pyshbundle.reader_replacer_csr.reader_replacer_csr
::: pyshbundle.reader_replacer_itsg.reader_replacer_itsg
::: pyshbundle.reader_replacer.reader

## Computing Ling Term Mean
::: pyshbundle.load_longterm_mean

## Physical and Geodetic Constants
::: pyshbundle.GRACEconstants
