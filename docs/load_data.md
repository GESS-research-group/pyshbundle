# Load Data

## GRACE Processing Centers

There are 3 major research centers which disseminate GRACE data. These are:

+   [University of Texas at Austin, Center for Space Research (CSR)](https://www.csr.utexas.edu/)
+   [Jet Propulsion Laboratory (JPL)](https://grace.jpl.nasa.gov/)
+   [GeoForschungsZentrum Potsdam (GFZ)](https://www.gfz-potsdam.de/en/)

## GRACE Data Levels

The `GRACE` data products are being developed, processed  and archieved in a shared Science Data System between the `Jet Propulsion Laboratory(JPL)`, the `University of Texas Center for Space Research (UT-CSR)` and `GeoForschungsZentrum Potsdam (GFZ)`.

+ Level 0:<br>
    The level-0 data are the result of the data reception, collection and decommutation by the Raw Data Center (RDC) of the Mission Operation System (MOS) located in Neustrelitz, Germany. The MOS receives twice per day using its Weilheim and Neustrelitz tracking antennae the science instrument and housekeeping data from each GRACE satellite which will be stored in two appropriate files in the level-0 rolling archive at DFD/Neustrelitz. The SDS retrieves these files and extracts and reformats the orresponding instrument and ancillary housekeeping data like GPS navigation solutions,space segment temperatures or thruster firing events. Level-0 products are available 24-hours after data reception.

+ Level 1:<br>
    The level-1 data are the preprocessed, time-tagged and normal-pointed instrument data. These are the K-band ranging, accelerometer, star camera and GPS data of both satellites. Additionally the preliminary orbits of both GRACE satellites will be generated. Level-1 data processing software is developed by JPL with support from GFZ (e.g. accelerometer data preprocessing). Processing of level-1 products is done primarily at JPL. An identical processing system (hardware/software) is installed at GFZ to serve as a backup system in case of hardware or network problems. This double implementation is necessary to guarantee the envisaged level-1 product delay of 5 days. All level-1 products are archived at JPL’s Physical Oceanography Distributed Active Data Center(PODAAC) and at GFZ’s Integrated System Data Center (ISDC) . Both archives are harmonized on a sub-daily timeframe.

+ Level 2: `Spherical Harmonic Coefficients` for the geopotential

    Level-2 data include the short term (30 days) and mean gravity field derived from calibrated and validated GRACE level-1 data products. This level also includes ancillary data sets (temperature and pressure fields, ocean bottom pressure, and hydrological data) which are necessary to eliminate time variabilities in gravity field solutions. Additionally the precise orbits of both GRACE satellites are generated. All level-2 products are archived at JPL’s PODAAC and at GFZs ISDC and are available 60 days after data taking. The level-2 processing software were developed independently by all three processing centres using already existing but completely independent software packages which were upgraded for GRACE specific tasks. Common data file interfaces guarantees a strong product validation. Routine processing is done at UTCSR and GFZ, while JPL only generate level-2 products at times for verification purposes.
    
+ Level 3: `Mascons`<br>
    consists of mass anomalies or other standardized products such as Monthly Ocean/Land Water Equivalent Thickness, Surface-Mass Anomaly. Similarly mass concentration blocks or `mascons` are also available.

+ Level 4: `Time Series`<br>
    Time-series of catchment level hydrological estimates of TWSA

`PySHBundle` provides the capability to obtain grided Total Water Storage Anomaly(TWSA) from Level 2 data.

## Pre-Processing of Data

### Older Approach
::: pyshbundle.read_GRACE_SH_paths
::: pyshbundle.reader_replacer_jpl.reader_replacer_jpl
::: pyshbundle.reader_replacer_csr.reader_replacer_csr
::: pyshbundle.reader_replacer_itsg.reader_replacer_itsg
::: pyshbundle.reader_replacer.reader

<br>

### Newer Approach
::: pyshbundle.io.read_jpl
::: pyshbundle.io.read_csr
::: pyshbundle.io.read_csr
::: pyshbundle.io.read_tn13
::: pyshbundle.io.read_tn14


## Computing Ling Term Mean
::: pyshbundle.load_longterm_mean

## Physical and Geodetic Constants
::: pyshbundle.GRACEconstants
