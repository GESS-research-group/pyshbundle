# Core functionality
The core functionality is processing spherical harmonics coefficients to mass change fields.

## Spherical Harmonic Analysis and Synthesis
`gsha(f, method: str, grid: str = None, lmax: int = -9999):`
::: pyshbundle.pysh_core.gsha
`gshs(field, quant = 'none', grd = 'mesh', n = -9999, h = 0, jflag = 1):`
::: pyshbundle.pysh_core.gshs
`PhaseCalc(fts, ffts):`
::: pyshbundle.pysh_core.PhaseCalc

## Intro to Grace Data Driven Correction
`GRACE_Data_Driven_Correction_Vishwakarma(F, cf, GaussianR, basins):`
::: pyshbundle.pysh_core.GRACE_Data_Driven_Correction_Vishwakarma

## Hydrological Applications with GRACE
`Basinaverage(temp, gs, shp_basin, basin_area):`
::: pyshbundle.hydro.Basinaverage
`TWSCalc(data, lmax: int, gs: float, r:float, m: int):`
::: pyshbundle.hydro.TWSCalc