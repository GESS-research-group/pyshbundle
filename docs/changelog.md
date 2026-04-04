# Changelog

## v1.3.2 — 2026-04-04

**Code Quality (JOSS Review)**

- Applied `ruff` auto-formatting across all source files for consistent style
- Applied `ruff` linting fixes (unused imports, deprecated APIs, style violations)
- Added docstring validation via pydocstyle D-rules (missing docstrings, formatting, raw strings for backslashes)
- Configured `mypy` for static type checking; suppressed unavoidable false positives from numpy/cartopy stubs
- Added complete function type annotations (parameter and return types) to all source files: `pysh_core.py`, `shutils.py`, `hydro.py`, `reshape_SH_coefficients.py`, `io.py`, `viz_utils.py`
- Added `from __future__ import annotations` for Python 3.9 compatibility with union type syntax

**Bug Fixes**

- Fixed `sc2cs.exit(...)` and `sc2cs.sc2cs(...)` in `reshape_SH_coefficients.py` — function name was mistakenly used as a module
- Fixed `lovenrPREM.exit(...)` in `GRACEpy.py` — same issue
- Fixed missing `import rioxarray` in `hydro.py` — caused `AttributeError: 'Dataset' object has no attribute 'rio'` in `Basinaverage`
- Replaced invalid escape sequences `"\s+"` with `r"\s+"` in `io.py` (SyntaxWarning in Python 3.14)
- Replaced deprecated `np.row_stack` with `np.vstack` in `deprecated_functions.py`

**Deprecations**

- Moved `replace_zonal_coeff` and `load_longterm_mean` from `io.py` to `deprecated/deprecated_functions.py`

**CI / Docs**

- Replaced `apr-fixes` branch with `joss-review-7550` in all CI workflow triggers
- Added concurrency group to `ci-docs-gh-pages` workflow to prevent simultaneous gh-pages push failures
- Updated PyPI classifiers: added OS (Windows, Linux, macOS), Python 3.10/3.11, and Topic (Scientific/Engineering, GIS, Physics)
- Added Documentation URL to PyPI metadata
- Updated `index.md` to reflect working pip installation

---

## v1.3.1

- Bug fixes and stability improvements

## v1.3.0

- Added support for ITSG data source in `io.py`
- Added `parse_itsg_file`, `parse_itsg_header` functions
- Added degree-1 and degree-2/3 coefficient replacement utilities (`extract_C10_11_replcmnt_coeff`, `extract_C20_replcmnt_coeff`, `extract_C30_replcmnt_coeff`)

## v1.2.0

- Added `GRACE_Data_Driven_Correction_Vishwakarma` for data-driven leakage correction
- Added `Basinaverage` and `TWSCalc` for basin-scale terrestrial water storage computation
- Added visualization utilities: `polar_plot`, `mapfield`, `sc_triplot`, `cs_sqplot`

## v1.0.0

- Initial public release
- Core spherical harmonic synthesis (`gshs`) and analysis (`gsha`) functions
- SH coefficient format converters: `clm2sc`, `clm2cs`, `cs2sc`, `sc2cs`, `klm2sc`
- Associated Legendre functions: `plm`, `iplm`
- Gaussian filter: `Gaussian`, `gfilter`
