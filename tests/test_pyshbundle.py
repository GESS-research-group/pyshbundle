# tests/test_validation.py

import pytest
import numpy as np
from .test_validation_pyshbundle import validation_pyshbundle, load_matlab_reference

# ── Constants ─────────────────────────────────────────────────────────────────

RMSE_THRESHOLD = 1e-3
NRMSE_THRESHOLD = 1e-5

# ── Fixtures ──────────────────────────────────────────────────────────────────


@pytest.fixture(scope="module")
def tws_results(tmp_path_factory):
    """Compute TWS once per test session.

    scope="module" avoids re-running the expensive computation for each test.

    NOTE: path_sh points to external data. In CI this requires either:
      - The data committed to the repo (not recommended for large datasets)
      - A data download step in the GitHub Actions workflow before tests run
    """
    import os

    PROJECT_ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
    base = os.environ.get("PYSHBUNDLE_DATA_DIR", os.path.join(PROJECT_ROOT, "data"))
    path_sh = os.path.join(base, "JPL_input")
    path_tn14 = os.path.join(
        PROJECT_ROOT, "pyshbundle", "data", "JPL_TN_files", "TN-14_C30_C20_GSFC_SLR.txt"
    )
    path_tn13 = os.path.join(
        PROJECT_ROOT, "pyshbundle", "data", "JPL_TN_files", "TN-13_GEOC_JPL_RL06.txt"
    )

    if not os.path.exists(path_sh):
        pytest.fail(
            f"JPL input data not found at {path_sh}. "
            "Set PYSHBUNDLE_DATA_DIR env var or provide data/JPL_input/."
        )

    tws_computed = validation_pyshbundle(path_sh, path_tn13, path_tn14, source="jpl")
    tws_reference = load_matlab_reference()
    return tws_computed, tws_reference


# ── Tests ─────────────────────────────────────────────────────────────────────


def test_tws_output_shape(tws_results):
    """Computed and reference TWS fields must have identical shape."""
    tws_computed, tws_reference = tws_results
    assert tws_computed.shape == tws_reference.shape, (
        f"Shape mismatch: computed {tws_computed.shape} "
        f"vs reference {tws_reference.shape}"
    )


def test_tws_output_dtype(tws_results):
    """Output must be float32 as specified."""
    tws_computed, _ = tws_results
    assert tws_computed.dtype == np.float32


def test_gridwise_rmse(tws_results):
    """Gridwise RMSE between computed and reference TWS must be below 1e-3.

    This is the primary accuracy criterion.
    """
    tws_computed, tws_reference = tws_results
    # tws_computed, tws_reference = tws_computed[0:100], tws_reference[0:100]
    diff = tws_reference - tws_computed
    gridwise_rmse = np.sqrt(np.mean(diff**2, axis=0))

    # Report the worst offending grid point on failure
    max_rmse = np.nanmax(gridwise_rmse)
    assert np.all(gridwise_rmse < RMSE_THRESHOLD), (
        f"Gridwise RMSE exceeded threshold. "
        f"Max RMSE = {max_rmse:.6e}, threshold = {RMSE_THRESHOLD:.6e}"
    )


def test_gridwise_nrmse(tws_results):
    """Gridwise NRMSE must be below 1e-5.

    Normalised by the std of the reference field along the time axis.
    """
    tws_computed, tws_reference = tws_results
    # tws_computed, tws_reference = tws_computed[0:100], tws_reference[0:100]
    diff = tws_reference - tws_computed
    gridwise_rmse = np.sqrt(np.mean(diff**2, axis=0))
    std_ref = np.std(tws_reference, axis=0)

    # Guard against division by zero at invariant grid points
    with np.errstate(invalid="ignore", divide="ignore"):
        gridwise_nrmse = np.where(std_ref > 0, gridwise_rmse / std_ref, 0.0)

    max_nrmse = np.nanmax(gridwise_nrmse)
    assert np.all(gridwise_nrmse < NRMSE_THRESHOLD), (
        f"Gridwise NRMSE exceeded threshold. "
        f"Max NRMSE = {max_nrmse:.6e}, threshold = {NRMSE_THRESHOLD:.6e}"
    )


def test_no_nan_in_output(tws_results):
    """TWS output must contain no NaN values."""
    tws_computed, _ = tws_results
    assert not np.any(np.isnan(tws_computed)), "NaN values found in computed TWS output"


def test_no_nan_in_reference(tws_results):
    """Reference data sanity check — if this fails, the .mat file is corrupt."""
    _, tws_reference = tws_results
    assert not np.any(np.isnan(tws_reference)), (
        "NaN values found in MATLAB reference data"
    )
