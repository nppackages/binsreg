import warnings
from types import SimpleNamespace

import numpy as np
import pytest

from binsreg import binsglm, binspwc, binsqreg, binsreg, binsregselect, binstest
from binsreg.funs import binsreg_pred


def _sample(seed=20260821, n=360):
    rng = np.random.default_rng(seed)
    x = rng.uniform(-1, 1, size=n)
    w = rng.normal(size=(n, 2))
    index = -0.3 + 0.8 * x + 0.25 * x**2 + 0.35 * w[:, 0] - 0.2 * w[:, 1]
    y = index + rng.normal(scale=0.4, size=n)
    probability = 1 / (1 + np.exp(-index))
    binary = rng.binomial(1, probability)
    by = np.arange(n) % 2
    return y, binary, x, w, by


def _assert_column_equal(left, right, column, rtol=1e-7, atol=1e-8):
    np.testing.assert_allclose(
        np.asarray(left[column], dtype=float),
        np.asarray(right[column], dtype=float),
        rtol=rtol,
        atol=atol,
        equal_nan=True,
    )


def test_binsreg_pred_treats_numeric_at_as_a_vector():
    model = SimpleNamespace(params=np.array([1.0, 2.0, 3.0]))
    basis = np.array([[1.0], [2.0]])

    flat_fit = binsreg_pred(basis, model, wvec=np.array([10.0, 20.0]))[0]
    column_fit = binsreg_pred(
        basis, model, wvec=np.array([10.0, 20.0]).reshape(-1, 1)
    )[0]

    np.testing.assert_allclose(flat_fit, [81.0, 82.0])
    np.testing.assert_allclose(column_fit, flat_fit)

    with pytest.raises(ValueError, match="Length of at"):
        binsreg_pred(basis, model, wvec=np.array([10.0]))


@pytest.mark.parametrize("deriv", [0, 1])
def test_binsglm_numeric_at_matches_centered_controls(deriv):
    _, y, x, w, _ = _sample()
    at = np.array([0.25, -0.4])
    common = {
        "dist": "Binomial",
        "link": "Logit",
        "deriv": deriv,
        "nbins": 6,
        "dots": (1, 1),
        "line": (2, 2),
        "linegrid": 3,
        "ci": None,
        "cb": None,
        "polyreg": 2,
        "polyreggrid": 3,
        "polyregcigrid": 2,
        "noplot": True,
        "vce": "HC1",
        "masspoints": "off",
    }

    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        direct = binsglm(y, x, w, at=at, **common).data_plot[0]
        centered = binsglm(y, x, w - at, at="zero", **common).data_plot[0]

    _assert_column_equal(direct.dots, centered.dots, "fit")
    _assert_column_equal(direct.line, centered.line, "fit")
    _assert_column_equal(direct.poly, centered.poly, "fit")
    _assert_column_equal(direct.polyci, centered.polyci, "polyci_l")
    _assert_column_equal(direct.polyci, centered.polyci, "polyci_r")


def test_binsqreg_numeric_at_matches_centered_controls():
    y, _, x, w, _ = _sample(n=240)
    at = np.array([0.25, -0.4])
    common = {
        "nbins": 6,
        "dots": (1, 1),
        "line": (1, 1),
        "linegrid": 3,
        "ci": None,
        "cb": None,
        "polyreg": None,
        "noplot": True,
        "masspoints": "off",
    }

    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        direct = binsqreg(y, x, w, at=at, **common).data_plot[0]
        centered = binsqreg(y, x, w - at, at="zero", **common).data_plot[0]

    _assert_column_equal(direct.dots, centered.dots, "fit", rtol=3e-6, atol=1e-6)
    _assert_column_equal(direct.line, centered.line, "fit", rtol=3e-6, atol=1e-6)


def test_binstest_numeric_at_matches_centered_controls():
    y, _, x, w, _ = _sample()
    at = np.array([0.25, -0.4])
    common = {
        "testmodelpoly": 2,
        "nbins": 10,
        "testmodel": (2, 2),
        "nsims": 39,
        "simsgrid": 20,
        "simsseed": 21,
        "masspoints": "off",
    }

    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        direct = binstest(y, x, w, at=at, **common)
        centered = binstest(y, x, w - at, at="zero", **common)

    np.testing.assert_allclose(direct.testpoly.stat, centered.testpoly.stat)
    np.testing.assert_allclose(direct.testpoly.pval, centered.testpoly.pval)


def test_binspwc_numeric_at_matches_centered_controls():
    y, _, x, w, by = _sample()
    at = np.array([0.25, -0.4])
    common = {
        "by": by,
        "bynbins": 6,
        "pselect": 1,
        "sselect": 1,
        "pwc": (2, 2),
        "nsims": 39,
        "simsgrid": 20,
        "simsseed": 21,
        "masspoints": "off",
    }

    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        direct = binspwc(y, x, w, at=at, **common)
        centered = binspwc(y, x, w - at, at="zero", **common)

    np.testing.assert_allclose(direct.tstat, centered.tstat)
    np.testing.assert_allclose(direct.pval, centered.pval)


@pytest.mark.parametrize("command", [binsreg, binsglm])
def test_small_sample_ci_with_multiple_controls_has_one_se_per_dot(command):
    rng = np.random.default_rng(20260822)
    n = 120
    x = np.linspace(0, 1, n)
    w = rng.normal(size=(n, 2))
    y = 1 + x + 0.4 * w[:, 0] - 0.2 * w[:, 1] + rng.normal(scale=0.1, size=n)
    cluster = np.arange(n) % 20

    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        result = command(
            y,
            x,
            w,
            at=np.array([0.2, -0.1]),
            dots=(0, 0),
            ci=(0, 0),
            cluster=cluster,
            noplot=True,
            masspoints="off",
        ).data_plot[0]

    assert len(result.dots) == len(result.ci)
    assert np.isfinite(result.ci["ci_l"]).all()
    assert np.isfinite(result.ci["ci_r"]).all()


def test_derivative_confidence_band_masks_nonsmooth_knots():
    y, _, x, w, _ = _sample()

    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        cb = binsreg(
            y,
            x,
            w,
            at=np.array([0.25, -0.4]),
            deriv=1,
            nbins=6,
            dots=(1, 1),
            cb=(1, 1),
            cbgrid=2,
            nsims=39,
            simsgrid=20,
            simsseed=21,
            noplot=True,
            masspoints="off",
        ).data_plot[0].cb

    knots = cb["isknot"] == 1
    assert knots.any()
    assert cb.loc[knots, "cb_l"].isna().all()
    assert cb.loc[knots, "cb_r"].isna().all()


def test_single_pselect_is_also_used_for_smoothness_selection():
    y, _, x, _, _ = _sample()

    result = binsregselect(
        y,
        x,
        pselect=[1],
        nbins=True,
        binsmethod="rot",
        masspoints="off",
    )

    assert result.prot_regul == 1
    assert result.srot_regul == 1
