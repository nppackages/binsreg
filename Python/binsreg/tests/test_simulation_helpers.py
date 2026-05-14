from __future__ import annotations

import importlib.util

import pytest


np = pytest.importorskip("numpy")


def reference_binsreg_pval(num, denom, rep, tstat=None, side=None, alpha=95, lp=np.inf):
    tvec = np.empty(rep)
    pval = np.nan
    if tstat is not None:
        pval = np.zeros(tstat.shape[0])
    cval = np.nan
    k = num.shape[1]

    for i in range(rep):
        eps = np.random.normal(size=k).reshape(-1, 1)
        tx = np.matmul(num, eps) / denom.reshape(-1, 1)

        if side is not None:
            if side == "two":
                if np.isinf(lp):
                    tvec[i] = np.max(np.abs(tx))
                else:
                    tvec[i] = (np.mean(np.abs(tx) ** lp)) ** (1 / lp)
            elif side == "left":
                tvec[i] = np.max(tx)
            elif side == "right":
                tvec[i] = np.min(tx)

        if tstat is not None:
            for j in range(tstat.shape[0]):
                if tstat[j, 1] == 1:
                    pval[j] += np.max(tx) >= tstat[j, 0]
                elif tstat[j, 1] == 2:
                    pval[j] += np.min(tx) <= tstat[j, 0]
                elif tstat[j, 1] == 3:
                    if np.isinf(lp):
                        pval[j] += np.max(np.abs(tx)) >= tstat[j, 0]
                    else:
                        pval[j] += (np.mean(np.abs(tx) ** lp) ** (1 / lp)) >= tstat[j, 0]

    if tstat is not None:
        pval = pval / rep
    if side is not None:
        cval = np.quantile(tvec, alpha / 100)
    return pval, cval


def reference_binspwc_pval(nummat1, nummat2, denom1, denom2, rep, tstat=None, testtype=None, lp=np.inf):
    pval = 0
    k1 = nummat1.shape[1]
    k2 = nummat2.shape[1]
    denom1 = denom1.reshape(len(denom1), -1)
    denom2 = denom2.reshape(len(denom2), -1)

    for _ in range(rep):
        eps1 = np.random.normal(size=k1).reshape(-1, 1)
        eps2 = np.random.normal(size=k2).reshape(-1, 1)
        tx = (np.matmul(nummat1, eps1) - np.matmul(nummat2, eps2)) / np.sqrt(denom1**2 + denom2**2)

        if testtype == "left":
            pval += np.max(tx) >= tstat
        elif testtype == "right":
            pval += np.min(tx) <= tstat
        else:
            if not np.isfinite(lp):
                pval += np.max(np.abs(tx)) >= tstat
            else:
                pval += np.mean(np.abs(tx) ** lp) ** (1 / lp) >= tstat

    return pval / rep


def test_binsreg_pval_matches_reference_loop() -> None:
    if importlib.util.find_spec("binsreg") is None:
        pytest.skip("binsreg package is not importable")

    from binsreg.funs import binsreg_pval

    rng = np.random.default_rng(20240513)
    num = rng.normal(size=(17, 5))
    denom = 0.75 + rng.uniform(size=17)
    tstat = np.array([[0.5, 1], [-0.25, 2], [1.5, 3]], dtype=float)

    np.random.seed(98765)
    expected_pval, expected_cval = reference_binsreg_pval(
        num, denom, rep=200, tstat=tstat, side="two", alpha=95, lp=np.inf
    )
    np.random.seed(98765)
    actual_pval, actual_cval = binsreg_pval(
        num, denom, rep=200, tstat=tstat, side="two", alpha=95, lp=np.inf
    )

    np.testing.assert_allclose(actual_pval, expected_pval, rtol=0, atol=0)
    np.testing.assert_allclose(actual_cval, expected_cval, rtol=0, atol=1e-14)


def test_binspwc_pval_matches_reference_loop() -> None:
    if importlib.util.find_spec("binsreg") is None:
        pytest.skip("binsreg package is not importable")

    from binsreg.funs import binspwc_pval

    rng = np.random.default_rng(20240513)
    nummat1 = rng.normal(size=(13, 4))
    nummat2 = rng.normal(size=(13, 3))
    denom1 = 0.75 + rng.uniform(size=13)
    denom2 = 0.75 + rng.uniform(size=13)

    for testtype, tstat, lp in (
        ("left", 0.5, np.inf),
        ("right", -0.5, np.inf),
        ("two-sided", 1.25, np.inf),
        ("two-sided", 1.1, 2),
    ):
        np.random.seed(98765)
        expected = reference_binspwc_pval(
            nummat1, nummat2, denom1, denom2, rep=200, tstat=tstat, testtype=testtype, lp=lp
        )
        np.random.seed(98765)
        actual = binspwc_pval(
            nummat1, nummat2, denom1, denom2, rep=200, tstat=tstat, testtype=testtype, lp=lp
        )

        np.testing.assert_allclose(actual, expected, rtol=0, atol=0)


def test_fast_qreg_fit_matches_statsmodels() -> None:
    if importlib.util.find_spec("binsreg") is None:
        pytest.skip("binsreg package is not importable")

    sm = pytest.importorskip("statsmodels.api")
    from binsreg.funs import binsreg_qreg_fit

    rng = np.random.default_rng(20240513)
    x = np.column_stack((np.ones(80), rng.normal(size=(80, 4))))
    beta = np.array([0.25, -0.5, 0.75, 0.1, -0.2])
    y = x @ beta + rng.standard_t(df=5, size=80)

    actual = binsreg_qreg_fit(y, x, q=0.5)
    expected = sm.QuantReg(y, x).fit(q=0.5)

    assert actual.iterations == expected.iterations
    np.testing.assert_allclose(actual.params, expected.params, rtol=0, atol=1e-10)
    np.testing.assert_allclose(actual.fittedvalues, expected.fittedvalues, rtol=0, atol=1e-10)
    np.testing.assert_allclose(actual.cov_params(), expected.cov_params(), rtol=0, atol=1e-10)
