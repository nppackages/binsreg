import numpy as np

from binsreg import binsglm


def test_binsglm_nolink_ci_cb_runs_for_binomial_logit():
    rng = np.random.default_rng(20260521)
    n = 400
    x = rng.uniform(size=n)
    w = rng.normal(size=(n, 1))
    p = 1 / (1 + np.exp(-(-0.25 + 0.8 * x + 0.3 * w[:, 0])))
    y = rng.binomial(1, p)

    out = binsglm(
        y,
        x,
        w,
        dist="Binomial",
        link="Logit",
        nbins=8,
        dots=(0, 0),
        line=(1, 1),
        ci=(1, 1),
        cb=(1, 1),
        nsims=20,
        simsseed=20260521,
        polyreg=None,
        noplot=True,
        nolink=True,
        vce="HC1",
        masspoints="off",
    )

    ci = out.data_plot[0].ci
    cb = out.data_plot[0].cb
    assert ci is not None
    assert cb is not None
    assert np.isfinite(ci["ci_l"].dropna()).all()
    assert np.isfinite(ci["ci_r"].dropna()).all()
    assert np.isfinite(cb["cb_l"].dropna()).all()
    assert np.isfinite(cb["cb_r"].dropna()).all()
