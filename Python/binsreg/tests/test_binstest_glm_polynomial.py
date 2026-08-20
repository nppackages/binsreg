import warnings

import numpy as np

from binsreg import binstest


def test_glm_polynomial_tests_compare_fits_on_requested_scale():
    rng = np.random.default_rng(20260820)
    n = 2000
    x = rng.uniform(-1, 1, size=n)
    w = rng.normal(size=(n, 2))
    index = -0.4 + 1.1 * x + 0.35 * w[:, 0] - 0.2 * w[:, 1]
    probability = 1 / (1 + np.exp(-index))
    y = rng.binomial(1, probability)

    common = {
        "w": w,
        "dist": "Binomial",
        "link": "Logit",
        "testmodelpoly": 1,
        "nbins": 16,
        "testmodel": (1, 1),
        "nsims": 99,
        "simsgrid": 20,
        "simsseed": 17,
        "masspoints": "off",
    }

    with warnings.catch_warnings():
        warnings.filterwarnings(
            "ignore",
            message="Testing procedures are valid when nbins is much larger",
        )
        response_level = binstest(y, x, **common)
        index_level = binstest(y, x, nolink=True, **common)
        response_derivative = binstest(y, x, deriv=1, **common)

    response_level_stat = float(np.asarray(response_level.testpoly.stat).squeeze())
    index_level_stat = float(np.asarray(index_level.testpoly.stat).squeeze())
    response_derivative_stat = float(
        np.asarray(response_derivative.testpoly.stat).squeeze()
    )

    assert np.isfinite(response_level_stat)
    assert np.isfinite(response_derivative_stat)
    assert response_level_stat < 10
    assert response_derivative_stat < 10
    assert abs(response_level_stat - index_level_stat) < 2
