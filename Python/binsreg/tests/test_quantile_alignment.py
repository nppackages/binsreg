import numpy as np

from binsreg.funs import binsreg_quantile, genKnot_qs


def test_binsreg_quantile_matches_stata_r_type2_values():
    x = np.array([1.0, 2.0, 3.0, 4.0])
    q = np.array([0.0, 0.10, 0.25, 0.30, 0.50, 0.75, 1.0])

    observed = binsreg_quantile(x, q)

    np.testing.assert_allclose(observed, [1.0, 1.0, 1.5, 2.0, 2.5, 3.5, 4.0])


def test_binsreg_quantile_axis_zero_matches_columnwise_type2_median():
    x = np.array(
        [
            [1.0, 10.0],
            [2.0, 20.0],
            [3.0, 30.0],
            [4.0, 40.0],
        ]
    )

    observed = binsreg_quantile(x, 0.5, axis=0)

    np.testing.assert_allclose(observed, [2.5, 25.0])


def test_genknot_qs_uses_type2_quantiles_with_ties():
    x = np.array([0.0, 0.0, 1.0, 2.0, 3.0, 100.0])

    observed = genKnot_qs(x, 4)

    np.testing.assert_allclose(observed, [0.0, 0.0, 1.5, 3.0, 100.0])
