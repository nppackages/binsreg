from __future__ import annotations

import importlib.util
import subprocess
import sys
from pathlib import Path

import pytest


REQUIRED_MODULES = ("numpy", "pandas", "scipy", "statsmodels", "plotnine")


def repo_root() -> Path:
    path = Path(__file__).resolve()
    for candidate in (path, *path.parents):
        if (candidate / "Python" / "binsreg" / "setup.cfg").is_file():
            return candidate
    raise RuntimeError("Could not find repository root.")


@pytest.mark.parametrize("module", REQUIRED_MODULES)
def test_required_dependency_available(module: str) -> None:
    if importlib.util.find_spec(module) is None:
        pytest.skip(f"{module} is not installed in this Python environment")


def test_numerical_regression_fixture() -> None:
    missing = [module for module in REQUIRED_MODULES if importlib.util.find_spec(module) is None]
    if missing:
        pytest.skip(f"Python scientific dependencies are not installed: {', '.join(missing)}")

    root = repo_root()
    fixture = root / "Python" / "binsreg" / "tests" / "fixtures" / "numerical-regression.json"
    if not fixture.is_file():
        pytest.skip("Python numerical fixture has not been generated yet")

    completed = subprocess.run(
        [sys.executable, "scripts/replication-python.py", "--check"],
        cwd=root,
        check=False,
    )
    assert completed.returncode == 0


def test_plot_data_report_raw_bin_counts() -> None:
    missing = [module for module in REQUIRED_MODULES if importlib.util.find_spec(module) is None]
    if missing:
        pytest.skip(f"Python scientific dependencies are not installed: {', '.join(missing)}")

    import numpy as np

    root = repo_root()
    src = root / "Python" / "binsreg" / "src"
    if str(src) not in sys.path:
        sys.path.insert(0, str(src))

    from binsreg import binsglm, binsqreg, binsreg

    x = np.linspace(0, 1, 90)
    y = 1 + 2 * x + 0.1 * np.sin(np.arange(x.size))
    expected = np.array([30, 30, 30])

    def check_counts(result) -> None:
        data = result.data_plot[0]
        np.testing.assert_array_equal(data.data_bin["n"].to_numpy(), expected)
        np.testing.assert_array_equal(
            data.dots["n"].to_numpy(),
            expected[data.dots["bin"].astype(int).to_numpy()],
        )

    common = dict(
        nbins=3,
        binspos="es",
        dots=(0, 0),
        line=(0, 0),
        ci=None,
        cb=None,
        masspoints="off",
        noplot=True,
    )
    check_counts(binsreg(y, x, **common))
    check_counts(binsqreg(y, x, **common))
    check_counts(binsglm(y, x, **common))


def test_noplot_keeps_requested_polynomial_data() -> None:
    missing = [module for module in REQUIRED_MODULES if importlib.util.find_spec(module) is None]
    if missing:
        pytest.skip(f"Python scientific dependencies are not installed: {', '.join(missing)}")

    import numpy as np

    root = repo_root()
    src = root / "Python" / "binsreg" / "src"
    if str(src) not in sys.path:
        sys.path.insert(0, str(src))

    from binsreg import binsglm, binsqreg, binsreg

    x = np.linspace(0, 1, 90)
    y = 1 + 2 * x + 0.1 * np.sin(np.arange(x.size))

    common = dict(
        nbins=3,
        binspos="es",
        dots=(0, 0),
        line=None,
        ci=None,
        cb=None,
        polyreg=2,
        polyreggrid=4,
        masspoints="off",
        noplot=True,
    )

    for fit in (binsreg(y, x, **common), binsqreg(y, x, **common), binsglm(y, x, **common)):
        assert fit.bins_plot is None
        assert fit.data_plot[0].poly is not None
        assert not fit.data_plot[0].poly.empty
        assert {"x", "fit", "n"}.issubset(fit.data_plot[0].poly.columns)
