# Python Numerical Regression Fixture

Generate `numerical-regression.json` only after validating the current Python
package outputs against the relevant replication files and cross-language
expectations:

```sh
python scripts/replication-python.py --write-fixture
```

Future package speedups should preserve this fixture unless a deliberate
methodological change is documented.
