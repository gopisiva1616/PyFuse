# Testing

PyFuse includes a helper command `pyfuse-test` to run pytest with HTML report output by default.

## Quick start

```bash
pyfuse-test
```

This runs pytest with verbose summary and creates a timestamped self-contained HTML report.

To print per-case output for passed integration comparisons too:

```bash
pyfuse-test --show-passed-details -m integration
```

## Accepted environment variables for integration tests

The integration test runtime reads these environment variables:

- `INPUT_BKPT`: path to input breakpoint file
- `TRUTH_OUTPUT`: path to truth output file
- `INPUT_FORMAT`: input format parser (for example `star`, `default`, `arriba`)
- `REFERENCE_GENOME`: FASTA path used for sequence-related annotation tests
- `OUTPUT_PATH`: base directory for test result output folders
- `CACHE_DIR`: cache directory used by resource manager
- `RESOURCE_PATH`: explicit resource bundle directory for tests
- `GENOME`: default genome family (`grch37` or `grch38`) when `RESOURCE_PATH` is not set
- `FUSION2VCF`: toggles alternate truth parsing mode when set to `True` or `TRUE`
- `perl`: enables legacy normalization branch in truth parser logic

## Examples

Run integration tests with defaults:

```bash
pyfuse-test -m integration
```

Run integration tests with custom input format and reference:

```bash
INPUT_FORMAT=star REFERENCE_GENOME=/data/ref/genome.fa pyfuse-test -m integration
```

Pin a specific resource bundle:

```bash
RESOURCE_PATH=/data/resources/custom_bundle pyfuse-test -m integration
```

Use a custom report name:

```bash
pyfuse-test --report pyfuse_pytest_report_custom.html
```

Pass any additional pytest options:

```bash
pyfuse-test -m integration -k annotation --maxfail=1
```

The `--show-passed-details` option sets `PYFUSE_SHOW_MATCH_DETAILS=1` and uses `--capture=tee-sys` with INFO logging so pass comparison messages are visible in terminal output and preserved in pytest HTML output.

## Troubleshooting

If `pyfuse-test` is not found, your current environment likely has an older install without the new entry point.

Reinstall in the active environment:

```bash
pip install -e .[dev]
```

Or use module invocation directly:

```bash
python -m pyfuse.test_cli -m integration
```
