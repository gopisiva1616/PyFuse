"""CLI helper for running PyFuse tests with HTML reporting."""

from __future__ import annotations

import argparse
import os
import sys
from datetime import datetime


def main() -> int:
    parser = argparse.ArgumentParser(
        prog="pyfuse-test",
        description="Run pytest for PyFuse and generate an HTML report by default.",
        epilog=(
            "Environment variables consumed by integration tests:\n"
            "  INPUT_BKPT, TRUTH_OUTPUT, INPUT_FORMAT, REFERENCE_GENOME, OUTPUT_PATH,\n"
            "  CACHE_DIR, RESOURCE_PATH, GENOME, FUSION2VCF, perl\n\n"
            "Examples:\n"
            "  pyfuse-test -m integration\n"
            "  pyfuse-test --show-passed-details -m integration\n"
            "  INPUT_FORMAT=star REFERENCE_GENOME=/data/ref.fa pyfuse-test -m integration\n"
            "  RESOURCE_PATH=/data/resources/custom_bundle pyfuse-test -m integration\n"
        ),
        formatter_class=argparse.RawTextHelpFormatter,
    )
    parser.add_argument(
        "--report",
        default=None,
        help="output HTML report path (default: timestamped pyfuse report filename)",
    )
    parser.add_argument(
        "--no-html",
        action="store_true",
        help="disable HTML report generation",
    )
    parser.add_argument(
        "--show-passed-details",
        action="store_true",
        help="print pass/match details for each case and include them in HTML output",
    )

    args, passthrough = parser.parse_known_args()

    try:
        import pytest
    except ImportError:
        print("ERROR: pytest is not installed. Install dev dependencies: pip install -e '.[dev]'", file=sys.stderr)
        return 2

    pytest_args = ["-vrs"]

    if args.show_passed_details:
        os.environ["PYFUSE_SHOW_MATCH_DETAILS"] = "1"
        if all(opt not in passthrough for opt in ["-s", "--capture=no", "--capture=tee-sys"]):
            pytest_args.append("--capture=tee-sys")
        if "--log-cli-level" not in " ".join(passthrough):
            pytest_args.extend(["--log-cli-level=INFO", "--log-level=INFO"])

    if not args.no_html:
        report_name = args.report or f"pyfuse_pytest_report_{datetime.now().strftime('%b%d-%Y_%H%M%S')}.html"
        pytest_args.extend([f"--html={report_name}", "--self-contained-html"])

    pytest_args.extend(passthrough)
    return int(pytest.main(pytest_args))


if __name__ == "__main__":
    raise SystemExit(main())
