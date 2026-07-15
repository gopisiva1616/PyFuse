import os
from datetime import datetime
from pathlib import Path

import pytest

# Keep numerical libs from oversubscribing threads in constrained test environments.
# This mirrors the runtime guard in cli.py but applies to direct pytest imports.
for _var in (
    "OPENBLAS_NUM_THREADS",
    "OMP_NUM_THREADS",
    "MKL_NUM_THREADS",
    "NUMEXPR_NUM_THREADS",
    "VECLIB_MAXIMUM_THREADS",
    "BLIS_NUM_THREADS",
):
    os.environ.setdefault(_var, "1")

from pyfuse.utils.common_utils import config, utils
from pyfuse.utils.resource_manger import ResourceManager

DEFAULT_PYTEST_INPUT = "input_27cases.txt"
DEFAULT_PYTEST_TRUTH_OUTPUT = "expected_output_27cases.xlsx"
DEFAULT_PYTEST_INPUT_FORMAT = "star"


@pytest.fixture(scope="session")
def fixture_root() -> Path:
    return Path(__file__).resolve().parent / "fixtures"


@pytest.fixture(scope="session")
def input_bkpt_path(fixture_root: Path) -> str:
    env_input = os.environ.get("INPUT_BKPT")
    if env_input:
        return os.path.abspath(env_input)
    return str((fixture_root / "input" / DEFAULT_PYTEST_INPUT).resolve())


@pytest.fixture(scope="session")
def truth_output_path(fixture_root: Path) -> str:
    env_truth = os.environ.get("TRUTH_OUTPUT")
    if env_truth:
        return os.path.abspath(env_truth)
    return str((fixture_root / "truth_output" / DEFAULT_PYTEST_TRUTH_OUTPUT).resolve())


@pytest.fixture(scope="session")
def input_format() -> str:
    return os.environ.get("INPUT_FORMAT", DEFAULT_PYTEST_INPUT_FORMAT)


@pytest.fixture(scope="session")
def reference_genome() -> str | None:
    ref = os.environ.get("REFERENCE_GENOME")
    return os.path.abspath(ref) if ref else None


@pytest.fixture(scope="session")
def output_dir(tmp_path_factory: pytest.TempPathFactory) -> str:
    env_out = os.environ.get("OUTPUT_PATH")
    if env_out:
        base = Path(os.path.abspath(env_out))
        out = base / f"pytest_results_{datetime.now().strftime('%d-%m-%Y_%I-%M-%S')}"
        out.mkdir(parents=True, exist_ok=True)
        return str(out)

    out = tmp_path_factory.mktemp("pyfuse_pytest_output")
    return str(out)


@pytest.fixture(scope="session")
def configured_runtime_resource_dir(output_dir: str) -> str:
    """Resolve config assets/resources and apply runtime config used by RunPyfuse."""
    res_manager = ResourceManager(cache_dir=os.environ.get("CACHE_DIR"))

    cfg_assets = res_manager.resolve_config_assets()
    config.update({k: str(v) for k, v in cfg_assets.items()})

    env_resource = os.environ.get("RESOURCE_PATH")
    if env_resource:
        bundle = res_manager.resolve(user_dir=os.path.abspath(env_resource))
    else:
        genome = os.environ.get("GENOME", "grch37").lower()
        default_bundle = f"default_{genome}"
        try:
            bundle = res_manager.resolve_from_cache(default_bundle=default_bundle)
        except FileNotFoundError as exc:
            pytest.skip(
                "No installed default resources found for test run. "
                "Set RESOURCE_PATH to a bundle directory or install one via "
                f"'pyfuse resources install --genome {genome}'. Details: {exc}"
            )

    config.update({k: str(v) for k, v in bundle.resource_files.items()})
    utils.out_path = output_dir
    utils.res_path = str(bundle.root_dir)
    return str(bundle.root_dir)
