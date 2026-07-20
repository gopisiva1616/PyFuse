


# src/pyfuse/cli.py
import os
import sys
import time
import shutil
import logging
import argparse
import yaml
from tempfile import mkdtemp
from pathlib import Path
from datetime import datetime

# Conservative thread defaults to avoid OpenBLAS/MKL process-limit failures
# in constrained environments. Users can override by setting env vars explicitly.
for _var in (
    "OPENBLAS_NUM_THREADS",
    "OMP_NUM_THREADS",
    "MKL_NUM_THREADS",
    "NUMEXPR_NUM_THREADS",
    "VECLIB_MAXIMUM_THREADS",
    "BLIS_NUM_THREADS",
):
    os.environ.setdefault(_var, "1")

from pyfuse import setup_logger, __version__
from pyfuse.utils.resource_manger import ResourceManager

logger = logging.getLogger()

def _normalized_user_command() -> str:
    cmd = ["pyfuse"]
    cmd.extend(sys.argv[1:])
    return " ".join(cmd)


def _resolved_arg_value(key: str, value):
    path_keys = {"settings", "fusion_bkpt", "out_path", "resource", "target_bed", "ref"}
    if value is None:
        return "None"
    if key in path_keys:
        return str(Path(str(value)).expanduser().resolve())
    return str(value)


def _build_run_manifest(args, resource_path: str, resource_source: str) -> str:
    arg_specs = [
        ("fusion_bkpt", ("-i", "--fusion-breakpoints", "--fus_bkpt")),
        ("out_path", ("-o",)),
        ("format", ("--input_format",)),
        ("genome", ("--genome",)),
        ("resource", ("-r", "--resource_path", "--res_path")),
        ("target_bed", ("-t", "--target_bed")),
        ("ref", ("-g", "--reference")),
        ("settings", ("-s", "--settings")),
        ("debug", ("-d",)),
        ("out_name", ("-n",)),
    ]

    arg_lines = []
    for name, flags in arg_specs:
        if _arg_was_provided(*flags):
            arg_lines.append(f"  - {name}: {_resolved_arg_value(name, getattr(args, name, None))}")

    if not arg_lines:
        arg_lines = ["  - (none)"]

    return (
        "\n====================================================\n"
        "PyFuse: Python Fusion Annotator\n"
        f"version: v{__version__}\n"
        f"Run timestamp: {datetime.now()}\n\n"
        f"Users command:  {_normalized_user_command()}\n\n"
        f"Resource path ({resource_source}): {resource_path}\n\n"
        "Input user arguments:\n"
        + "\n".join(arg_lines)
        + "\n====================================================\n"
    )


class RequiredHelpFormatter(argparse.ArgumentDefaultsHelpFormatter, argparse.RawTextHelpFormatter):
    """Append a required marker and avoid noisy '(default: None)' in help output."""

    def _get_help_string(self, action):
        help_text = action.help or ""

        if '%(default)' not in help_text and action.default not in (None, argparse.SUPPRESS):
            help_text = f"{help_text} (default: {action.default})"

        if action.required and action.option_strings and '[required]' not in help_text.lower():
            help_text = f"[required] {help_text}"

        return help_text


def _ensure_cmd_logger(debug: bool = False) -> None:
    level = logging.DEBUG if debug else logging.INFO
    if not logger.handlers:
        logging.basicConfig(level=level, format="%(asctime)s %(levelname)s [%(module)s]: %(message)s")
    logger.setLevel(level)


def _load_settings_file(settings_path: str | None) -> dict[str, object]:
    if not settings_path:
        return {}
    path = Path(settings_path).expanduser()
    if not path.is_file():
        return {}
    data = yaml.safe_load(path.read_text(encoding="utf-8")) or {}
    return data if isinstance(data, dict) else {}


def _extract_settings_path_from_argv() -> str | None:
    argv = sys.argv[1:]
    for idx, token in enumerate(argv):
        if token in {"-s", "--settings"} and idx + 1 < len(argv):
            return argv[idx + 1]
        if token.startswith("--settings="):
            return token.split("=", 1)[1]
    return None


def _load_help_settings() -> dict[str, object]:
    # Help text should reflect actual values from packaged defaults and optional user settings override.
    packaged_settings = Path(__file__).resolve().parent / "config" / "settings.yaml"
    cfg: dict[str, object] = {}
    cfg.update(_load_settings_file(str(packaged_settings)))
    cfg.update(_load_settings_file(_extract_settings_path_from_argv()))
    return cfg


def _display_genome(genome: str) -> str:
    mapping = {
        "grch37": "GRCh37",
        "grch38": "GRCh38",
    }
    return mapping.get(genome, genome)


def _resolve_resource_root(args, settings_cfg: dict[str, object] | None = None):
    if getattr(args, "resource_root", None):
        return args.resource_root
    env_root = os.getenv("PYFUSE_RESOURCE_ROOT")
    if env_root:
        return env_root
    if settings_cfg and settings_cfg.get("resource_root"):
        return str(settings_cfg.get("resource_root"))
    return None


def _arg_was_provided(*flags: str) -> bool:
    argv = sys.argv[1:]
    for token in argv:
        for flag in flags:
            if token == flag or token.startswith(f"{flag}="):
                return True
    return False


def _confirm_yes_no(prompt_message: str) -> bool:
    while True:
        try:
            answer = input(prompt_message).strip().lower()
        except EOFError:
            return False
        if answer in {"y", "yes"}:
            return True
        if answer in {"n", "no"}:
            return False
        print("Please answer yes or no.")


def _resources_install(args) -> int:
    from pyfuse.prep_resources import extract_gene_info_from_gtf as prep_res

    settings_cfg = _load_settings_file(getattr(args, "settings", None))
    res_manager = ResourceManager(
        cache_dir=args.cache_dir,
        resource_root=_resolve_resource_root(args, settings_cfg),
    )
    default_bundle = f"default_{args.genome}"
    display_genome = _display_genome(args.genome)

    # Runtime precedence for install inputs:
    # CLI option > custom settings.yaml > packaged defaults inside prep script.
    gtf_key = f"default_refseq_gtf_url_{args.genome}"
    annot_key = f"default_refseq_assembly_summary_url_{args.genome}"
    gtf_input = args.gtf or settings_cfg.get(gtf_key)
    annot_input = args.annot_summary or settings_cfg.get(annot_key)
    mane_input = args.mane_file or settings_cfg.get("default_mane_url")

    if args.genome == "grch37" and mane_input:
        logger.warning(
            "Custom MANE input provided for %s. MANE is GRCh38-centered; "
            "the resource prep step will ask for explicit confirmation before using back-mapped GRCh37 MANE inputs.",
            display_genome,
        )

    existing_versions = [p.name for p in res_manager.list_cached_versions(default_bundle)]
    if existing_versions:
        latest_version = existing_versions[-1]
        proceed = _confirm_yes_no(
            "Bundle has already been set up for the selected genome "
            f"{display_genome} (latest: {latest_version}). "
            "Do you still want to create a new version? [yes/no]: "
        )
        if not proceed:
            logger.info("Skipping new resource bundle creation for %s", display_genome)
            return 0

    ts = datetime.now().strftime("%Y%m%d%H%M%S")
    out_name = f"{default_bundle}_{ts}"
    logger.info("Preparing new managed resource bundle for %s", display_genome)

    staging_root = Path(mkdtemp(prefix="pyfuse_resource_bootstrap_", dir=str(res_manager.cache_base)))
    prep_args = argparse.Namespace(
        genome=args.genome,
        source=args.source,
        gtf=str(gtf_input) if gtf_input else None,
        annot_summary=str(annot_input) if annot_input else None,
        out_path=str(staging_root),
        out_name=out_name,
        roi=None,
        annot_delim=' ',
        long_trans=False,
        req_annot_cols=['transcript_id', 'gene_id', 'gene', 'exon_number', 'exonCount'],
        mane_file=str(mane_input) if mane_input else None,
    )

    try:
        prep_res.main(prep_args)
        generated_dir = staging_root / "custom_resource_files" / out_name
        if not generated_dir.is_dir():
            raise FileNotFoundError(f"Expected generated resource directory not found: {generated_dir}")

        dest_dir = res_manager.resources_cache_base / default_bundle / f"v{ts}"
        dest_dir.parent.mkdir(parents=True, exist_ok=True)
        if dest_dir.exists():
            shutil.rmtree(dest_dir, ignore_errors=True)
        shutil.copytree(generated_dir, dest_dir)

        # Validate copied bundle via existing ResourceManager checks.
        _ = res_manager.resolve(user_dir=str(dest_dir))

        logger.info("Installed resource bundle for %s at %s", display_genome, dest_dir)
        return 0
    finally:
        shutil.rmtree(staging_root, ignore_errors=True)


def _resources_list(args) -> int:
    settings_cfg = _load_settings_file(getattr(args, "settings", None))
    res_manager = ResourceManager(
        cache_dir=args.cache_dir,
        resource_root=_resolve_resource_root(args, settings_cfg),
    )
    for bundle in ["default_grch37", "default_grch38"]:
        versions = [p.name for p in res_manager.list_cached_versions(bundle)]
        display_bundle = _display_genome(bundle.removeprefix("default_"))
        if versions:
            print(f"{display_bundle}: {', '.join(versions)}")
        else:
            print(f"{display_bundle}: (none installed)")
    return 0


def _resources_path(args) -> int:
    settings_cfg = _load_settings_file(getattr(args, "settings", None))
    res_manager = ResourceManager(
        cache_dir=args.cache_dir,
        resource_root=_resolve_resource_root(args, settings_cfg),
    )
    bundle = f"default_{args.genome}"
    display_genome = _display_genome(args.genome)
    try:
        rb = res_manager.resolve_from_cache(default_bundle=bundle)
        print(str(rb.root_dir))
        return 0
    except FileNotFoundError:
        print(f"No installed bundle found for {display_genome}. Run: pyfuse resources install --genome {args.genome}")
        return 1


def _resources_verify(args) -> int:
    settings_cfg = _load_settings_file(getattr(args, "settings", None))
    res_manager = ResourceManager(
        cache_dir=args.cache_dir,
        resource_root=_resolve_resource_root(args, settings_cfg),
    )
    bundle = f"default_{args.genome}"
    display_genome = _display_genome(args.genome)
    try:
        rb = res_manager.resolve_from_cache(default_bundle=bundle)
        print(f"OK: {display_genome} -> {rb.root_dir}")
        return 0
    except Exception as exc:
        print(f"FAILED: {display_genome} verification error: {exc}")
        return 1


def main():
    help_cfg = _load_help_settings()
    gtf37_default = str(help_cfg.get("default_refseq_gtf_url_grch37", "(not configured)"))
    gtf38_default = str(help_cfg.get("default_refseq_gtf_url_grch38", "(not configured)"))
    asm37_default = str(help_cfg.get("default_refseq_assembly_summary_url_grch37", "(not configured)"))
    asm38_default = str(help_cfg.get("default_refseq_assembly_summary_url_grch38", "(not configured)"))
    mane_default = str(help_cfg.get("default_mane_url", "(not configured)"))

    # pyfuse main options
    tool_des = (
        "\nDescription: Tool to annotate gene fusion\n"
        "\nResource setup quick start:\n"
        "  pyfuse resources install --genome grch37\n"
        "  pyfuse resources install -h   # see URL/file override options\n"
    )
    parser = argparse.ArgumentParser(description=tool_des, add_help=False, formatter_class=RequiredHelpFormatter)
    parser.add_argument('-h',
                        action=_HelpAction)
    parser.add_argument('-d', dest='debug',
                        action="store_true",
                        help="flag to enable debug log messages")
    parser.add_argument('-s', '--settings',
                        type=(os.path.abspath),
                        help='custom settings yaml file (optional)') 

    subparsers = parser.add_subparsers(help='MODES', dest="command", required=True)

    # subparser to run complete annotation program
    app_des = "Description: Mode to run entire annotation program"
    annotator_parser = subparsers.add_parser("annotator", description=app_des, formatter_class=RequiredHelpFormatter)
    annotator_parser.add_argument('-i', '--fusion-breakpoints', 
                            dest='fusion_bkpt',
                            required=True,
                            type=(os.path.abspath),
                            help='input breakpoint file')
    annotator_parser.add_argument('-n', dest='out_name',
                        help="output name")
    annotator_parser.add_argument('-o', dest='out_path',
                        type=(os.path.abspath),
                        help="output path",
                        required=True)
    annotator_parser.add_argument("--genome", choices=["grch37", "grch38"], default="grch37",
                    help="Which managed default resource bundle to use (ignored if --resource_path is provided).")
    annotator_parser.add_argument('-r', '--resource_path', '--res_path',
                            dest='resource',
                            type=(os.path.abspath),
                            help='exact PyFuse resource bundle directory to use for this run')
    annotator_parser.add_argument('-t', '--target_bed',
                            dest='target_bed',
                            type=(os.path.abspath),
                            help='bed file to filter fusion')
    annotator_parser.add_argument('--input_format',
                            required=True,
                            dest='format',
                            default='default',
                            help='format of input fusion_bkpt',
                            choices=['default', 'star', 'arriba', 'tophat',
                                     'fusion_catcher', 'longgf', 'fusion_inspector'])
    annotator_parser.add_argument('-g', '--reference',
                            dest='ref',
                            type=(os.path.abspath),
                            default=None,
                            help='reference genome in .fasta format for fusion sequence annotation')
    annotator_parser.add_argument(
                            '--fusion-plot-mode',
                            dest='fusion_plot_mode',
                            choices=['embed', 'external'],
                            default='external',
                            help='how fusion visualizations are stored in report output (default: external; use embed for one self-contained HTML report)')


    # subparser to manage local resource cache
    res_des = (
        "Description: build/install and manage local PyFuse resource bundles. "
        "Use 'pyfuse resources install -h' to see --gtf/--annot-summary/--mane-file overrides."
    )
    resources_parser = subparsers.add_parser("resources", description=res_des, formatter_class=RequiredHelpFormatter)
    resources_sub = resources_parser.add_subparsers(dest="resources_command", required=True)

    list_parser = resources_sub.add_parser("list", formatter_class=RequiredHelpFormatter,
                                           help="list cached default resource bundles")
    list_parser.add_argument("--cache-dir", default=None,
                             help="cache base directory (defaults to platform user cache)")
    list_parser.add_argument("--resource-root", default=None,
                             help="persistent root for managed resources (defaults to platform user data dir)")

    path_parser = resources_sub.add_parser("path", formatter_class=RequiredHelpFormatter,
                                           help="print latest cached path for a genome bundle")
    path_parser.add_argument("--genome", choices=["grch37", "grch38"], default="grch37")
    path_parser.add_argument("--cache-dir", default=None,
                             help="cache base directory (defaults to platform user cache)")
    path_parser.add_argument("--resource-root", default=None,
                             help="persistent root for managed resources (defaults to platform user data dir)")

    verify_parser = resources_sub.add_parser("verify", formatter_class=RequiredHelpFormatter,
                                             help="verify latest cached bundle manifest/files")
    verify_parser.add_argument("--genome", choices=["grch37", "grch38"], default="grch37")
    verify_parser.add_argument("--cache-dir", default=None,
                               help="cache base directory (defaults to platform user cache)")
    verify_parser.add_argument("--resource-root", default=None,
                               help="persistent root for managed resources (defaults to platform user data dir)")

    install_parser = resources_sub.add_parser(
        "install",
        formatter_class=RequiredHelpFormatter,
        help="bootstrap and install a bundle (supports URL/file overrides)",
        description=(
            "Create or update a managed PyFuse resource bundle for one genome.\n"
            "This is typically a one-time setup per genome, but you can create new versions anytime.\n"
            "If a bundle already exists for that genome, PyFuse prompts before creating a new version."
        ),
    )
    install_parser.add_argument("--genome", choices=["grch37", "grch38"], required=True)
    install_parser.add_argument("--source", choices=["refseq-gtf", "ucsc-gtf"], default="refseq-gtf",
                                help="source mode used for bootstrap")
    install_parser.add_argument(
        "--gtf",
        default=None,
        help=(
            "optional local path or URL to override default GTF.\n"
            "Defaults from config/settings.yaml:\n"
            f"  grch37: {gtf37_default}\n"
            f"  grch38: {gtf38_default}"
        ),
    )
    install_parser.add_argument(
        "--annot-summary",
        default=None,
        help=(
            "optional local path or URL to override default assembly summary report.\n"
            "Defaults from config/settings.yaml:\n"
            f"  grch37: {asm37_default}\n"
            f"  grch38: {asm38_default}"
        ),
    )
    install_parser.add_argument(
        "--mane-file", "--mane-url",
        dest="mane_file",
        default=None,
        help=(
            "optional MANE input override (local file path or URL).\n"
            f"Default from config/settings.yaml (default_mane_url): {mane_default}\n"
            "For grch37, this is treated as a back-mapped resource and requires explicit\n"
            "interactive confirmation during install."
        ),
    )
    install_parser.add_argument("--cache-dir", default=None,
                                help="cache base directory (defaults to platform user cache)")
    install_parser.add_argument("--resource-root", default=None,
                                help="persistent root for managed resources (defaults to platform user data dir)")
    install_parser.epilog = (
        "One-time setup and versioning:\n"
        "  - Run once per genome to create baseline bundles (default_grch37, default_grch38).\n"
        "  - Re-running install for the same genome creates a new version after confirmation.\n\n"
        "Examples:\n"
        "  pyfuse resources install --genome grch37\n"
        "  pyfuse resources install --genome grch38\n"
        "  pyfuse resources install --genome grch38 --gtf <url_or_file> --annot-summary <url_or_file>\n"
        "  pyfuse resources install --genome grch38 --mane-file <url_or_file>\n"
        "  pyfuse resources install --genome grch37 --mane-file <backmapped_mane_url_or_file>\n\n"
        "After install:\n"
        "  pyfuse resources list\n"
        "  pyfuse resources path --genome grch37\n"
        "  pyfuse resources verify --genome grch37"
    )

    args = parser.parse_args()
    stime = time.time()

    if args.command == 'resources':
        _ensure_cmd_logger(debug=getattr(args, "debug", False))
        if args.resources_command == 'list':
            return _resources_list(args)
        if args.resources_command == 'path':
            return _resources_path(args)
        if args.resources_command == 'verify':
            return _resources_verify(args)
        if args.resources_command == 'install':
            return _resources_install(args)
        return 1

    # Utilities needed by run modes.
    from pyfuse.utils.common_utils import utils, config

    # setting up output path for run modes
    out_path = os.path.join(args.out_path, 'pyfuse_output_' + datetime.now().strftime("%d-%m-%Y_%I-%M-%S"))
    if not os.path.exists(out_path):
        os.mkdir(out_path)
   

    logger.info("CLI command: " + _normalized_user_command(), extra={"console": False})

    # ------------- setting up logger levels -----------------

    setup_logger(logger, out_path)
    utils.out_path = out_path

    if args.debug:
        logger.info('DEBUG enabled by user')
        logger.setLevel(logging.DEBUG)
    else:
        logger.setLevel(logging.INFO)

    # ----------- running the program based on the input arguments ---------
    if args.command == 'annotator':
        import pandas as pd
        from pyfuse.annotator import RunPyfuse
        from pyfuse.utils.common_utils import utils, config

        resource_root = os.getenv("PYFUSE_RESOURCE_ROOT") or config.get("resource_root")
        res_manager = ResourceManager(resource_root=resource_root)

        display_genome = _display_genome(args.genome)

        if not _arg_was_provided("--genome"):
            logger.info("Genome not provided by user; using default genome: '%s'", display_genome)
        else:
            logger.info("Using user-provided genome: '%s'", display_genome)

        # ----------- loading default settings -----------------
        cfg_assets = res_manager.resolve_config_assets()
        config.update({k: str(v) for k, v in cfg_assets.items()})


        if args.settings:
            user_cfg = utils.load_yaml(Path(args.settings).expanduser())
            config.update(user_cfg)  # shallow merge; can do deep merge if needed


        # loading the required resource files into config dict
        default_bundle = f"default_{args.genome}"
        if args.resource:
            logger.info(f"Resource path provided by user; using: {args.resource}")
            bundle = res_manager.resolve(user_dir=args.resource)
            resource_source = "user-provided"
        else:
            logger.info(
                "Resource path not provided by user; selecting latest installed "
                f"resource bundle for genome '{display_genome}'"
            )
            try:
                bundle = res_manager.resolve_from_cache(default_bundle=default_bundle)
                resource_source = "default-selected"
                logger.info(
                    f"Selected latest resource bundle for genome '{display_genome}': {bundle.root_dir}"
                )
            except FileNotFoundError:
                logger.error(
                    "No installed resource bundle found for "
                    f"{display_genome}. Run: pyfuse resources install --genome {args.genome}."
                )
                return 1

            print(_build_run_manifest(args, str(bundle.root_dir), resource_source))
        files_dict = {k: str(v) for k, v in bundle.resource_files.items()}
        config.update(files_dict)

        # Run the program
        target_df = pd.DataFrame()
        input_bkpt_file = utils.check_path(args.fusion_bkpt, msg=f'input file -'
                                           f'{args.fusion_bkpt} do not exist')
        if args.target_bed:
            target_bed_path = utils.check_path(args.target_bed,
                                               msg=f'input file - {args.target_bed} do not exist')
            target_df = pd.read_csv(str(target_bed_path), sep='\t', header=None)
            target_df = target_df.iloc[:, :3]
            target_df.columns = ['chr', 'start', 'end']

        pyfuse_obj = RunPyfuse(input_bkpt_file, target_df,
                           args.format, ref=args.ref, genome=args.genome, fusion_plot_mode=args.fusion_plot_mode)
        pyfuse_obj.run_pipeline()

    etime = time.time()
    logger.info(f'Time elapsed: {((etime - stime) / 60):.2f} minutes')


# filter messages lower than level (exclusive)
class MaxLevelFilter(logging.Filter):
    def __init__(self, level):
        self.level = level

    def filter(self, record):
        return record.levelno < self.level


class _HelpAction(argparse._HelpAction):
    ''' Formats subparser and its help menu '''

    def __call__(self, parser, namespace, values, option_string=None):
        parser.print_help()

        # retrieve subparsers from parser
        subparsers_actions = [
            action for action in parser._actions
            if isinstance(action, argparse._SubParsersAction)]
        # there will probably only be one subparser_action,
        # but better save than sorry
        print('\n\n\n+++++++++++++ Help for individual modes ++++++++++++\n')
        for subparsers_action in subparsers_actions:
            # get all subparsers and print help
            for choice, subparser in subparsers_action.choices.items():
                print("Mode: '{}'".format(choice), )
                print(subparser.format_help(), '\n\n')

        parser.exit()