# Setup and Usage

This page is a quickstart for first-time users. It summarizes setup and links to authoritative resource references.

## Help menus

Use help menus at each level:

```bash
pyfuse --help
pyfuse annotator --help
pyfuse resources --help
pyfuse resources install --help
```

## Install package

```bash
pip install pyfuse
```

## Resource setup is required

PyFuse annotator requires an installed resource bundle. Setup is genome-specific.

Typical install time:

- Resource building can take around 30 minutes to 1 hour depending on selected genome and options.
- PyFuse reports total elapsed time at the end of both `resources` and `annotator` runs.

You should install at least once for each genome build you plan to run:

```bash
pyfuse resources install --genome grch37
pyfuse resources install --genome grch38
```

If you only work in one genome build, install only that build.

## Resource commands you should know

```bash
pyfuse resources list
pyfuse resources verify --genome grch37
pyfuse resources path --genome grch38
```

- list: shows installed versions per genome
- verify: checks manifest/files for latest version of that genome
- path: prints the exact latest bundle path selected for that genome

## How annotator picks resources

Default behavior:

1. PyFuse selects latest installed bundle for `--genome`.
2. If no bundle exists for that genome, PyFuse exits with install guidance.

Override behavior:

1. If `--resource_path` is provided, that exact bundle is used for that run.

Examples:

```bash
pyfuse annotator -i path/to/fusions.tsv -o output_dir --input_format star --genome grch37
pyfuse annotator -i path/to/fusions.tsv -o output_dir --input_format arriba --resource_path path/to/custom_bundle
```

## Resource root and settings overrides

Resource root precedence for resources commands (install/list/path/verify):

1. --resource-root
2. PYFUSE_RESOURCE_ROOT environment variable
3. resource_root in settings file
4. platform default user data directory

Settings file behavior:

- pass a settings file with -s/--settings
- settings can define `resource_root` and default install inputs (`default_refseq_gtf_url_grch37`, `default_refseq_gtf_url_grch38`, `default_refseq_assembly_summary_url_grch37`, `default_refseq_assembly_summary_url_grch38`, `default_mane_url`)
- explicit CLI flags override settings values for that run

Example:

```bash
pyfuse -s /path/custom_settings.yaml resources install --genome grch38
pyfuse -s /path/custom_settings.yaml resources install --genome grch38 --gtf /data/custom.gtf.gz
```

## Read Next

- [Installation](installation.md)
- [Resource Setup Overview](resources-overview.md)
- [Resource Install Options and Examples](resources-install-options.md)
- [Resource Selection Behavior](resource-selection.md)
- [MANE Usage](mane-transcripts.md)
- [Annotator Usage](annotator-usage.md)
- [Output Files](outputs.md)
- [Annotation Columns](annotation-columns.md)
