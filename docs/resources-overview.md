# Resource Setup Overview

PyFuse requires resource bundles to run annotation. Resource setup is genome-specific and should be done at least once per genome build you use.

PyFuse uses publicly available RefSeq GTF inputs to extract gene/transcript mappings, coordinate context, and supporting resource tables required for fusion computation and annotation.

## First-time user help menu

Use these help menus before running commands:

```bash
pyfuse --help
pyfuse resources --help
pyfuse resources install --help
```

## What Setup Produces

Resource setup (`pyfuse resources install`) builds a versioned bundle that contains the files required by `pyfuse annotator`.

At install time, PyFuse:

1. resolves GTF, assembly summary, and optional MANE inputs
2. accepts each input as either URL or local file path
3. downloads URL inputs into staging cache
4. generates normalized resource files
5. installs a versioned bundle under managed storage

## One-Time Setup

Run once per genome you use:

```bash
pyfuse resources install --genome grch37
pyfuse resources install --genome grch38
```

Managed bundle pattern:

- `<resource_root>/resources/default_grch37/v<timestamp>`
- `<resource_root>/resources/default_grch38/v<timestamp>`

After installation, inspect generated bundles with:

```bash
pyfuse resources list
pyfuse resources path --genome grch37
pyfuse resources verify --genome grch37
```

## Input Strategy for `pyfuse resources install` (Summary)

- Typical inputs are RefSeq-style GTF, assembly summary, and optional MANE transcript summary.
- PyFuse ships with default URLs for GRCh37 and GRCh38 in packaged `settings.yaml`.
- Inputs can be overridden by command-line options (`--gtf`, `--annot-summary`, `--mane-file`) or custom settings passed with `-s/--settings`.

Order of input resolution:
- Command line values override custom settings values for that run.
- If not provided by CLI, values can come from keys in a custom settings file passed via `-s/--settings`.
- If not present in custom settings, packaged defaults are used.

Full per-option behavior and exact default URLs are documented in [resources-install-options.md](resources-install-options.md).

## Resource Storage

By default, managed resources are installed under the platform user data root.
On Linux this is typically:

```bash
~/.local/share/pyfuse/resources/
```

Managed root can be overridden in this order for `pyfuse resources` subcommands (`install`, `list`, `path`, `verify`):

1. CLI `--resource-root`
2. environment variable `PYFUSE_RESOURCE_ROOT`
3. settings key `"resource_root"` from a custom settings.yaml file passed to `-s/--settings` argument
4. platform default user data directory

Note: annotator mode does not expose a `--resource-root` option.

## Settings Template

Use [settings.eg.yaml](settings.eg.yaml) as a template for custom settings files.

For complete option-by-option behavior, defaults, and examples, see [resources-install-options.md](resources-install-options.md).