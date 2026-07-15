# PyFuse

<p align="center">
  <img src="docs/logo/PyFuse_logo.png" alt="PyFuse logo" width="420">
</p>

PyFuse (Python Fusion Annotator) is a command-line tool for annotating gene fusion breakpoints with exon context, frame status, optional sequence context, and report outputs (Excel, HTML, VCF).

## Full Documentation

Use the MkDocs user manual for complete usage details:

- Documentation site: https://gopisiva1616.github.io/PyFuse/
- Docs source index: [docs/index.md](docs/index.md)

Recommended reading order:

1. [Description](docs/description.md)
2. [Installation](docs/installation.md)
3. [Resource Setup Overview](docs/resources-overview.md)
4. [Resource Install Options and Examples](docs/resources-install-options.md)
5. [MANE Usage](docs/mane-transcripts.md)
6. [Annotator Usage](docs/annotator-usage.md)
7. [Output Files](docs/outputs.md)
8. [Annotation Columns](docs/annotation-columns.md)

## Quick Start

```bash
pyfuse resources install --genome grch37
pyfuse annotator -i /path/to/fusions.tsv -o /path/to/output --input_format star --genome grch37
```

## Installation
From Git (latest source)

```bash
git clone https://github.com/gopisiva1616/PyFuse.git
cd PyFuse
python -m venv .venv
source .venv/bin/activate
pip install -U pip
pip install -e .
```

From PyPI:

```bash
pip install pyfuse
```

### Conda 

```bash
conda install -c conda-forge pyfuse
```

## Requirements

- Python 3.10+
- `bedtools` available on PATH (required by `pybedtools` workflows)

## CLI overview

```bash
pyfuse -h
```

Top-level modes:

- `annotator`: run complete fusion annotation
- `resources`: build/install and manage resource bundles

## Annotator mode

```bash
pyfuse annotator \
  -i /path/to/breakpoints.tsv \
  -o /path/to/output_dir \
  --input_format star \
  --genome grch37
```

Common options:

- `-i, --fusion-breakpoints`: input breakpoint file (required)
- `-o`: output directory (required)
- `--input_format`: one of input caller format `default`, `star`, `arriba`, `tophat`, `fusion_catcher`, `longgf`, `fusion_inspector` (required)
- `--genome`: one of `grch37`, `grch38` (default: `grch37`)
- `-r, --resource_path`: exact resource bundle directory to use for this run; default chooses from inbuilt resource budle directory
- `-t, --target_bed`: BED file used for filtering
- `-g, --reference`: FASTA file for fusion sequence annotation
- `-s, --settings`: optional custom settings YAML
- `-d`: enable debug logs

Notes:

- For normal usage, install bundles with `pyfuse resources install` and run `pyfuse annotator` without resource flags.
- Use `--resource_path` only when you want to pin an exact bundle directory for a run.
- Advanced managed-root control is via `PYFUSE_RESOURCE_ROOT` env var or `resource_root` in settings.

## Resource bootstrap and cache

PyFuse can manage local cached resource bundles (one-time setup, reuse later).

```bash
pyfuse resources list
pyfuse resources install --genome grch37
pyfuse resources verify --genome grch37
pyfuse resources path --genome grch37
```

If default source URLs become unavailable, override them directly during install:

```bash
# Override with custom URLs
pyfuse resources install --genome grch38 \
  --gtf https://example.org/custom_grch38.gtf.gz \
  --annot-summary https://example.org/custom_grch38_assembly_report.txt \
  --mane-file https://ftp.ncbi.nlm.nih.gov/refseq/MANE/MANE_human/current/

# Or use local files (no download for these inputs)
pyfuse resources install --genome grch37 \
  --gtf /data/ref/GRCh37_latest_genomic.gtf.gz \
  --annot-summary /data/ref/GRCh37_latest_assembly_report.txt
```

To see all install override options in CLI help:

```bash
pyfuse resources install -h
```

Notes:

- `install` is the primary setup command. It builds and installs a local bundle using defaults from settings or user-provided URLs/files.
- `annotator` uses `--resource_path` if provided; otherwise it uses installed bundles under managed resource storage.
- Packaged default resource fallback is disabled. If a bundle is missing, run `pyfuse resources install --genome <grch37|grch38>`.
- Managed resource root precedence is: `--resource-root` > `PYFUSE_RESOURCE_ROOT` env var > `resource_root` in settings file > platform default data dir.
- Temporary build/download staging uses cache storage (`--cache-dir`), while installed resources are kept in persistent resource storage.

## Output files

For each run, PyFuse writes a timestamped output folder containing:

- `pyfuse_output.xlsx`
- `pyfuse_fusion_annotation.html`
- `pyfuse_output.vcf`
- `pyfuse_fusion_summary.txt`
- `excluded_breakpoints.txt`
- `pyfuse_<timestamp>.log`

## License and attribution

- PyFuse code: GNU GPLv3 (see `LICENSE.md`)
- Bundled third-party web/font assets: see `THIRD_PARTY_NOTICES.md`
- Resource provenance and use restrictions depend on upstream data providers

Attribution request:

- PyFuse is open source under GPLv3 copyleft.
- If PyFuse contributes to your analysis, manuscript, report, presentation, or derivative tool, please credit the project and acknowledge the contributing author/institution .
- Citation is strongly encouraged for academic and scientific use.


## Testing

Run tests with HTML report generation:

```bash
pyfuse-test
```

Integration tests accept runtime environment variables:

- `INPUT_BKPT`, `TRUTH_OUTPUT`, `INPUT_FORMAT`, `REFERENCE_GENOME`, `OUTPUT_PATH`
- `CACHE_DIR`, `RESOURCE_PATH`, `GENOME`, `FUSION2VCF`, `perl`

Example:

```bash
INPUT_FORMAT=star REFERENCE_GENOME=/data/ref/genome.fa pyfuse-test -m integration
```

If `pyfuse-test` is not found, reinstall in your active environment:

```bash
pip install -e .[dev]
```

Fallback invocation:

```bash
python -m pyfuse.test_cli -m integration
```

## Citation

If PyFuse contributes to your analysis or publication, please cite it.
Use machine-readable citation metadata in `CITATION.cff`.
Repository URL: https://github.com/gopisiva1616/PyFuse.git

## Data and Compliance Notes

- Prefer default NCBI RefSeq-based resource generation for a conservative compliance posture.
- Treat optional external datasets (for example COSMIC/GTEx or custom UCSC-derived inputs) as user-managed licensing responsibility.
- Keep provenance records (source URL, date, checksum) for generated resource files.

For practical release guidance and checklists, see `docs/publishing-and-compliance.md`.
For bundled asset licenses/notices, see `THIRD_PARTY_NOTICES.md`.

