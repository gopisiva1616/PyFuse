# Resource Install Options and Examples

Use `pyfuse resources install` to build/install a managed bundle.

Start with help menus:

```bash
pyfuse --help
pyfuse resources --help
pyfuse resources install --help
```

## Operational Notes

- Resource building can take about 30 minutes to 1 hour depending on genome, network speed, and override options.

## Command shape

```bash
pyfuse resources install --genome {grch37|grch38} [options]
```

## Option Summary

- `--genome`: required genome build
- `--source`: `refseq-gtf` (default) or `ucsc-gtf`
- `--gtf`: optional GTF input override (`<url_or_file>`)
- `--annot-summary`: optional assembly summary override (`<url_or_file>`)
- `--mane-file`: optional MANE override (`<url_or_file>`)
- `--cache-dir`: optional cache location for staging/downloads
- `--resource-root`: optional managed resource root

Top-level option shared by all subcommands:

- `-s, --settings`: optional custom settings YAML

## Input files and formats

For installer inputs (`--gtf`, `--annot-summary`, `--mane-file`), accepted forms are:

- URL (downloaded and staged)
- local file path (used directly)

Runtime order:
During runtime, inputs are chosen in the following order until corresponding values are found:

1. explicit CLI option
2. value from custom settings passed with `-s/--settings`
3. packaged defaults from `src/pyfuse/config/settings.yaml`

## Per-Option Behavior

#### `--gtf`
Default behavior: Uses genome-specific RefSeq GTF URL in packaged settings.

  - `grch37`: [NCBI RefSeq GRCh37 genomic GTF](https://ftp.ncbi.nlm.nih.gov/refseq/H_sapiens/annotation/GRCh37_latest/refseq_identifiers/GRCh37_latest_genomic.gtf.gz)

  - `grch38`: [NCBI RefSeq GRCh38 genomic GTF](https://ftp.ncbi.nlm.nih.gov/refseq/H_sapiens/annotation/GRCh38_latest/refseq_identifiers/GRCh38_latest_genomic.gtf.gz)

Override methods:

- CLI: `--gtf <url_or_file>` or 

- Custom settings YAML: set `"default_refseq_gtf_url_grch37"` and `"default_refseq_gtf_url_grch38"` in a file passed with `-s/--settings`


### `--annot-summary`
Default behavior: Uses genome-specific RefSeq assembly summary URL in packaged settings.

- `grch37`: [NCBI RefSeq GRCh37 assembly report](https://ftp.ncbi.nlm.nih.gov/refseq/H_sapiens/annotation/GRCh37_latest/refseq_identifiers/GRCh37_latest_assembly_report.txt)

- `grch38`: [NCBI RefSeq GRCh38 assembly report](https://ftp.ncbi.nlm.nih.gov/refseq/H_sapiens/annotation/GRCh38_latest/refseq_identifiers/GRCh38_latest_assembly_report.txt)

Override methods:
- CLI: `--annot-summary <url_or_file>`
- Custom settings YAML: set `"default_refseq_assembly_summary_url_grch37"` and `"default_refseq_assembly_summary_url_grch38"` in a file passed with `-s/--settings`


### `--mane-file`
Default behavior:
Uses `"default_mane_url"` in packaged settings.

- Default value: [NCBI MANE current directory](https://ftp.ncbi.nlm.nih.gov/refseq/MANE/MANE_human/current/)
- MANE directory URL resolves to latest summary file automatically.
- Direct summary URL ending with `.summary.txt.gz` is also accepted.
- For default GRCh37 setup, MANE is skipped.

Override methods:
- CLI: `--mane-file <url_or_file>`
- Custom settings YAML: set `"default_mane_url"` in a file passed with `-s/--settings`

GRCh37 note:
- If a non-default MANE input is provided, PyFuse warns and asks for explicit confirmation.

Compatibility alias accepted by CLI: `--mane-url`.

Run help to view exact defaults:

```bash
pyfuse resources install -h
```

## Examples

### Default GRCh37 install

```bash
pyfuse resources install --genome grch37
```

### Default GRCh38 install

```bash
pyfuse resources install --genome grch38
```

### List what is installed

```bash
pyfuse resources list
```

### Verify latest installed bundle for a genome

```bash
pyfuse resources verify --genome grch37
pyfuse resources verify --genome grch38
```

### Print selected latest bundle path for a genome

```bash
pyfuse resources path --genome grch37
pyfuse resources path --genome grch38
```

### Override with custom local files

```bash
pyfuse resources install --genome grch38 \
  --gtf /data/ref/GRCh38_latest_genomic.gtf.gz \
  --annot-summary /data/ref/GRCh38_latest_assembly_report.txt
```

### Install and store managed bundles under a custom root

```bash
pyfuse resources install --genome grch38 --resource-root /data/pyfuse_resources
```

### Use a settings file, then override one value from CLI

```bash
pyfuse -s /data/settings.custom.yaml resources install --genome grch38
pyfuse -s /data/settings.custom.yaml resources install --genome grch38 --gtf /data/custom.gtf.gz
```

### Override with custom URLs

```bash
pyfuse resources install --genome grch38 \
  --gtf https://example.org/custom.gtf.gz \
  --annot-summary https://example.org/custom_assembly_report.txt
```

### GRCh37 with explicit back-mapped MANE URL

```bash
pyfuse resources install --genome grch37 --mane-file https://example.org/backmapped_mane/
```

PyFuse will warn and ask for explicit confirmation before using MANE for GRCh37.

## Resource Root Precedence

For resources subcommands (`install`, `list`, `path`, `verify`):

1. `--resource-root`
2. `PYFUSE_RESOURCE_ROOT`
3. `resource_root` in settings YAML
4. platform default user data directory
