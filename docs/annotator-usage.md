# Annotator Usage

`pyfuse annotator` runs the full pipeline from input breakpoints to final report outputs.

For command discovery, use:

```bash
pyfuse --help
pyfuse annotator --help
```

## Basic command

```bash
pyfuse annotator \
  -i /path/to/fusions.tsv \
  -o /path/to/output_dir \
  --input_format star \
  --genome grch37
```

This command uses the latest installed managed bundle for GRCh37 unless `--resource_path` is provided.

## Core options

- `-i, --fusion-breakpoints`: input breakpoint file
- `-o`: output directory
- `--input_format`: caller format parser
- `--genome`: `grch37` or `grch38` (default `grch37`)

## Optional arguments

- `-r, --resource_path`: use exact resource bundle path
- `-g, --reference`: FASTA reference for fusion nucleotide/peptide sequence
- `-t, --target_bed`: BED filter for target regions
- `-s, --settings`: custom settings YAML
- `-d`: debug logging

## Resource resolution behavior

### Default behavior

1. if `--resource_path` is not provided, PyFuse selects the latest managed bundle for `--genome` from `<resource_root>/resources/default_grch*/v*` (Linux default root is typically `~/.local/share/pyfuse`)
2. if no managed bundle exists for that genome, PyFuse exits and asks you to install resources

### Override behavior

1. `--resource_path` pins an exact bundle directory for that run
2. this bypasses managed latest-version auto-selection

Install if missing:

```bash
pyfuse resources install --genome grch37
pyfuse resources install --genome grch38
pyfuse resources install --help
```

## More examples

### Use default bundle by genome selection

```bash
pyfuse annotator \
  -i /data/fusions.star.tsv \
  -o /data/results \
  --input_format star \
  --genome grch38
```

### Use explicit custom resource path (bypasses genome default selection)

```bash
pyfuse annotator \
  -i /data/fusions.arriba.tsv \
  -o /data/results \
  --input_format arriba \
  --resource_path /data/resources/custom_grch38_bundle
```

### Include reference FASTA and target BED

```bash
pyfuse annotator \
  -i /data/fusions.tsv \
  -o /data/results \
  --input_format fusion_inspector \
  --genome grch37 \
  --reference /data/ref/genome.fa \
  --target_bed /data/panel_targets.bed
```

### Use settings file defaults and override one value

```bash
pyfuse -s /data/settings.custom.yaml annotator \
  -i /data/fusions.tsv \
  -o /data/results \
  --input_format star \
  --genome grch38
```

```bash
pyfuse -s /data/settings.custom.yaml annotator \
  -i /data/fusions.tsv \
  -o /data/results \
  --input_format star \
  --genome grch38 \
  --resource_path /data/resources/override_bundle
```
