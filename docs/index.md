# PyFuse User Manual

PyFuse (Python Fusion Annotator) is a one-pass gene fusion annotation workflow with two modes:

- `resources`: one-time, genome-specific bundle setup and management
- `annotator`: end-to-end fusion annotation using managed bundles

## Recommended reading order

1. Description
2. Installation
3. Resource Setup Overview
4. Resource Install Options and Examples
5. MANE Usage
6. Annotator Usage
7. Output Files
8. Annotation Columns

## Quick start

```bash
pyfuse resources install --genome grch37
pyfuse annotator -i /path/to/fusions.tsv -o /path/to/output --input_format star --genome grch37
```

## UI controls

Use the floating `Font size` control (A-, A, A+) in the docs site to adjust text size.
