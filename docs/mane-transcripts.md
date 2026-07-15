# MANE Transcript Usage

## Default behavior

- For GRCh38 installs, PyFuse enables MANE resource generation by default.
- For GRCh37 installs, MANE is skipped unless a custom MANE input is explicitly provided.

## Why this matters

Official MANE resources are GRCh38-centered. Using MANE with GRCh37 usually implies back-mapped resources and should be treated carefully.


Reference: https://www.ncbi.nlm.nih.gov/refseq/MANE/

## GRCh37 safety confirmation

If you provide `--mane-file` with `--genome grch37`, PyFuse:

1. warns that MANE is GRCh38-centered
2. asks for interactive yes/no confirmation
3. skips MANE generation if you answer no

## Output columns added by MANE

When MANE resource is available, PyFuse adds:

- `5'_MANE_status`
- `3'_MANE_status`
- `Gene_function`

When MANE resource is unavailable, these fields are filled with `.`.
