# Annotation Column Reference

This page summarizes columns that may appear in PyFuse outputs.

## Core exon annotation columns

- `Fusion_id`: sequential fusion identifier inside a run
- `5'-3'Gene_Partners`: normalized fusion partner gene pair
- `Fusion_Position`: positional relationship of breakpoints across partners
- `Fusion_Annotation`: fusion class label from exon context logic
- `5'_Exon_Annotation`: 5-prime exon-level annotation (`strand|loc|gene|transcript`)
- `3'_Exon_Annotation`: 3-prime exon-level annotation (`strand|loc|gene|transcript`)
- `5'co-ordinate`: genomic coordinate for 5-prime breakpoint
- `3'co-ordinate`: genomic coordinate for 3-prime breakpoint
- `Distance_between_breakpoints`: same-chromosome genomic distance or `NA`
- `Genome`: input genome selected for the run (for example `GRCh37`, `GRCh38`)

## Frame annotation columns

- `Frame_5p`: computed coding frame state of 5-prime partner
- `Frame_3p`: computed coding frame state of 3-prime partner
- `Frame_Status`: fusion-level frame interpretation (for example in-frame/out-of-frame classes)

## MANE annotation columns (when MANE resource is available)

- `5'_MANE_status`: MANE class assignment(s) for 5-prime transcript
- `3'_MANE_status`: MANE class assignment(s) for 3-prime transcript
- `Gene_function`: combined function labels (`5p|3p`) from MANE resource map

## Sequence annotation columns (when `--reference` is provided)

- `Fusion_nucleotide_sequence`: assembled nucleotide sequence around fusion junction
- `Fusion_peptide_sequence`: translated peptide sequence from assembled nucleotide sequence

## Optional external-resource columns (when GTEx/COSMIC resources are enabled)

- `Present_in_COSMIC`
- `Histology`
- `Present_in_GTEX`
- `Average_Expression`
- `Number_of_Tissues_that_contain_fusion`
- `Tissue_Names`

## Optional blacklist columns

When `black_list` resource is configured, PyFuse appends all columns from that blacklist table after merge on `5'-3'Gene_Partners`.
