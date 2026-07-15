# Output Files

Each run creates a timestamped output folder under your `-o` path.

## Files

- `pyfuse_output.xlsx`: full annotation table
- `pyfuse_fusion_annotation.html`: interactive HTML table with sorting/filter/search
- `pyfuse_output.vcf`: BND-style VCF representation
- `pyfuse_fusion_summary.txt`: fusion summary counts
- `excluded_breakpoints.txt`: entries filtered during preprocessing
- `pyfuse_<timestamp>.log`: runtime logs

## HTML report behavior

In the HTML report:

- gene partners are hyperlinked to GeneCards
- coordinates/transcripts are hyperlinked to UCSC and RefSeq locations
- DataTables controls support filtering and search
