import pandas as pd

from pyfuse.prep_resources.build_gene_feature_table import build_gene_feature_table


def test_build_gene_feature_table_includes_promoter_and_codons():
    gtf_df = pd.DataFrame(
        {
            "chrom": ["chr1", "chr1", "chr1", "chr1", "chr1"],
            "region": ["transcript", "exon", "exon", "start_codon", "stop_codon"],
            "start": [1000, 1000, 1400, 1010, 1490],
            "end": [1500, 1100, 1500, 1012, 1492],
            "strand": ["+", "+", "+", "+", "+"],
            "gene": ["TMPRSS2", "TMPRSS2", "TMPRSS2", "TMPRSS2", "TMPRSS2"],
            "transcript": ["NM_001", "NM_001", "NM_001", "NM_001", "NM_001"],
            "exon_number": [".", "1", "2", ".", "."],
        }
    )

    out = build_gene_feature_table(gtf_df)

    assert "promoter" in set(out["feature_type"])
    assert "start_codon" in set(out["feature_type"])
    assert "stop_codon" in set(out["feature_type"])
    assert "gene_body" in set(out["feature_type"])
