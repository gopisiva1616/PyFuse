from pathlib import Path

import pandas as pd

from pyfuse.annotators.fusion_plot_links import add_fusion_plot_links
from pyfuse.utils.common_utils import config


def test_add_fusion_plot_links_creates_html_pages(tmp_path):
    feature_df = pd.DataFrame(
        {
            "chrom": ["chr21", "chr21", "chr21", "chr21", "chr21", "chr21", "chr21", "chr21"],
            "start": [1000, 950, 1020, 980, 3000, 2950, 3020, 3070],
            "end": [1100, 1200, 1022, 999, 3100, 3200, 3022, 3090],
            "strand": ["+", "+", "+", "+", "-", "-", "-", "-"],
            "gene": ["TMPRSS2", "TMPRSS2", "TMPRSS2", "TMPRSS2", "ERG", "ERG", "ERG", "ERG"],
            "transcript": ["NM_1", "NM_1", "NM_1", "NM_1", "NM_2", "NM_2", "NM_2", "NM_2"],
            "feature_type": ["exon", "promoter", "start_codon", "five_prime_utr", "exon", "promoter", "start_codon", "three_prime_utr"],
            "feature_label": ["E1", "promoter", "start_codon", "five_prime_utr", "E1", "promoter", "start_codon", "three_prime_utr"],
            "exon_number": ["1", ".", ".", ".", "1", ".", ".", "."],
        }
    )
    feature_path = tmp_path / "gene_feature_table.tsv.gz"
    feature_df.to_csv(feature_path, sep="\t", index=False, compression="gzip")

    original = config.get("gene_feature_table")
    config["gene_feature_table"] = str(feature_path)

    try:
        fusion_df = pd.DataFrame(
            {
                "5'-3'Gene_Partners": ["TMPRSS2-ERG"],
                "5'co-ordinate": ["chr21:1050"],
                "3'co-ordinate": ["chr21:3050"],
            }
        )

        out, embedded = add_fusion_plot_links(fusion_df, tmp_path, mode="external")

        assert "Fusion Visualization" in out.columns
        assert "href=\"fusion_plots/" in out.loc[0, "Fusion Visualization"]
        assert embedded == {}
        plot_dir = tmp_path / "fusion_plots"
        assert plot_dir.is_dir()
        assert any(p.suffix == ".html" for p in plot_dir.iterdir())
    finally:
        if original is None:
            config.pop("gene_feature_table", None)
        else:
            config["gene_feature_table"] = original


def test_add_fusion_plot_links_uses_external_mode_by_default(tmp_path):
    feature_df = pd.DataFrame(
        {
            "chrom": ["chr21", "chr21", "chr21", "chr21", "chr21", "chr21", "chr21", "chr21"],
            "start": [1000, 950, 1020, 980, 3000, 2950, 3020, 3070],
            "end": [1100, 1200, 1022, 999, 3100, 3200, 3022, 3090],
            "strand": ["+", "+", "+", "+", "-", "-", "-", "-"],
            "gene": ["TMPRSS2", "TMPRSS2", "TMPRSS2", "TMPRSS2", "ERG", "ERG", "ERG", "ERG"],
            "transcript": ["NM_1", "NM_1", "NM_1", "NM_1", "NM_2", "NM_2", "NM_2", "NM_2"],
            "feature_type": ["exon", "promoter", "start_codon", "five_prime_utr", "exon", "promoter", "start_codon", "three_prime_utr"],
            "feature_label": ["E1", "promoter", "start_codon", "five_prime_utr", "E1", "promoter", "start_codon", "three_prime_utr"],
            "exon_number": ["1", ".", ".", ".", "1", ".", ".", "."],
        }
    )
    feature_path = tmp_path / "gene_feature_table.tsv.gz"
    feature_df.to_csv(feature_path, sep="\t", index=False, compression="gzip")

    original = config.get("gene_feature_table")
    config["gene_feature_table"] = str(feature_path)

    try:
        fusion_df = pd.DataFrame(
            {
                "5'-3'Gene_Partners": ["TMPRSS2-ERG"],
                "5'co-ordinate": ["chr21:1050"],
                "3'co-ordinate": ["chr21:3050"],
            }
        )

        out, embedded = add_fusion_plot_links(fusion_df, tmp_path)

        assert "Fusion Visualization" in out.columns
        assert "href=\"fusion_plots/" in out.loc[0, "Fusion Visualization"]
        assert embedded == {}
        assert (tmp_path / "fusion_plots").is_dir()
    finally:
        if original is None:
            config.pop("gene_feature_table", None)
        else:
            config["gene_feature_table"] = original
