from __future__ import annotations

from pathlib import Path

import pandas as pd


def _normalize_region_token(value: object) -> str:
    token = str(value or "").strip().lower().replace("-", "_").replace(" ", "_")
    if token in {"5utr", "5_utr", "fiveprimeutr", "fiveprime_utr", "utr5"}:
        return "five_prime_utr"
    if token in {"3utr", "3_utr", "threeprimeutr", "threeprime_utr", "utr3"}:
        return "three_prime_utr"
    if token in {"utr", "untranslated_region"}:
        return "utr"
    return token


def _normalize_feature_frame(gtf_df: pd.DataFrame) -> pd.DataFrame:
    df = gtf_df.copy()
    required = {"chrom", "start", "end", "strand", "region", "gene", "transcript"}
    missing = [c for c in required if c not in df.columns]
    if missing:
        missing_txt = ", ".join(missing)
        raise ValueError(f"Gene feature extraction is missing required columns: {missing_txt}")

    df = df.dropna(subset=["chrom", "start", "end", "strand", "region", "gene", "transcript"])
    df["start"] = pd.to_numeric(df["start"], errors="coerce")
    df["end"] = pd.to_numeric(df["end"], errors="coerce")
    df = df.dropna(subset=["start", "end"])
    df["start"] = df["start"].astype(int)
    df["end"] = df["end"].astype(int)
    return df


def _promoter_interval(tx_start: int, tx_end: int, strand: str, upstream_bp: int, downstream_bp: int) -> tuple[int, int]:
    if strand == "-":
        start = max(1, tx_end - downstream_bp)
        end = tx_end + upstream_bp
    else:
        start = max(1, tx_start - upstream_bp)
        end = tx_start + downstream_bp
    return start, end


def _derive_utrs_from_exon_and_cds(exons: pd.DataFrame, cds: pd.DataFrame) -> tuple[pd.DataFrame, pd.DataFrame]:
    """Infer transcript-level 5' and 3' UTR intervals from exon/CDS boundaries."""
    cols = ["chrom", "start", "end", "strand", "gene", "transcript"]
    if exons.empty or cds.empty:
        return pd.DataFrame(columns=cols), pd.DataFrame(columns=cols)

    cds_bounds = (
        cds.groupby(["chrom", "strand", "gene", "transcript"], as_index=False)
        .agg(cds_start=("start", "min"), cds_end=("end", "max"))
    )

    merged = exons[cols].merge(cds_bounds, on=["chrom", "strand", "gene", "transcript"], how="inner")
    if merged.empty:
        return pd.DataFrame(columns=cols), pd.DataFrame(columns=cols)

    utr5_rows: list[dict[str, object]] = []
    utr3_rows: list[dict[str, object]] = []

    for _, row in merged.iterrows():
        exon_start = int(row["start"])
        exon_end = int(row["end"])
        cds_start = int(row["cds_start"])
        cds_end = int(row["cds_end"])
        strand = str(row["strand"])

        left_start = exon_start
        left_end = min(exon_end, cds_start - 1)
        right_start = max(exon_start, cds_end + 1)
        right_end = exon_end

        base = {
            "chrom": row["chrom"],
            "strand": strand,
            "gene": row["gene"],
            "transcript": row["transcript"],
        }

        if left_start <= left_end:
            left = {**base, "start": left_start, "end": left_end}
            if strand == "-":
                utr3_rows.append(left)
            else:
                utr5_rows.append(left)

        if right_start <= right_end:
            right = {**base, "start": right_start, "end": right_end}
            if strand == "-":
                utr5_rows.append(right)
            else:
                utr3_rows.append(right)

    utr5_df = pd.DataFrame(utr5_rows, columns=cols)
    utr3_df = pd.DataFrame(utr3_rows, columns=cols)
    if not utr5_df.empty:
        utr5_df = utr5_df.drop_duplicates().reset_index(drop=True)
    if not utr3_df.empty:
        utr3_df = utr3_df.drop_duplicates().reset_index(drop=True)
    return utr5_df, utr3_df


def build_gene_feature_table(
    gtf_df: pd.DataFrame,
    *,
    upstream_bp: int = 2000,
    downstream_bp: int = 200,
) -> pd.DataFrame:
    """Build a plotting-ready gene feature table from parsed GTF rows.

    Output columns:
      chrom, start, end, strand, gene, transcript,
      feature_type, feature_label, exon_number
    """

    df = _normalize_feature_frame(gtf_df)

    df["region_norm"] = df["region"].map(_normalize_region_token)

    transcripts = df[df["region_norm"] == "transcript"].copy()
    exons = df[df["region_norm"] == "exon"].copy()
    cds = df[df["region_norm"] == "cds"].copy()
    start_codons = df[df["region_norm"] == "start_codon"].copy()
    stop_codons = df[df["region_norm"] == "stop_codon"].copy()
    five_prime_utrs = df[df["region_norm"] == "five_prime_utr"].copy()
    three_prime_utrs = df[df["region_norm"] == "three_prime_utr"].copy()
    generic_utrs = df[df["region_norm"] == "utr"].copy()

    explicit_utr_tx = set(
        pd.concat(
            [
                five_prime_utrs[["chrom", "strand", "gene", "transcript"]],
                three_prime_utrs[["chrom", "strand", "gene", "transcript"]],
                generic_utrs[["chrom", "strand", "gene", "transcript"]],
            ],
            ignore_index=True,
        )
        .drop_duplicates()
        .itertuples(index=False, name=None)
    )

    derived_utr5, derived_utr3 = _derive_utrs_from_exon_and_cds(exons, cds)
    if explicit_utr_tx and (not derived_utr5.empty or not derived_utr3.empty):
        if not derived_utr5.empty:
            keep_mask5 = ~derived_utr5[["chrom", "strand", "gene", "transcript"]].apply(tuple, axis=1).isin(explicit_utr_tx)
            derived_utr5 = derived_utr5[keep_mask5].reset_index(drop=True)
        if not derived_utr3.empty:
            keep_mask3 = ~derived_utr3[["chrom", "strand", "gene", "transcript"]].apply(tuple, axis=1).isin(explicit_utr_tx)
            derived_utr3 = derived_utr3[keep_mask3].reset_index(drop=True)

    if not transcripts.empty:
        tx_core = transcripts[["chrom", "start", "end", "strand", "gene", "transcript"]].copy()
    else:
        # Fallback: infer transcript extents from exon ranges when transcript rows are absent.
        tx_core = (
            exons.groupby(["chrom", "strand", "gene", "transcript"], as_index=False)
            .agg(start=("start", "min"), end=("end", "max"))
        )

    tx_core["feature_type"] = "transcript"
    tx_core["feature_label"] = tx_core["transcript"].astype(str)
    tx_core["exon_number"] = "."

    promoter = tx_core[["chrom", "strand", "gene", "transcript", "start", "end"]].copy()
    promoter[["start", "end"]] = promoter.apply(
        lambda row: _promoter_interval(
            int(row["start"]),
            int(row["end"]),
            str(row["strand"]),
            upstream_bp,
            downstream_bp,
        ),
        axis=1,
        result_type="expand",
    )
    promoter["feature_type"] = "promoter"
    promoter["feature_label"] = "promoter"
    promoter["exon_number"] = "."

    if "exon_number" not in exons.columns:
        exons["exon_number"] = "."
    exon_out = exons[["chrom", "start", "end", "strand", "gene", "transcript", "exon_number"]].copy()
    exon_out["feature_type"] = "exon"
    exon_out["feature_label"] = "E" + exon_out["exon_number"].astype(str)

    sc_out = start_codons[["chrom", "start", "end", "strand", "gene", "transcript"]].copy()
    sc_out["feature_type"] = "start_codon"
    sc_out["feature_label"] = "start_codon"
    sc_out["exon_number"] = "."

    stc_out = stop_codons[["chrom", "start", "end", "strand", "gene", "transcript"]].copy()
    stc_out["feature_type"] = "stop_codon"
    stc_out["feature_label"] = "stop_codon"
    stc_out["exon_number"] = "."

    utr_cols = ["chrom", "start", "end", "strand", "gene", "transcript"]

    utr5_out = five_prime_utrs[utr_cols].copy()
    utr5_out["feature_type"] = "five_prime_utr"
    utr5_out["feature_label"] = "five_prime_utr"
    utr5_out["exon_number"] = "."

    if not derived_utr5.empty:
        derived_utr5_out = derived_utr5[utr_cols].copy()
        derived_utr5_out["feature_type"] = "five_prime_utr"
        derived_utr5_out["feature_label"] = "five_prime_utr"
        derived_utr5_out["exon_number"] = "."
        utr5_out = pd.concat([utr5_out, derived_utr5_out], ignore_index=True)

    utr3_out = three_prime_utrs[utr_cols].copy()
    utr3_out["feature_type"] = "three_prime_utr"
    utr3_out["feature_label"] = "three_prime_utr"
    utr3_out["exon_number"] = "."

    if not derived_utr3.empty:
        derived_utr3_out = derived_utr3[utr_cols].copy()
        derived_utr3_out["feature_type"] = "three_prime_utr"
        derived_utr3_out["feature_label"] = "three_prime_utr"
        derived_utr3_out["exon_number"] = "."
        utr3_out = pd.concat([utr3_out, derived_utr3_out], ignore_index=True)

    utr_out = generic_utrs[utr_cols].copy()
    utr_out["feature_type"] = "utr"
    utr_out["feature_label"] = "utr"
    utr_out["exon_number"] = "."

    gene_body = (
        tx_core.groupby(["chrom", "strand", "gene"], as_index=False)
        .agg(start=("start", "min"), end=("end", "max"))
    )
    gene_body["transcript"] = "."
    gene_body["feature_type"] = "gene_body"
    gene_body["feature_label"] = gene_body["gene"].astype(str)
    gene_body["exon_number"] = "."

    out = pd.concat([gene_body, promoter, tx_core, exon_out, utr5_out, utr3_out, utr_out, sc_out, stc_out], ignore_index=True)
    out = out[["chrom", "start", "end", "strand", "gene", "transcript", "feature_type", "feature_label", "exon_number"]]
    out = out.sort_values(by=["chrom", "gene", "transcript", "start", "end", "feature_type"]).reset_index(drop=True)
    return out


def write_gene_feature_table(gtf_df: pd.DataFrame, out_dir: str | Path, filename: str) -> Path:
    out_dir = Path(out_dir)
    out_dir.mkdir(parents=True, exist_ok=True)
    out_path = out_dir / filename
    table = build_gene_feature_table(gtf_df)
    table.to_csv(out_path, sep="\t", index=False, compression="gzip")
    return out_path
