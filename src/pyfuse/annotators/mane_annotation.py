#!/usr/bin/env python3

"""Annotate fusion breakpoints with MANE transcript status"""

import logging
import pandas as pd
from pyfuse.utils.common_utils import utils, config

logger = logging.getLogger(__name__)


def _normalize_tx_id(tx):
    """Strip version suffix from a transcript accession (e.g. NM_001234.5 -> NM_001234)."""
    tx = str(tx).strip()
    if not tx or tx.lower() in ('nan', '.'):
        return ''
    return tx.split('.')[0]


def _extract_transcript(annotation):
    """Return the transcript ID from an annotation string 'strand|loc|gene|transcript'."""
    if isinstance(annotation, str) and '|' in annotation:
        return annotation.split('|')[-1]
    return ''


def _build_mane_lookup(mane_file):
    """
    Build per-transcript lookup dicts from the mane_status resource file.

    Returns
    -------
    status_map : dict
        {version-stripped transcript ID: semicolon-joined MANE class names}
        '.' when the transcript has no active MANE class.
    function_map : dict
        {version-stripped transcript ID: gene_function string}
    """
    mane_df = pd.read_csv(mane_file, sep='\t', dtype=str, compression='infer')

    fixed_cols = {'gene_symbol', 'gene_function', 'transcript'}
    class_cols = [c for c in mane_df.columns if c not in fixed_cols]

    status_map = {}
    function_map = {}

    for _, row in mane_df.iterrows():
        tx = _normalize_tx_id(row.get('transcript', ''))
        if not tx:
            continue
        active = [c for c in class_cols if str(row.get(c, '0')).strip() == '1']
        status_map[tx] = ';'.join(active) if active else '.'
        gf = row.get('gene_function', '.')
        function_map[tx] = str(gf) if pd.notna(gf) else '.'

    return status_map, function_map


def add_mane(annotated_df):
    """
    Annotate a fusion dataframe with MANE transcript status.

    Adds three new columns:
      - ``5'_MANE_status``  : MANE class(es) for the 5' fusion transcript
                              (e.g. 'MANE Select', 'MANE Plus Clinical', or '.')
      - ``3'_MANE_status``  : MANE class(es) for the 3' fusion transcript
      - ``Gene_function``   : gene-function descriptions for both partners,
                              formatted as '<5p_function>|<3p_function>'

    MANE status resource is loaded from ``config['mane_status']``. When the
    resource is unavailable all three columns are filled with '.'.

    Parameters
    ----------
    annotated_df : pd.DataFrame
        Fusion annotation dataframe containing ``5'_Exon_Annotation`` and
        ``3'_Exon_Annotation`` columns whose values follow the format
        ``strand|loc|gene|transcript``.

    Returns
    -------
    pd.DataFrame
        Input dataframe with the three new columns appended.
    """
    logger.info("-- Adding MANE transcript annotation")

    if 'mane_status' not in config:
        logger.warning(
            "MANE status resource not found in config; "
            "skipping MANE annotation and filling columns with '.'"
        )
        annotated_df = annotated_df.copy()
        annotated_df["5'_MANE_status"] = '.'
        annotated_df["3'_MANE_status"] = '.'
        annotated_df["Gene_function"] = '.'
        return annotated_df

    mane_file = config['mane_status']
    status_map, function_map = _build_mane_lookup(mane_file)

    def lookup_status(annotation):
        tx = _normalize_tx_id(_extract_transcript(annotation))
        return status_map.get(tx, '.') if tx else '.'

    def lookup_function(annotation):
        tx = _normalize_tx_id(_extract_transcript(annotation))
        return function_map.get(tx, '.') if tx else '.'

    annotated_df = annotated_df.copy()
    annotated_df["5'_MANE_status"] = annotated_df["5'_Exon_Annotation"].apply(lookup_status)
    annotated_df["3'_MANE_status"] = annotated_df["3'_Exon_Annotation"].apply(lookup_status)
    annotated_df["Gene_function"] = (
        annotated_df["5'_Exon_Annotation"].apply(lookup_function)
        + '|'
        + annotated_df["3'_Exon_Annotation"].apply(lookup_function)
    )

    logger.info(
        "MANE annotation complete. Matched 5' transcripts: %d / %d; "
        "3' transcripts: %d / %d",
        (annotated_df["5'_MANE_status"] != '.').sum(), len(annotated_df),
        (annotated_df["3'_MANE_status"] != '.').sum(), len(annotated_df),
    )
    return annotated_df

