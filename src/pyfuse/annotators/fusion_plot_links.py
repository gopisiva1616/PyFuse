from __future__ import annotations

from pathlib import Path
import base64
import html
import json
import logging

import pandas as pd

from pyfuse.utils.common_utils import config

logger = logging.getLogger(__name__)


FEATURE_COLORS = {
    "gene_body": "#c7d2fe",
    "transcript": "#bfdbfe",
    "exon": "#38bdf8",
    "five_prime_utr": "#1d4ed8",
    "three_prime_utr": "#1d4ed8",
    "utr": "#1d4ed8",
    "promoter": "#f59e0b",
    "start_codon": "#22c55e",
    "stop_codon": "#ef4444",
}


def _normalize_feature_type(value: object) -> str:
    raw = str(value or "unknown").strip().lower()
    token = raw.replace("'", "").replace("-", "_").replace(" ", "_")
    if token in {"5utr", "5_utr", "fiveprimeutr", "five_prime_utr", "utr5", "fiveprime_utr"}:
        return "five_prime_utr"
    if token in {"3utr", "3_utr", "threeprimeutr", "three_prime_utr", "utr3", "threeprime_utr"}:
        return "three_prime_utr"
    if token in {"utr", "untranslated_region"}:
        return "utr"
    return token


def _insert_as_eighth_column(df: pd.DataFrame, col_name: str, values: object) -> pd.DataFrame:
    out = df.copy()
    if col_name in out.columns:
        out = out.drop(columns=[col_name])
    insert_at = min(7, len(out.columns))
    out.insert(insert_at, col_name, values)
    return out


def _extract_transcript_hint(row: pd.Series, side: str) -> str:
    # Prefer explicit transcript columns when available, then parse exon-annotation format.
    candidates = [
        f"{side}'_Transcript",
        f"{side}'_transcript",
        f"{side}' Transcript",
        f"{side}' transcript",
    ]
    for col in candidates:
        if col in row.index:
            val = str(row.get(col, "")).strip()
            if val and val.lower() != "nan" and val != ".":
                return val

    annot_col = "5'_Exon_Annotation" if side == "5" else "3'_Exon_Annotation"
    annot = str(row.get(annot_col, "")).strip()
    if not annot or annot.lower() in {"nan", "na"}:
        return "ALL"

    parts = [p.strip() for p in annot.split("|") if p.strip()]
    if len(parts) >= 4:
        return parts[-1]
    return "ALL"


def _parse_gene_pair(value: object) -> tuple[str | None, str | None]:
    text = str(value or "").strip()
    if not text:
        return None, None
    for sep in ["-", "~", "|"]:
        if sep in text:
            parts = [p.strip() for p in text.split(sep) if p.strip()]
            if len(parts) >= 2:
                return parts[0], parts[1]
    return None, None


def _parse_coord(value: object) -> tuple[str, int] | None:
    text = str(value or "").strip()
    if not text or text.lower() == "nan":
        return None
    parts = text.split(":")
    if len(parts) < 2:
        return None
    chrom = parts[0]
    try:
        pos = int(parts[1])
    except ValueError:
        return None
    return chrom, pos


def _feature_records(feature_df: pd.DataFrame) -> list[dict[str, object]]:
    rows: list[dict[str, object]] = []
    for _, row in feature_df.iterrows():
        rows.append(
            {
                "chrom": str(row.get("chrom", ".")),
                "start": int(row["start"]),
                "end": int(row["end"]),
                "strand": str(row.get("strand", ".")),
                "gene": str(row.get("gene", ".")),
                "transcript": str(row.get("transcript", ".")),
                "feature_type": _normalize_feature_type(row.get("feature_type", "unknown")),
                "feature_label": str(row.get("feature_label", "unknown")),
                "exon_number": str(row.get("exon_number", ".")),
            }
        )
    return rows


def _transcript_options(feature_df: pd.DataFrame) -> list[str]:
    tx = sorted({str(x) for x in feature_df.get("transcript", pd.Series(dtype=str)).dropna().tolist() if str(x) not in {"", ".", "nan"}})
    return ["ALL"] + tx


def _render_plot_html(
    gene1: str,
    gene2: str,
    bp1: tuple[str, int],
    bp2: tuple[str, int],
    df1: pd.DataFrame,
    df2: pd.DataFrame,
    default_tx_5: str = "ALL",
    default_tx_3: str = "ALL",
) -> str:
    features_1 = _feature_records(df1)
    features_2 = _feature_records(df2)
    transcripts_1 = _transcript_options(df1)
    transcripts_2 = _transcript_options(df2)
    features_1_json = json.dumps(features_1)
    features_2_json = json.dumps(features_2)
    transcripts_1_json = json.dumps(transcripts_1)
    transcripts_2_json = json.dumps(transcripts_2)
    feature_colors_json = json.dumps(FEATURE_COLORS)
    bp1_json = json.dumps({"chrom": bp1[0], "pos": bp1[1]})
    bp2_json = json.dumps({"chrom": bp2[0], "pos": bp2[1]})
    feature_order_json = json.dumps(["promoter", "gene_body", "transcript", "exon", "five_prime_utr", "three_prime_utr", "utr", "start_codon", "stop_codon"])
    default_tx_5_json = json.dumps(default_tx_5 or "ALL")
    default_tx_3_json = json.dumps(default_tx_3 or "ALL")

    return f"""<!DOCTYPE html>
<html>
<head>
    <meta charset=\"utf-8\" />
    <title>{html.escape(gene1)}::{html.escape(gene2)} fusion view</title>
    <style>
        body {{
            margin: 0;
            font-family: 'Segoe UI', Tahoma, sans-serif;
            background: radial-gradient(circle at 0% 0%, #eef6ff, #f8fafc 42%, #e2e8f0 100%);
            color: #0f172a;
            padding: 20px;
        }}
        .card {{
            background: #ffffff;
            border-radius: 14px;
            box-shadow: 0 12px 28px rgba(15, 23, 42, 0.12);
            padding: 16px 18px;
            margin-bottom: 14px;
        }}
        .title {{ font-size: 20px; font-weight: 700; margin-bottom: 8px; }}
        .meta {{ color: #334155; font-size: 13px; margin-bottom: 8px; }}
        .controls {{
            display: grid;
            gap: 10px;
            grid-template-columns: repeat(auto-fit, minmax(260px, 1fr));
            margin-top: 6px;
        }}
        .control-block {{
            background: #f8fafc;
            border: 1px solid #dbeafe;
            border-radius: 10px;
            padding: 10px;
        }}
        .control-block.thin {{
            padding: 6px 8px;
            border-radius: 8px;
        }}
        .control-block label {{ display: block; font-size: 12px; color: #334155; margin-bottom: 4px; }}
        .control-block select, .control-block input[type=range] {{ width: 100%; }}
        .note {{ margin-top: 8px; font-size: 12px; color: #334155; background: #eef6ff; border: 1px solid #dbeafe; border-radius: 8px; padding: 7px 10px; }}
        .legend {{ display: flex; flex-wrap: wrap; gap: 8px; margin-top: 8px; }}
        .chip {{ border-radius: 999px; padding: 4px 9px; font-size: 12px; font-weight: 600; color: #0f172a; display: inline-flex; align-items: center; gap: 6px; }}
        .chip input {{ margin: 0; }}
        .btn {{ border: 1px solid #cbd5e1; background: #f8fafc; color: #0f172a; border-radius: 8px; padding: 6px 10px; font-size: 12px; font-weight: 600; cursor: pointer; }}
        .btn.active {{ background: #dbeafe; color: #1d4ed8; border-color: #93c5fd; }}
        .panel-title {{ font-size: 14px; font-weight: 700; margin: 0; }}
        .panel-head {{ display: flex; align-items: center; justify-content: flex-start; gap: 12px; margin: 2px 0 10px; }}
        .track-wrap {{ border-radius: 12px; border: 1px solid #dbeafe; background: #fbfdff; padding: 8px; margin-top: 10px; overflow-x: hidden; cursor: grab; user-select: none; }}
        .track-wrap.dragging {{ cursor: grabbing; }}
        svg {{ width: 100%; display: block; }}
        .details {{
            margin-top: 10px;
            border-radius: 10px;
            border: 1px dashed #bfdbfe;
            background: #f8fbff;
            padding: 10px;
            font-size: 12px;
            color: #334155;
        }}
        .details strong {{ color: #0f172a; }}
        .hover-box {{
            position: fixed;
            pointer-events: none;
            background: rgba(15, 23, 42, 0.92);
            color: #f8fafc;
            border-radius: 8px;
            padding: 8px 10px;
            font-size: 12px;
            max-width: 420px;
            z-index: 999;
            opacity: 0;
            transition: opacity .08s linear;
        }}
    </style>
</head>
<body>
    <div class=\"card\">
        <div class=\"title\">Fusion structure: {html.escape(gene1)} :: {html.escape(gene2)}</div>
        <div class=\"meta\">5-prime breakpoint: {html.escape(bp1[0])}:{bp1[1]} | 3-prime breakpoint: {html.escape(bp2[0])}:{bp2[1]}</div>
        <div class=\"controls\">
            <div class=\"control-block\">
                <label for=\"tx5\">Gene A transcript view</label>
                <select id=\"tx5\"></select>
            </div>
            <div class=\"control-block\">
                <label for=\"tx3\">Gene B transcript view</label>
                <select id=\"tx3\"></select>
            </div>
            <div class=\"control-block thin\">
                <label for="zoomPairPct">Block1 Zoom (%): <span id="zoomPairPctValue">0</span></label>
                <input id="zoomPairPct" type="range" min="0" max="1900" step="10" value="0" />
            </div>
            <div class=\"control-block thin\">
                <label for="zoomFusedPct">Block2 Zoom (%): <span id="zoomFusedPctValue">0</span></label>
                <input id="zoomFusedPct" type="range" min="0" max="1900" step="10" value="0" />
            </div>
        </div>
        <div class=\"legend\">
            <label class=\"chip\" style=\"background:{FEATURE_COLORS['gene_body']}\"><input data-ftype=\"gene_body\" type=\"checkbox\"/>gene body</label>
            <label class=\"chip\" style=\"background:{FEATURE_COLORS['transcript']}\"><input data-ftype=\"transcript\" type=\"checkbox\"/>transcript</label>
            <label class=\"chip\" style=\"background:{FEATURE_COLORS['exon']}\"><input data-ftype=\"exon\" type=\"checkbox\" checked/>exon</label>
            <label class=\"chip\" style=\"background:{FEATURE_COLORS['five_prime_utr']}\"><input data-ftype=\"five_prime_utr\" type=\"checkbox\" checked/>5' UTR</label>
            <label class=\"chip\" style=\"background:{FEATURE_COLORS['three_prime_utr']}\"><input data-ftype=\"three_prime_utr\" type=\"checkbox\" checked/>3' UTR</label>
            <label class=\"chip\" style=\"background:{FEATURE_COLORS['utr']}\"><input data-ftype=\"utr\" type=\"checkbox\" checked/>UTR</label>
            <label class=\"chip\" style=\"background:{FEATURE_COLORS['promoter']}\"><input data-ftype=\"promoter\" type=\"checkbox\" checked/>promoter</label>
            <label class=\"chip\" style=\"background:{FEATURE_COLORS['start_codon']}\"><input data-ftype=\"start_codon\" type=\"checkbox\" checked/>start codon</label>
            <label class=\"chip\" style=\"background:{FEATURE_COLORS['stop_codon']}\"><input data-ftype=\"stop_codon\" type=\"checkbox\" checked/>stop codon</label>
        </div>
        <div class=\"details\" id=\"details\"><strong>Hover a feature</strong> to inspect type, transcript, and coordinates.</div>
        <div class=\"note\">Use mouse wheel over each block to zoom that block only. Click and drag in a block to pan left/right. Exon labels are shown below the track; promoter/UTR/start/stop labels are shown above.</div>
    </div>

    <div class=\"card\">
        <div class=\"panel-head\">
            <div class=\"panel-title\">Block 1: Gene A and Gene B side-by-side</div>
            <button id=\"toggleLabelsBtn\" class=\"btn active\" type=\"button\">Annotation labels: ON</button>
        </div>
        <div id=\"pairWrap\" class=\"track-wrap\"><svg id=\"pairView\" height=\"250\"></svg></div>
    </div>

    <div class=\"card\">
        <div class=\"panel-title\">Block 2: Fused product (Gene A 5-prime segment + Gene B 3-prime segment)</div>
        <div id=\"fusedWrap\" class=\"track-wrap\"><svg id=\"fusedView\" height=\"170\"></svg></div>
    </div>

    <div class=\"hover-box\" id=\"hoverBox\"></div>

    <script>
        const features5 = {features_1_json};
        const features3 = {features_2_json};
        const transcriptOptions5 = {transcripts_1_json};
        const transcriptOptions3 = {transcripts_2_json};
        const featureColors = {feature_colors_json};
        const featureOrder = {feature_order_json};
        const bp5 = {bp1_json};
        const bp3 = {bp2_json};
        const defaultTx5 = {default_tx_5_json};
        const defaultTx3 = {default_tx_3_json};

        const svgNs = "http://www.w3.org/2000/svg";
        const canvasW = 1700;
        const pairH = 250;
        const fusedH = 170;
        const pairLeftX0 = 30;
        const pairLeftX1 = 810;
        const pairRightX0 = 890;
        const pairRightX1 = 1670;
        const pairSplitX = 850;
        const fusedX0 = 50;
        const fusedJunctionX = 850;
        const fusedX1 = 1650;

        const defaultEnabledTypes = new Set(["exon", "five_prime_utr", "three_prime_utr", "utr", "promoter", "start_codon", "stop_codon"]);
        const enabledTypes = new Set(defaultEnabledTypes);
        const panStatePair = {{ g5: 0, g3: 0 }};
        const panStateFused = {{ g5: 0, g3: 0 }};
        let labelsOn = true;
        let dragState = null;

        const tx5El = document.getElementById("tx5");
        const tx3El = document.getElementById("tx3");
        const zoomPairEl = document.getElementById("zoomPairPct");
        const zoomPairValueEl = document.getElementById("zoomPairPctValue");
        const zoomFusedEl = document.getElementById("zoomFusedPct");
        const zoomFusedValueEl = document.getElementById("zoomFusedPctValue");
        const toggleLabelsBtn = document.getElementById("toggleLabelsBtn");
        const pairWrapEl = document.getElementById("pairWrap");
        const fusedWrapEl = document.getElementById("fusedWrap");
        const hoverBox = document.getElementById("hoverBox");
        const detailsEl = document.getElementById("details");

        function fillSelect(el, options) {{
            options.forEach((opt) => {{
                const node = document.createElement("option");
                node.value = opt;
                node.textContent = opt;
                el.appendChild(node);
            }});
        }}

        function setPreferredTranscript(el, preferred) {{
            if (!preferred || preferred === "ALL") {{
                el.value = "ALL";
                return;
            }}
            const opts = Array.from(el.options).map((o) => o.value);
            if (opts.includes(preferred)) {{
                el.value = preferred;
                return;
            }}
            const pv = preferred.split(".")[0];
            const versionless = opts.find((o) => o.split(".")[0] === pv);
            el.value = versionless || "ALL";
        }}

        function mkSvg(tag, attrs) {{
            const node = document.createElementNS(svgNs, tag);
            Object.entries(attrs).forEach(([k, v]) => node.setAttribute(k, String(v)));
            return node;
        }}

        function clamp(value, minV, maxV) {{
            return Math.min(maxV, Math.max(minV, value));
        }}

        function xScale(value, minV, maxV, x0, x1) {{
            if (maxV <= minV) return x0;
            return ((value - minV) / (maxV - minV)) * (x1 - x0) + x0;
        }}

        function xScaleDirectional(value, minV, maxV, x0, x1, reverse = false) {{
            if (maxV <= minV) return x0;
            const t = (value - minV) / (maxV - minV);
            const u = reverse ? (1 - t) : t;
            return u * (x1 - x0) + x0;
        }}

        function filteredFeatures(features, selectedTx) {{
            return features.filter((f) => {{
                if (!enabledTypes.has(f.feature_type)) return false;
                if (selectedTx === "ALL") return true;
                return f.transcript === selectedTx || f.transcript === ".";
            }});
        }}

        function domainAroundBreakpoint(features, breakpoint, selectedTx, zoomFactor, panFraction) {{
            const byTx = features.filter((f) => selectedTx === "ALL" || f.transcript === selectedTx || f.transcript === ".");
            let baseMin = breakpoint.pos;
            let baseMax = breakpoint.pos;
            byTx.forEach((f) => {{
                if (f.start < baseMin) baseMin = f.start;
                if (f.end > baseMax) baseMax = f.end;
            }});
            const span = Math.max(1, baseMax - baseMin);
            const pad = Math.max(200, Math.round(span * 0.08));
            baseMin -= pad;
            baseMax += pad;

            const totalSpan = Math.max(1, baseMax - baseMin);
            const viewSpan = Math.max(120, totalSpan / Math.max(1, zoomFactor));
            const halfWindow = viewSpan / 2;
            const baseCenter = (baseMin + baseMax) / 2;
            const maxPanShift = Math.max(0, (totalSpan - viewSpan) / 2);
            const center = baseCenter + clamp(panFraction, -1, 1) * maxPanShift;

            let minV = Math.floor(center - halfWindow);
            let maxV = Math.ceil(center + halfWindow);
            if (minV < baseMin) {{
                maxV += (baseMin - minV);
                minV = baseMin;
            }}
            if (maxV > baseMax) {{
                minV -= (maxV - baseMax);
                maxV = baseMax;
            }}
            if (maxV <= minV) maxV = minV + 1;
            return [minV, maxV];
        }}

        function retainedBaseRange(features, breakpoint, selectedTx, keepLower) {{
            const byTx = features.filter((f) => selectedTx === "ALL" || f.transcript === selectedTx || f.transcript === ".");
            let minR = Number.POSITIVE_INFINITY;
            let maxR = Number.NEGATIVE_INFINITY;

            byTx.forEach((f) => {{
                const rs = keepLower ? f.start : Math.max(f.start, breakpoint.pos);
                const re = keepLower ? Math.min(f.end, breakpoint.pos) : f.end;
                if (re > rs) {{
                    if (rs < minR) minR = rs;
                    if (re > maxR) maxR = re;
                }}
            }});

            if (!Number.isFinite(minR) || !Number.isFinite(maxR) || maxR <= minR) {{
                minR = breakpoint.pos - 100;
                maxR = breakpoint.pos + 100;
            }}
            return [minR, maxR];
        }}

        function domainFromBaseRange(baseMin, baseMax, zoomFactor, panFraction) {{
            const span = Math.max(1, baseMax - baseMin);
            const pad = Math.max(30, Math.round(span * 0.04));
            const paddedMin = baseMin - pad;
            const paddedMax = baseMax + pad;
            const totalSpan = Math.max(1, paddedMax - paddedMin);
            const viewSpan = Math.max(60, totalSpan / Math.max(1, zoomFactor));
            const halfWindow = viewSpan / 2;
            const baseCenter = (paddedMin + paddedMax) / 2;
            const maxPanShift = Math.max(0, (totalSpan - viewSpan) / 2);
            const center = baseCenter + clamp(panFraction, -1, 1) * maxPanShift;

            let minV = Math.floor(center - halfWindow);
            let maxV = Math.ceil(center + halfWindow);
            if (minV < paddedMin) {{
                maxV += (paddedMin - minV);
                minV = paddedMin;
            }}
            if (maxV > paddedMax) {{
                minV -= (maxV - paddedMax);
                maxV = paddedMax;
            }}
            if (maxV <= minV) maxV = minV + 1;
            return [minV, maxV];
        }}

        function retainedWindowDomain(features, breakpoint, selectedTx, zoomFactor, panFraction, keepLower) {{
            const [baseMin, baseMax] = retainedBaseRange(features, breakpoint, selectedTx, keepLower);
            return domainFromBaseRange(baseMin, baseMax, zoomFactor, panFraction);
        }}

        function panForAnchor(baseMin, baseMax, zoomFactor, ratio, anchorValue) {{
            const span = Math.max(1, baseMax - baseMin);
            const pad = Math.max(30, Math.round(span * 0.04));
            const paddedMin = baseMin - pad;
            const paddedMax = baseMax + pad;
            const totalSpan = Math.max(1, paddedMax - paddedMin);
            const viewSpan = Math.max(60, totalSpan / Math.max(1, zoomFactor));
            const baseCenter = (paddedMin + paddedMax) / 2;
            const maxPanShift = Math.max(0, (totalSpan - viewSpan) / 2);
            if (maxPanShift <= 0) return 0;

            const center = anchorValue - (clamp(ratio, 0, 1) - 0.5) * viewSpan;
            return clamp((center - baseCenter) / maxPanShift, -1, 1);
        }}

        function adjustedPanForRetainedZoom(features, breakpoint, selectedTx, oldZoomPct, newZoomPct, oldPan, ratio, keepLower) {{
            const oldZoom = 1 + Math.max(0, oldZoomPct) / 100;
            const newZoom = 1 + Math.max(0, newZoomPct) / 100;
            const [baseMin, baseMax] = retainedBaseRange(features, breakpoint, selectedTx, keepLower);
            const oldDomain = domainFromBaseRange(baseMin, baseMax, oldZoom, oldPan);
            const anchor = oldDomain[0] + clamp(ratio, 0, 1) * (oldDomain[1] - oldDomain[0]);
            return panForAnchor(baseMin, baseMax, newZoom, ratio, anchor);
        }}

        function showHover(evt, feature) {{
            hoverBox.style.opacity = "1";
            hoverBox.style.left = `${{evt.clientX + 12}}px`;
            hoverBox.style.top = `${{evt.clientY + 12}}px`;
            const exonTxt = (feature.exon_number && feature.exon_number !== "." && feature.exon_number !== "nan") ? ` | exon: E${{feature.exon_number}}` : "";
            hoverBox.innerHTML = `<strong>${{feature.feature_type}}</strong><br>${{feature.gene}} | ${{feature.transcript}}${{exonTxt}}<br>${{feature.chrom}}:${{feature.start}}-${{feature.end}}`;
            detailsEl.innerHTML = `<strong>Selected feature</strong><br>Type: ${{feature.feature_type}}<br>Label: ${{feature.feature_label}}<br>Gene: ${{feature.gene}}<br>Transcript: ${{feature.transcript}}<br>Exon: ${{(feature.exon_number && feature.exon_number !== "." && feature.exon_number !== "nan") ? `E${{feature.exon_number}}` : "NA"}}<br>Coordinates: ${{feature.chrom}}:${{feature.start}}-${{feature.end}}`;
        }}

        function hideHover() {{
            hoverBox.style.opacity = "0";
        }}

        function inferDirection(features) {{
            let plus = 0;
            let minus = 0;
            features.forEach((f) => {{
                if (f.strand === "+") plus += 1;
                else if (f.strand === "-") minus += 1;
            }});
            return minus > plus ? -1 : 1;
        }}

        function featureHeight(featureType, context) {{
            const baseExon = context === "fused" ? 22 : 20;
            if (featureType === "exon") return baseExon;
            if (featureType === "promoter" || featureType === "five_prime_utr" || featureType === "three_prime_utr" || featureType === "utr") return baseExon - 6;
            return baseExon - 4;
        }}

        function featureTag(feature) {{
            if (feature.feature_type === "exon") {{
                if (feature.exon_number && feature.exon_number !== "." && feature.exon_number !== "nan") return `Ex${{feature.exon_number}}`;
                return "Ex";
            }}
            if (feature.feature_type === "promoter") return "Promotor";
            if (feature.feature_type === "five_prime_utr") return "5UTR";
            if (feature.feature_type === "three_prime_utr") return "3UTR";
            if (feature.feature_type === "utr") return "UTR";
            if (feature.feature_type === "start_codon") return "Start";
            if (feature.feature_type === "stop_codon") return "Stop";
            if (feature.feature_type === "transcript") return "Tx";
            if (feature.feature_type === "gene_body") return "Gene";
            return feature.feature_type;
        }}

        function drawLabelCallout(svg, fx0, fx1, y, h, label, isExon) {{
            const cx = (fx0 + fx1) / 2;
            const anchorY = isExon ? y + h / 2 : y - h / 2;
            const textY = isExon ? anchorY + 10 : anchorY - 8;

            const text = mkSvg("text", {{
                x: cx,
                y: textY,
                "font-size": 8,
                fill: "#1f2937",
                "font-weight": 600,
                "text-anchor": "middle",
            }});
            text.textContent = label;
            svg.appendChild(text);

            const leader = mkSvg("line", {{
                x1: cx,
                y1: isExon ? textY - 6 : textY + 3,
                x2: cx,
                y2: anchorY,
                stroke: "#64748b",
                "stroke-width": 0.8,
            }});
            svg.appendChild(leader);

            const dir = isExon ? -1 : 1;
            const head = mkSvg("path", {{
                d: `M ${{cx - 2}} ${{anchorY - 2 * dir}} L ${{cx}} ${{anchorY}} L ${{cx + 2}} ${{anchorY - 2 * dir}}`,
                stroke: "#64748b",
                "stroke-width": 0.8,
                fill: "none",
                "stroke-linecap": "round",
                "stroke-linejoin": "round",
            }});
            svg.appendChild(head);
        }}

        function drawTerminalDirectionArrows(svg, x0, x1, y, direction) {{
            function drawArrow(cx, dir) {{
                const p = mkSvg("path", {{
                    d: `M ${{cx - 4 * dir}} ${{y - 3}} L ${{cx}} ${{y}} L ${{cx - 4 * dir}} ${{y + 3}}`,
                    stroke: "#334155",
                    "stroke-width": 1,
                    fill: "none",
                    "stroke-linecap": "round",
                    "stroke-linejoin": "round",
                }});
                svg.appendChild(p);
            }}

            if (direction > 0) {{
                drawArrow(x0 + 10, 1);
                drawArrow(x1 - 10, 1);
            }} else {{
                drawArrow(x0 + 10, -1);
                drawArrow(x1 - 10, -1);
            }}
        }}

        function drawDirectionArrow(svg, cx, y, direction, scale = 1) {{
            const p = mkSvg("path", {{
                d: `M ${{cx - 4 * direction * scale}} ${{y - 3 * scale}} L ${{cx}} ${{y}} L ${{cx - 4 * direction * scale}} ${{y + 3 * scale}}`,
                stroke: "#334155",
                "stroke-width": 1,
                fill: "none",
                "stroke-linecap": "round",
                "stroke-linejoin": "round",
            }});
            svg.appendChild(p);
        }}

        function mergedIntervals(intervals) {{
            if (!intervals.length) return [];
            const sorted = intervals.slice().sort((a, b) => a[0] - b[0]);
            const merged = [sorted[0].slice()];
            for (let i = 1; i < sorted.length; i += 1) {{
                const cur = sorted[i];
                const last = merged[merged.length - 1];
                if (cur[0] <= last[1] + 1) {{
                    last[1] = Math.max(last[1], cur[1]);
                }} else {{
                    merged.push(cur.slice());
                }}
            }}
            return merged;
        }}

        function drawIntronDirectionArrows(svg, features, minV, maxV, x0, x1, y, direction, reverse = false) {{
            const exons = mergedIntervals(
                features
                    .filter((f) => f.feature_type === "exon")
                    .map((f) => [Math.max(minV, f.start), Math.min(maxV, f.end)])
                    .filter((iv) => iv[1] > iv[0])
            );

            const introns = [];
            let cursor = minV;
            exons.forEach((iv) => {{
                if (iv[0] > cursor) introns.push([cursor, iv[0]]);
                cursor = Math.max(cursor, iv[1]);
            }});
            if (cursor < maxV) introns.push([cursor, maxV]);

            introns.forEach((iv) => {{
                const pxA = xScaleDirectional(iv[0], minV, maxV, x0, x1, reverse);
                const pxB = xScaleDirectional(iv[1], minV, maxV, x0, x1, reverse);
                const ix0 = Math.min(pxA, pxB);
                const ix1 = Math.max(pxA, pxB);
                if (ix1 - ix0 < 64) return;

                const spacing = 150;
                for (let cx = ix0 + 28; cx <= ix1 - 28; cx += spacing) {{
                    drawDirectionArrow(svg, cx, y, direction, 1);
                }}
            }});
        }}

        function drawSingleTrack(svg, features, domain, x0, x1, y, label, breakpoint, colorBaseline) {{
            const [minV, maxV] = domain;
            const direction = inferDirection(features);
            const baseline = mkSvg("line", {{ x1: x0, y1: y, x2: x1, y2: y, stroke: colorBaseline, "stroke-width": 7, "stroke-linecap": "round" }});
            svg.appendChild(baseline);
            drawIntronDirectionArrows(svg, features, minV, maxV, x0, x1, y, direction);
            drawTerminalDirectionArrows(svg, x0, x1, y, direction);

            const geneLabel = mkSvg("text", {{ x: (x0 + x1) / 2, y: y + 96, "font-size": 12, fill: "#0f172a", "font-weight": 700, "text-anchor": "middle" }});
            geneLabel.textContent = label;
            svg.appendChild(geneLabel);

            const fiveX = direction > 0 ? x0 : x1 - 16;
            const threeX = direction > 0 ? x1 - 16 : x0;
            const five = mkSvg("text", {{ x: fiveX, y: y + 30, "font-size": 10, fill: "#334155", "font-weight": 700 }});
            five.textContent = "5'";
            svg.appendChild(five);
            const three = mkSvg("text", {{ x: threeX, y: y + 30, "font-size": 10, fill: "#334155", "font-weight": 700 }});
            three.textContent = "3'";
            svg.appendChild(three);

            const orderedFeatures = [
                ...features.filter((f) => f.feature_type === "exon"),
                ...features.filter((f) => f.feature_type !== "exon"),
            ];

            orderedFeatures.forEach((f) => {{
                const s = Math.max(f.start, minV);
                const e = Math.min(f.end, maxV);
                if (e <= s) return;
                const fx0 = xScale(s, minV, maxV, x0, x1);
                const fx1 = xScale(e, minV, maxV, x0, x1);
                const w = Math.max(2, fx1 - fx0);
                const h = featureHeight(f.feature_type, "pair");
                const rect = mkSvg("rect", {{
                    x: fx0,
                    y: y - h / 2,
                    width: w,
                    height: h,
                    rx: 3,
                    fill: featureColors[f.feature_type] || "#94a3b8",
                    "fill-opacity": 0.95,
                    stroke: "#0f172a",
                    "stroke-opacity": 0.2,
                    "stroke-width": 0.8,
                }});
                rect.addEventListener("mousemove", (evt) => showHover(evt, f));
                rect.addEventListener("mouseleave", hideHover);
                svg.appendChild(rect);

                if (labelsOn && w >= 8) {{
                    drawLabelCallout(svg, fx0, fx1, y, h, featureTag(f), f.feature_type === "exon");
                }}
            }});

            const bpPos = Math.max(minV, Math.min(maxV, breakpoint.pos));
            const bpX = xScale(bpPos, minV, maxV, x0, x1);
            const bpLine = mkSvg("line", {{ x1: bpX, y1: y - 30, x2: bpX, y2: y + 30, stroke: "#0ea5e9", "stroke-width": 2.2, "stroke-dasharray": "5 4" }});
            svg.appendChild(bpLine);
            const bpTxt = mkSvg("text", {{ x: bpX + 4, y: y - 30, "font-size": 10, fill: "#0ea5e9", "font-weight": 700 }});
            bpTxt.textContent = `bp:${{breakpoint.chrom}}:${{breakpoint.pos}}`;
            svg.appendChild(bpTxt);
            return bpX;
        }}

        function renderPairBlock(f5, f3, tx5, tx3, zoomFactor) {{
            const svg = document.getElementById("pairView");
            while (svg.firstChild) svg.removeChild(svg.firstChild);
            svg.setAttribute("viewBox", `0 0 ${{canvasW}} ${{pairH}}`);

            const y = 110;
            const d5 = domainAroundBreakpoint(features5, bp5, tx5, zoomFactor, panStatePair.g5);
            const d3 = domainAroundBreakpoint(features3, bp3, tx3, zoomFactor, panStatePair.g3);

            svg.appendChild(mkSvg("rect", {{ x: 0, y: 0, width: canvasW, height: pairH, rx: 12, fill: "#f8fbff" }}));

            const windowText = mkSvg("text", {{ x: 20, y: 24, "font-size": 12, fill: "#334155" }});
            windowText.textContent = `Window A: ${{d5[0]}}-${{d5[1]}} | Window B: ${{d3[0]}}-${{d3[1]}}`;
            svg.appendChild(windowText);

            const bpX5 = drawSingleTrack(svg, f5, d5, pairLeftX0, pairLeftX1, y, `{html.escape(gene1)} ({html.escape(bp1[0])})`, bp5, "#7c93d6");
            const bpX3 = drawSingleTrack(svg, f3, d3, pairRightX0, pairRightX1, y, `{html.escape(gene2)} ({html.escape(bp2[0])})`, bp3, "#9ca3af");

            const link = mkSvg("path", {{
                d: `M ${{bpX5}} ${{y + 34}} C ${{bpX5 + 80}} ${{y + 74}}, ${{bpX3 - 80}} ${{y + 74}}, ${{bpX3}} ${{y + 34}}`,
                fill: "none",
                stroke: "#2563eb",
                "stroke-width": 2.6,
                "stroke-linecap": "round",
            }});
            svg.appendChild(link);
        }}

        function renderFusedBlock(f5, f3, tx5, tx3, zoomFactor) {{
            const svg = document.getElementById("fusedView");
            while (svg.firstChild) svg.removeChild(svg.firstChild);
            svg.setAttribute("viewBox", `0 0 ${{canvasW}} ${{fusedH}}`);

            const y = 94;
            const dir5 = inferDirection(f5);
            const dir3 = inferDirection(f3);
            const keepLower5 = dir5 > 0;
            const keepLower3 = dir3 < 0;
            const d5 = retainedWindowDomain(features5, bp5, tx5, zoomFactor, panStateFused.g5, keepLower5);
            const d3 = retainedWindowDomain(features3, bp3, tx3, zoomFactor, panStateFused.g3, keepLower3);
            let leftDomainMin = d5[0];
            let leftDomainMax = d5[1];
            let rightDomainMin = d3[0];
            let rightDomainMax = d3[1];
            const reverseLeft = dir5 < 0;
            const reverseRight = dir3 < 0;
            if (leftDomainMax <= leftDomainMin) leftDomainMax = leftDomainMin + 1;
            if (rightDomainMax <= rightDomainMin) rightDomainMax = rightDomainMin + 1;

            svg.appendChild(mkSvg("rect", {{ x: 0, y: 0, width: canvasW, height: fusedH, rx: 12, fill: "#f8fbff" }}));

            const desc = mkSvg("text", {{ x: 20, y: 22, "font-size": 12, fill: "#334155" }});
            desc.textContent = "Fused product track (left: Gene A retained to breakpoint, right: Gene B retained from breakpoint)";
            svg.appendChild(desc);

            svg.appendChild(mkSvg("line", {{ x1: fusedX0, y1: y, x2: fusedJunctionX, y2: y, stroke: "#5b7fd1", "stroke-width": 8, "stroke-linecap": "round" }}));
            svg.appendChild(mkSvg("line", {{ x1: fusedJunctionX, y1: y, x2: fusedX1, y2: y, stroke: "#64748b", "stroke-width": 8, "stroke-linecap": "round" }}));
            // Fused transcript is shown in a single 5'->3' direction left-to-right.
            drawIntronDirectionArrows(svg, f5, leftDomainMin, leftDomainMax, fusedX0, fusedJunctionX, y, 1, reverseLeft);
            drawIntronDirectionArrows(svg, f3, rightDomainMin, rightDomainMax, fusedJunctionX, fusedX1, y, 1, reverseRight);
            drawTerminalDirectionArrows(svg, fusedX0, fusedJunctionX, y, 1);
            drawTerminalDirectionArrows(svg, fusedJunctionX, fusedX1, y, 1);

            const five = mkSvg("text", {{ x: fusedX0, y: y + 34, "font-size": 10, fill: "#334155", "font-weight": 700 }});
            five.textContent = "5' retained";
            svg.appendChild(five);
            const three = mkSvg("text", {{ x: fusedX1 - 64, y: y + 34, "font-size": 10, fill: "#334155", "font-weight": 700 }});
            three.textContent = "3' retained";
            svg.appendChild(three);

            const orderedF5 = [
                ...f5.filter((f) => f.feature_type === "exon"),
                ...f5.filter((f) => f.feature_type !== "exon"),
            ];
            const orderedF3 = [
                ...f3.filter((f) => f.feature_type === "exon"),
                ...f3.filter((f) => f.feature_type !== "exon"),
            ];

            orderedF5.forEach((f) => {{
                if (f.end <= leftDomainMin || f.start >= leftDomainMax) return;
                const s = Math.max(f.start, leftDomainMin);
                const e = Math.min(f.end, leftDomainMax);
                if (e <= s) return;
                const px0 = xScaleDirectional(s, leftDomainMin, leftDomainMax, fusedX0, fusedJunctionX, reverseLeft);
                const px1 = xScaleDirectional(e, leftDomainMin, leftDomainMax, fusedX0, fusedJunctionX, reverseLeft);
                const fx0 = Math.min(px0, px1);
                const fx1 = Math.max(px0, px1);
                const w = Math.max(2, fx1 - fx0);
                const h = featureHeight(f.feature_type, "fused");
                const rect = mkSvg("rect", {{ x: fx0, y: y - h / 2, width: w, height: h, rx: 3, fill: featureColors[f.feature_type] || "#94a3b8", stroke: "#0f172a", "stroke-opacity": 0.2, "stroke-width": 0.8 }});
                rect.addEventListener("mousemove", (evt) => showHover(evt, f));
                rect.addEventListener("mouseleave", hideHover);
                svg.appendChild(rect);
                if (labelsOn && w >= 8) drawLabelCallout(svg, fx0, fx1, y, h, featureTag(f), f.feature_type === "exon");
            }});

            orderedF3.forEach((f) => {{
                if (f.start >= rightDomainMax || f.end <= rightDomainMin) return;
                const s = Math.max(f.start, rightDomainMin);
                const e = Math.min(f.end, rightDomainMax);
                if (e <= s) return;
                const px0 = xScaleDirectional(s, rightDomainMin, rightDomainMax, fusedJunctionX, fusedX1, reverseRight);
                const px1 = xScaleDirectional(e, rightDomainMin, rightDomainMax, fusedJunctionX, fusedX1, reverseRight);
                const fx0 = Math.min(px0, px1);
                const fx1 = Math.max(px0, px1);
                const w = Math.max(2, fx1 - fx0);
                const h = featureHeight(f.feature_type, "fused");
                const rect = mkSvg("rect", {{ x: fx0, y: y - h / 2, width: w, height: h, rx: 3, fill: featureColors[f.feature_type] || "#94a3b8", stroke: "#0f172a", "stroke-opacity": 0.2, "stroke-width": 0.8 }});
                rect.addEventListener("mousemove", (evt) => showHover(evt, f));
                rect.addEventListener("mouseleave", hideHover);
                svg.appendChild(rect);
                if (labelsOn && w >= 8) drawLabelCallout(svg, fx0, fx1, y, h, featureTag(f), f.feature_type === "exon");
            }});

            svg.appendChild(mkSvg("line", {{ x1: fusedJunctionX, y1: y - 30, x2: fusedJunctionX, y2: y + 30, stroke: "#0ea5e9", "stroke-width": 2.4, "stroke-dasharray": "5 4" }}));
            const jText = mkSvg("text", {{ x: fusedJunctionX + 6, y: y - 34, "font-size": 11, fill: "#0ea5e9", "font-weight": 700 }});
            jText.textContent = "Fusion junction";
            svg.appendChild(jText);
        }}

        function renderAll() {{
            const zoomPairFactor = 1 + Math.max(0, parseInt(zoomPairEl.value, 10) / 100);
            const zoomFusedFactor = 1 + Math.max(0, parseInt(zoomFusedEl.value, 10) / 100);
            const f5 = filteredFeatures(features5, tx5El.value);
            const f3 = filteredFeatures(features3, tx3El.value);
            renderPairBlock(f5, f3, tx5El.value, tx3El.value, zoomPairFactor);
            renderFusedBlock(f5, f3, tx5El.value, tx3El.value, zoomFusedFactor);
        }}

        function applyZoom() {{
            zoomPairValueEl.textContent = String(parseInt(zoomPairEl.value, 10));
            zoomFusedValueEl.textContent = String(parseInt(zoomFusedEl.value, 10));
            renderAll();
        }}

        function attachWheelZoom() {{
            pairWrapEl.addEventListener("wheel", (evt) => {{
                evt.preventDefault();
                const current = parseInt(zoomPairEl.value, 10);
                const step = evt.deltaY < 0 ? 10 : -10;
                const next = clamp(current + step, parseInt(zoomPairEl.min, 10), parseInt(zoomPairEl.max, 10));
                if (next === current) return;
                zoomPairEl.value = String(next);
                applyZoom();
            }}, {{ passive: false }});

            fusedWrapEl.addEventListener("wheel", (evt) => {{
                evt.preventDefault();
                const current = parseInt(zoomFusedEl.value, 10);
                const step = evt.deltaY < 0 ? 10 : -10;
                const next = clamp(current + step, parseInt(zoomFusedEl.min, 10), parseInt(zoomFusedEl.max, 10));
                if (next === current) return;

                const fusedSvg = document.getElementById("fusedView");
                const cx = clamp(eventToCanvasX(evt, fusedSvg), fusedX0, fusedX1);
                const leftRatio = clamp((cx - fusedX0) / Math.max(1, fusedJunctionX - fusedX0), 0, 1);
                const rightRatio = clamp((cx - fusedJunctionX) / Math.max(1, fusedX1 - fusedJunctionX), 0, 1);
                const f5 = filteredFeatures(features5, tx5El.value);
                const f3 = filteredFeatures(features3, tx3El.value);
                const dir5 = inferDirection(f5);
                const dir3 = inferDirection(f3);
                const keepLower5 = dir5 > 0;
                const keepLower3 = dir3 < 0;
                panStateFused.g5 = adjustedPanForRetainedZoom(f5, bp5, tx5El.value, current, next, panStateFused.g5, leftRatio, keepLower5);
                panStateFused.g3 = adjustedPanForRetainedZoom(f3, bp3, tx3El.value, current, next, panStateFused.g3, rightRatio, keepLower3);

                zoomFusedEl.value = String(next);
                applyZoom();
            }}, {{ passive: false }});
        }}

        function eventToCanvasX(evt, svg) {{
            const rect = svg.getBoundingClientRect();
            const rel = (evt.clientX - rect.left) / Math.max(1, rect.width);
            return rel * canvasW;
        }}

        function trackFromX(x) {{
            return x < pairSplitX ? "g5" : "g3";
        }}

        function attachPanHandlers() {{
            const pairSvg = document.getElementById("pairView");
            const fusedSvg = document.getElementById("fusedView");

            function startDrag(evt, svg, isFused) {{
                const x = eventToCanvasX(evt, svg);
                const target = isFused ? (x < fusedJunctionX ? "g5" : "g3") : trackFromX(x);
                const state = isFused ? panStateFused : panStatePair;
                dragState = {{ target, startX: evt.clientX, startPan: state[target], isFused }};
                (isFused ? fusedWrapEl : pairWrapEl).classList.add("dragging");
            }}

            pairSvg.addEventListener("mousedown", (evt) => startDrag(evt, pairSvg, false));
            fusedSvg.addEventListener("mousedown", (evt) => startDrag(evt, fusedSvg, true));

            window.addEventListener("mousemove", (evt) => {{
                if (!dragState) return;
                const dx = evt.clientX - dragState.startX;
                const deltaNorm = dx / 320;
                const state = dragState.isFused ? panStateFused : panStatePair;
                state[dragState.target] = clamp(dragState.startPan - deltaNorm, -1, 1);
                renderAll();
            }});

            window.addEventListener("mouseup", () => {{
                pairWrapEl.classList.remove("dragging");
                fusedWrapEl.classList.remove("dragging");
                dragState = null;
            }});
        }}

        fillSelect(tx5El, transcriptOptions5);
        fillSelect(tx3El, transcriptOptions3);
        setPreferredTranscript(tx5El, defaultTx5);
        setPreferredTranscript(tx3El, defaultTx3);

        tx5El.addEventListener("change", renderAll);
        tx3El.addEventListener("change", renderAll);
        zoomPairEl.addEventListener("input", applyZoom);
        zoomFusedEl.addEventListener("input", applyZoom);
        toggleLabelsBtn.addEventListener("click", () => {{
            labelsOn = !labelsOn;
            toggleLabelsBtn.textContent = labelsOn ? "Annotation labels: ON" : "Annotation labels: OFF";
            toggleLabelsBtn.classList.toggle("active", labelsOn);
            renderAll();
        }});

        document.querySelectorAll("input[data-ftype]").forEach((cb) => {{
            cb.addEventListener("change", () => {{
                const type = cb.getAttribute("data-ftype");
                if (cb.checked) enabledTypes.add(type); else enabledTypes.delete(type);
                renderAll();
            }});
        }});

        renderAll();
        attachWheelZoom();
        attachPanHandlers();
    </script>
</body>
</html>
"""


def _sanitize_filename(text: str) -> str:
    safe = "".join(ch if ch.isalnum() or ch in {"_", "-", "."} else "_" for ch in text)
    return safe.strip("_") or "fusion"


def _build_plot_html(
    row: pd.Series,
    gene_feature_df: pd.DataFrame,
 ) -> str | None:
    gene_pair = row.get("5'-3'Gene_Partners")
    coord5 = row.get("5'co-ordinate")
    coord3 = row.get("3'co-ordinate")

    gene1, gene2 = _parse_gene_pair(gene_pair)
    bp5 = _parse_coord(coord5)
    bp3 = _parse_coord(coord3)

    if not gene1 or not gene2 or bp5 is None or bp3 is None:
        return None

    g1_df = gene_feature_df[gene_feature_df["gene"].astype(str) == str(gene1)]
    g2_df = gene_feature_df[gene_feature_df["gene"].astype(str) == str(gene2)]
    if g1_df.empty or g2_df.empty:
        return None

    default_tx_5 = _extract_transcript_hint(row, "5")
    default_tx_3 = _extract_transcript_hint(row, "3")

    return _render_plot_html(
        gene1,
        gene2,
        bp5,
        bp3,
        g1_df,
        g2_df,
        default_tx_5=default_tx_5,
        default_tx_3=default_tx_3,
    )


def _build_plot_page(
    row: pd.Series,
    gene_feature_df: pd.DataFrame,
    out_dir: Path,
    row_idx: int,
) -> str | None:
    plot_html = _build_plot_html(row, gene_feature_df)
    if plot_html is None:
        return None

    gene1, gene2 = _parse_gene_pair(row.get("5'-3'Gene_Partners"))
    bp5 = _parse_coord(row.get("5'co-ordinate"))
    bp3 = _parse_coord(row.get("3'co-ordinate"))
    out_name = _sanitize_filename(f"{row_idx}_{gene1}_{gene2}_{bp5[0]}_{bp5[1]}_{bp3[0]}_{bp3[1]}.html")
    out_file = out_dir / out_name
    out_file.write_text(plot_html, encoding="utf-8")
    return out_file.name


def add_fusion_plot_links(
    df: pd.DataFrame,
    report_dir: str | Path,
    mode: str = "external",
) -> tuple[pd.DataFrame, dict[str, str]]:
    output_col = "Fusion Visualization"
    resource_key = "gene_feature_table"
    logger.info("Generating interactive figures for fusion events")
    mode = str(mode or "embed").strip().lower()
    if mode not in {"embed", "external"}:
        raise ValueError(f"Unsupported fusion plot mode: {mode}")
    if resource_key not in config:
        logger.info("Fusion plot links skipped: resource key '%s' not present in config", resource_key)
        return _insert_as_eighth_column(df, output_col, "NA"), {}

    resource_path = Path(str(config[resource_key]))
    if not resource_path.is_file():
        logger.warning("Fusion plot links skipped: resource file not found: %s", resource_path)
        return _insert_as_eighth_column(df, output_col, "NA"), {}

    req_cols = {"5'-3'Gene_Partners", "5'co-ordinate", "3'co-ordinate"}
    if not req_cols.issubset(df.columns):
        logger.warning("Fusion plot links skipped: required fusion columns are missing")
        return _insert_as_eighth_column(df, output_col, "NA"), {}

    gene_feature_df = pd.read_csv(resource_path, sep="\t", compression="infer", dtype=str)
    for c in ["start", "end"]:
        gene_feature_df[c] = pd.to_numeric(gene_feature_df[c], errors="coerce")
    gene_feature_df = gene_feature_df.dropna(subset=["start", "end", "gene"]).copy()
    gene_feature_df["start"] = gene_feature_df["start"].astype(int)
    gene_feature_df["end"] = gene_feature_df["end"].astype(int)
    feature_types = {
        _normalize_feature_type(x)
        for x in gene_feature_df.get("feature_type", pd.Series(dtype=str)).dropna().astype(str).tolist()
    }
    if not (feature_types & {"five_prime_utr", "three_prime_utr", "utr"}):
        raise ValueError(
            "UTR features were not found in gene_feature_table. "
            "Regenerate resources with UTR parsing enabled before plotting fusions."
        )

    report_dir = Path(report_dir)
    plots_dir = report_dir / "fusion_plots"
    if mode == "external":
        plots_dir.mkdir(parents=True, exist_ok=True)

    output = df.copy()
    links: list[str] = []
    embedded_plots: dict[str, str] = {}

    for idx, row in output.iterrows():
        if mode == "external":
            rel_file = _build_plot_page(row, gene_feature_df, plots_dir, int(idx))
            if rel_file:
                href = f"fusion_plots/{rel_file}"
                links.append(
                    f'<a href="{href}" target="_blank" '
                    'style="display:inline-block;padding:2px 8px;border-radius:6px;'
                    'background:#e8f1ff;color:#1e40af;font-weight:600;text-decoration:none;" '
                    '>View Fusion</a>'
                )
            else:
                links.append("NA")
            continue

        plot_html = _build_plot_html(row, gene_feature_df)
        if plot_html:
            plot_key = f"fusion_plot_{int(idx)}"
            embedded_plots[plot_key] = base64.b64encode(plot_html.encode("utf-8")).decode("ascii")
            links.append(
                f'<a href="#" onclick="return openEmbeddedFusionPlot(\'{plot_key}\')" target="_blank" '
                'style="display:inline-block;padding:2px 8px;border-radius:6px;'
                'background:#e8f1ff;color:#1e40af;font-weight:600;text-decoration:none;" '
                '>View Fusion</a>'
            )
        else:
            links.append("NA")

    return _insert_as_eighth_column(output, output_col, links), embedded_plots
