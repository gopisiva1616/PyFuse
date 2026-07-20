import base64
import json
import re
from pathlib import Path
from importlib.resources import files
import pandas as pd
from jinja2 import Template, Environment, select_autoescape
from .common_utils import utils, config

def _read_text_asset(path: str) -> str:
    with open(path, "r", encoding="utf-8") as file_handle:
        return file_handle.read()


def _encode_file_as_data_url(path: str, mime_type: str) -> str:
    with open(path, "rb") as file_handle:
        encoded = base64.b64encode(file_handle.read()).decode("ascii")
    return f"data:{mime_type};base64,{encoded}"


def _safe_inline_js(js_text: str) -> str:
    return re.sub(r"</script", r"<\\/script", js_text, flags=re.IGNORECASE)


def _safe_inline_json(json_text: str) -> str:
    return json_text.replace("</", "<\\/")


def _resolve_report_logo_path() -> str | None:
    configured_logo = config.get("report_logo")
    if not configured_logo:
        return None

    configured_path = Path(str(configured_logo)).expanduser()
    if configured_path.is_file():
        return str(configured_path)

    settings_dir = Path(str(files("pyfuse.config")))
    resolved_from_settings = (settings_dir / configured_path).resolve()
    if resolved_from_settings.is_file():
        return str(resolved_from_settings)

    return None

def _apply_frame_status_styling(df: pd.DataFrame) -> list:
    """Apply CSS styling to Frame_Status column values."""
    data = df.values.tolist()
    if 'Frame_Status' in df.columns:
        frame_col_idx = df.columns.get_loc('Frame_Status')
        for row in data:
            if frame_col_idx < len(row):
                val = row[frame_col_idx]
                if isinstance(val, str):
                    if 'In-Frame' in val:
                        row[frame_col_idx] = f'<span class="cell-in-frame">{val}</span>'
                    elif 'Out-of-Frame' in val:
                        row[frame_col_idx] = f'<span class="cell-out-of-frame">{val}</span>'
    return data

def write_df2html(
    df: pd.DataFrame,
    output_html: str,
    title: str,
    embedded_fusion_plots: dict[str, str] | None = None,
):

    embedded_fusion_plots = embedded_fusion_plots or {}

    logo_path = _resolve_report_logo_path()

    assets = {
        "font_regular_data": _encode_file_as_data_url(config["font_regular"], "font/woff2"),
        "font_bold_data": _encode_file_as_data_url(config["font_bold"], "font/woff2"),
        "jquery_js_text": _safe_inline_js(_read_text_asset(config["jquery"])),
        "datatables_js_text": _safe_inline_js(_read_text_asset(config["datatables_js"])),
        "datatables_css_text": _read_text_asset(config["datatables_css"]),
        "report_logo_data": _encode_file_as_data_url(logo_path, "image/png") if logo_path else None,
        "embedded_fusion_plots_json": _safe_inline_json(json.dumps(embedded_fusion_plots)),
    }

    env = Environment(autoescape=False)
    html_template = env.from_string("""
    <!DOCTYPE html>
    <html>
    <head>
        <meta charset="utf-8">
        <title>PyFuse Results</title>
        <style>
            @font-face {
                font-family: 'Inter';
                font-weight: 400;
                src: url('{{ font_regular_data }}') format('woff2');
            }
            @font-face {
                font-family: 'Inter';
                font-weight: 600;
                src: url('{{ font_bold_data }}') format('woff2');
            }
            body {
                font-family: 'Inter', sans-serif;
                padding: 30px;
                background-color: #f6f8fa;
                color: #333;
            }
            h2 {
                text-align: center;
                font-weight: 600;
                color: #1f2937;
                margin-bottom: 6px;
            }
            .report-header {
                text-align: center;
                margin-bottom: 10px;
            }
            .report-logo {
                display: block;
                height: auto;
                margin: 0 auto;
                width: 18%;
                transform: scaleY(0.93);
                transform-origin: top center;
            }
            table.dataTable {
                border-collapse: collapse;
                box-shadow: 0 2px 12px rgba(0,0,0,0.05);
                border-radius: 10px;
                overflow: hidden;
                background: white;
            }
            table.dataTable thead th {
                background-color: #66a7ed;
                color: #1f2937;
                font-weight: 600;
                padding: 10px;
                border-bottom: 1px solid #d1d5db;
                border-right: 1px solid #dce6f2;
            }
            table.dataTable thead th:last-child { border-right: none; }
            table.dataTable td {
                padding: 8px 10px;
                border-bottom: 1px solid #f1f1f1;
            }
            table.dataTable tbody td {
                border-right: 1px solid #eef3f8;
            }
            table.dataTable tbody td:last-child { border-right: none; }
            table.dataTable tbody tr td:nth-child(even) {
                background-color: #fbfdff;
            }
            a {
                color: #4B7EC4;
                text-decoration: underline;
            }
            a:hover {
                color: #2C5AA0;
            }
            thead input {
                width: 100%;
                box-sizing: border-box;
                padding: 6px 8px;
                font-size: 12px;
                border: 1px solid #ccc;
                border-radius: 5px;
                background-color: #f9fafb;
            }
            .dataTables_wrapper .dataTables_length,
            .dataTables_wrapper .dataTables_filter,
            .dataTables_wrapper .dataTables_info,
            .dataTables_wrapper .dataTables_paginate {
                font-size: 13px;
                margin-top: 10px;
            }
            .cell-in-frame {
                background-color: #d4edda !important;
            }
            .cell-out-of-frame {
                background-color: #f8d7da !important;
            }
            {{ datatables_css_text }}
        </style>
    </head>
    <body>
        <div class="report-header">
            {% if report_logo_data %}
            <img class="report-logo" src="{{ report_logo_data }}" alt="PyFuse logo" />
            {% endif %}
        </div>
        <table id="data-table" class="display" style="width:100%">
            <thead>
                <tr>
                    {% for col in columns %}
                    <th>{{ col }}</th>
                    {% endfor %}
                </tr>
                <tr>
                    {% for col in columns %}
                    <th><input type="text" placeholder="Search {{ col }}" /></th>
                    {% endfor %}
                </tr>
            </thead>
            <tbody>
                {% for row in data %}
                <tr>
                    {% for val in row %}
                    <td>{{ val }}</td>
                    {% endfor %}
                </tr>
                {% endfor %}
            </tbody>
        </table>

        <script>{{ jquery_js_text }}</script>
        <script>{{ datatables_js_text }}</script>
        <script>
            const embeddedFusionPlots = {{ embedded_fusion_plots_json }};

            function decodeBase64Utf8(encoded) {
                const binary = window.atob(encoded);
                const bytes = Uint8Array.from(binary, (ch) => ch.charCodeAt(0));
                return new TextDecoder('utf-8').decode(bytes);
            }

            function openEmbeddedFusionPlot(plotKey) {
                const encoded = embeddedFusionPlots[plotKey];
                if (!encoded) {
                    return false;
                }

                const popup = window.open('', '_blank');
                if (!popup) {
                    return false;
                }

                popup.document.open();
                popup.document.write(decodeBase64Utf8(encoded));
                popup.document.close();
                return false;
            }

            $(document).ready(function () {
                var table = $('#data-table').DataTable({
                    orderCellsTop: true,
                    fixedHeader: true,
                    pageLength: 30
                });

                $('#data-table thead tr:eq(1) th').each(function (i) {
                    $('input', this).on('keyup change clear', function () {
                        if (table.column(i).search() !== this.value) {
                            table.column(i).search(this.value).draw();
                        }
                    });
                });
            });
        </script>
    </body>
    </html>
    """)

    html_str = html_template.render(
        title=title,
        columns=list(df.columns),
        data=_apply_frame_status_styling(df),
        **assets
    )

    with open(output_html, "w", encoding="utf-8") as f:
        f.write(html_str)
