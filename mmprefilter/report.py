"""Interactive report. Sliders choose the six cutoffs and can be saved to a file."""

import json
import os

from .config import THRESHOLD_KEYS
from .structure import row_matches


KEEP_FIELDS = [
    "locus",
    "n_exons",
    "n_introns",
    "max_exon_bp",
    "max_intron_bp",
    "avg_exon_bp",
    "avg_intron_bp",
]


def _slim(rows):
    slim = []
    for row in rows:
        slim.append({key: row.get(key, "") for key in KEEP_FIELDS})
    return slim


def _pass_counts(rows, thresholds):
    matched = [row for row in rows if row_matches(row, thresholds)]
    loci = {(row.get("locus") or "").strip() for row in matched}
    loci.discard("")
    all_loci = {(row.get("locus") or "").strip() for row in rows}
    all_loci.discard("")
    return len(rows), len(all_loci), len(matched), len(loci)


def write_report(output_dir, thresholds, run_rows, comparison_rows, qc_lines):
    text_lines = [
        "mmprefilter report",
        "",
        "Open report.html to move the cutoffs. Download thresholds.json from that",
        "page, then set filter.thresholds_file to the saved file and rerun from",
        "the filter stage. A thresholds_file replaces filter.thresholds.",
        "",
    ]
    if thresholds:
        text_lines.append("Cutoffs currently in the config: {0}".format(thresholds))
    else:
        text_lines.append("No cutoffs in the config yet. Choose them in report.html.")
    text_lines.append("")
    text_lines.extend(qc_lines)
    text_lines.append("")
    if thresholds:
        for run_name, rows in run_rows:
            total, n_loci, kept, kept_loci = _pass_counts(rows, thresholds)
            text_lines.append(
                "{0}: {1} rows, {2} loci; passing the config cutoffs: {3} rows, {4} loci".format(
                    run_name, total, n_loci, kept, kept_loci
                )
            )
    if comparison_rows:
        compared = coding_agree = differ_intron = 0
        for row in comparison_rows:
            if not row.get("gff_cds_bp") or not row.get("msa_coding_bp"):
                continue
            compared += 1
            if row.get("coding_bp_match") == "yes":
                coding_agree += 1
            if row.get("gff_max_intron_bp") != row.get("msa_max_intron_bp"):
                differ_intron += 1
        text_lines.extend(
            [
                "",
                "Compared {0} loci. Coding length agrees for {1}. Longest-intron length differs for {2}.".format(
                    compared, coding_agree, differ_intron
                ),
                "GFF intron lengths are genomic distances. Alignment intron lengths are N runs.",
                "The filter uses the alignment.",
            ]
        )

    text_path = os.path.join(output_dir, "report.txt")
    with open(text_path, "w") as handle:
        handle.write("\n".join(text_lines) + "\n")

    payload = {
        "thresholds": {key: thresholds[key] for key in THRESHOLD_KEYS} if thresholds else None,
        "datasets": [{"name": name, "rows": _slim(rows)} for name, rows in run_rows],
    }
    html_path = os.path.join(output_dir, "report.html")
    with open(html_path, "w") as handle:
        handle.write(_page(payload, "\n".join(text_lines)))
    return text_path, html_path


def _page(payload, preamble):
    data = json.dumps(payload).replace("<", "\\u003c")
    preamble = (
        preamble.replace("&", "&amp;").replace("<", "&lt;").replace(">", "&gt;")
    )
    return """<!DOCTYPE html>
<html lang="en">
<head>
<meta charset="utf-8"/>
<title>mmprefilter cutoffs</title>
<style>
  body {{ font-family: Georgia, serif; margin: 1.5rem auto; max-width: 1100px; color: #1c1c1c; line-height: 1.45; }}
  h1 {{ font-size: 1.5rem; }}
  h2 {{ font-size: 1.15rem; margin-top: 2rem; }}
  .note {{ color: #444; max-width: 50rem; }}
  .slider-row {{ margin: 0.8rem 0 1rem; }}
  .slider-row label {{ display: block; font-family: ui-sans-serif, sans-serif; font-size: 0.9rem; }}
  .slider-row input[type="range"] {{ width: 100%; }}
  button {{ font: 1rem ui-sans-serif, sans-serif; padding: 0.45rem 0.8rem; cursor: pointer; }}
  .counts {{ display: grid; grid-template-columns: repeat(auto-fit, minmax(220px, 1fr)); gap: 0.75rem; }}
  .card {{ border: 1px solid #ccc; padding: 0.8rem 1rem; background: #f7f4ee; }}
  .card strong {{ font-size: 1.25rem; }}
  .grid {{ display: grid; grid-template-columns: repeat(auto-fit, minmax(280px, 1fr)); gap: 0.5rem 1rem; }}
  figure {{ margin: 0; }}
  figcaption {{ font-family: ui-sans-serif, sans-serif; font-size: 0.8rem; }}
  pre {{ font-family: ui-monospace, monospace; font-size: 0.8rem; background: #f4f1ea; padding: 0.8rem; white-space: pre-wrap; }}
</style>
</head>
<body>
<h1>Choose the cutoffs</h1>
<p class="note">
  These plots are the alignment tables this run just wrote. Move a slider and the
  counts update. A row is kept only when every slider passes. An empty cell does
  not fail the row. Download the file, set <code>filter.thresholds_file</code> to
  that path, and rerun from the filter stage. You do not type the six numbers.
</p>
<pre>{preamble}</pre>
<h2>Sliders</h2>
<div id="sliders"></div>
<p><button id="download" type="button">Download thresholds.json</button></p>
<pre id="json"></pre>
<h2>Passing rows</h2>
<div id="counts" class="counts"></div>
<h2>Distributions</h2>
<div id="charts"></div>
<script>
var DATA = {data};
(function () {{
  var METRICS = [
    {{ key: "n_exons", cutoff: "n_exons_gt", op: "gt", integer: true, label: "n_exons >" }},
    {{ key: "n_introns", cutoff: "n_introns_lt", op: "lt", integer: true, label: "n_introns <" }},
    {{ key: "max_exon_bp", cutoff: "max_exon_bp_gt", op: "gt", integer: true, label: "max_exon_bp >" }},
    {{ key: "max_intron_bp", cutoff: "max_intron_bp_lt", op: "lt", integer: true, label: "max_intron_bp <" }},
    {{ key: "avg_exon_bp", cutoff: "avg_exon_bp_gt", op: "gt", integer: false, label: "avg_exon_bp >" }},
    {{ key: "avg_intron_bp", cutoff: "avg_intron_bp_lt", op: "lt", integer: false, label: "avg_intron_bp <" }}
  ];

  function num(row, key) {{
    var text = row[key];
    if (text === undefined || text === null || String(text).trim() === "") return NaN;
    var value = parseFloat(text);
    return isFinite(value) ? value : NaN;
  }}

  function rowMatches(row, t) {{
    for (var i = 0; i < METRICS.length; i++) {{
      var metric = METRICS[i];
      var value = num(row, metric.key);
      if (!isFinite(value)) continue;
      if (metric.op === "gt" && !(value > t[metric.cutoff])) return false;
      if (metric.op === "lt" && !(value < t[metric.cutoff])) return false;
    }}
    return true;
  }}

  function rangeFor(key) {{
    var low = Infinity, high = -Infinity, any = false;
    DATA.datasets.forEach(function (ds) {{
      ds.rows.forEach(function (row) {{
        var value = num(row, key);
        if (!isFinite(value)) return;
        any = true;
        if (value < low) low = value;
        if (value > high) high = value;
      }});
    }});
    if (!any) return {{ min: 0, max: 1 }};
    if (low === high) {{ low -= 1; high += 1; }}
    return {{ min: low, max: high }};
  }}

  var thresholds = {{}};
  METRICS.forEach(function (metric) {{
    var span = rangeFor(metric.key);
    var pad = metric.integer ? 1 : Math.max((span.max - span.min) * 0.05, 0.01);
    var min = metric.op === "gt" ? span.min - pad : span.min;
    var max = metric.op === "lt" ? span.max + pad : span.max;
    var initial = metric.op === "gt" ? min : max;
    if (DATA.thresholds && DATA.thresholds[metric.cutoff] !== undefined) {{
      initial = Number(DATA.thresholds[metric.cutoff]);
      if (initial < min) min = initial;
      if (initial > max) max = initial;
    }}
    metric.min = min;
    metric.max = max;
    metric.step = metric.integer ? 1 : Math.max((max - min) / 500, 0.01);
    thresholds[metric.cutoff] = initial;
  }});

  var canvases = [];
  var sliderBox = document.getElementById("sliders");
  METRICS.forEach(function (metric) {{
    var row = document.createElement("div");
    row.className = "slider-row";
    var label = document.createElement("label");
    var readout = document.createElement("span");
    label.appendChild(document.createTextNode(metric.label + " "));
    label.appendChild(readout);
    var input = document.createElement("input");
    input.type = "range";
    input.min = String(metric.min);
    input.max = String(metric.max);
    input.step = String(metric.step);
    input.value = String(thresholds[metric.cutoff]);
    function show() {{
      var value = Number(input.value);
      if (metric.integer) value = Math.round(value);
      thresholds[metric.cutoff] = value;
      readout.textContent = metric.integer ? String(value) : value.toFixed(2);
      draw();
    }}
    input.addEventListener("input", show);
    row.appendChild(label);
    row.appendChild(input);
    sliderBox.appendChild(row);
    show();
  }});

  var chartBox = document.getElementById("charts");
  DATA.datasets.forEach(function (ds) {{
    var section = document.createElement("section");
    var heading = document.createElement("h2");
    heading.textContent = ds.name;
    section.appendChild(heading);
    var grid = document.createElement("div");
    grid.className = "grid";
    METRICS.forEach(function (metric) {{
      var figure = document.createElement("figure");
      var caption = document.createElement("figcaption");
      caption.textContent = metric.key;
      var canvas = document.createElement("canvas");
      canvas.width = 360;
      canvas.height = 140;
      figure.appendChild(caption);
      figure.appendChild(canvas);
      grid.appendChild(figure);
      canvases.push({{ canvas: canvas, metric: metric, rows: ds.rows }});
    }});
    section.appendChild(grid);
    chartBox.appendChild(section);
  }});

  function drawHistogram(canvas, values, cutoff, greater) {{
    var ctx = canvas.getContext("2d");
    var w = canvas.width, h = canvas.height;
    ctx.clearRect(0, 0, w, h);
    var left = 28, right = 8, top = 8, bottom = 22;
    var plotW = w - left - right, plotH = h - top - bottom;
    if (!values.length) return;
    var low = Math.min(Math.min.apply(null, values), cutoff);
    var high = Math.max(Math.max.apply(null, values), cutoff);
    if (low === high) high = low + 1;
    var bins = 16;
    var span = (high - low) / bins;
    var counts = [];
    for (var i = 0; i < bins; i++) counts.push(0);
    values.forEach(function (value) {{
      var index = Math.floor((value - low) / span);
      if (index >= bins) index = bins - 1;
      if (index < 0) index = 0;
      counts[index] += 1;
    }});
    var peak = Math.max.apply(null, counts) || 1;
    counts.forEach(function (count, index) {{
      var midpoint = low + (index + 0.5) * span;
      var pass = greater ? midpoint > cutoff : midpoint < cutoff;
      ctx.fillStyle = pass ? "#2f6f4e" : "#c5cdd6";
      var barH = plotH * (count / peak);
      var x = left + index * (plotW / bins);
      ctx.fillRect(x + 1, top + plotH - barH, Math.max(1, plotW / bins - 2), barH);
    }});
    var cutX = left + ((cutoff - low) / (high - low)) * plotW;
    ctx.strokeStyle = "#9a3412";
    ctx.beginPath();
    ctx.moveTo(cutX, top);
    ctx.lineTo(cutX, top + plotH);
    ctx.stroke();
    ctx.strokeStyle = "#1c1c1c";
    ctx.beginPath();
    ctx.moveTo(left, top + plotH);
    ctx.lineTo(left + plotW, top + plotH);
    ctx.stroke();
  }}

  function draw() {{
    canvases.forEach(function (item) {{
      var values = [];
      item.rows.forEach(function (row) {{
        var value = num(row, item.metric.key);
        if (isFinite(value)) values.push(value);
      }});
      drawHistogram(
        item.canvas,
        values,
        thresholds[item.metric.cutoff],
        item.metric.op === "gt"
      );
    }});
    var counts = document.getElementById("counts");
    counts.innerHTML = "";
    DATA.datasets.forEach(function (ds) {{
      var loci = {{}};
      var kept = 0;
      var keptLoci = {{}};
      ds.rows.forEach(function (row) {{
        if (row.locus) loci[row.locus] = true;
        if (!rowMatches(row, thresholds)) return;
        kept += 1;
        if (row.locus) keptLoci[row.locus] = true;
      }});
      var card = document.createElement("div");
      card.className = "card";
      card.innerHTML = "<div>" + ds.name + "</div><strong>" +
        Object.keys(keptLoci).length + " loci</strong><div>" +
        kept + " of " + ds.rows.length + " rows pass</div><div>" +
        Object.keys(loci).length + " loci in the table</div>";
      counts.appendChild(card);
    }});
    var exported = {{}};
    METRICS.forEach(function (metric) {{
      var value = thresholds[metric.cutoff];
      exported[metric.cutoff] = metric.integer ? Math.round(value) : Math.round(value * 100000) / 100000;
    }});
    document.getElementById("json").textContent = JSON.stringify(exported, null, 2);
    document.getElementById("download").onclick = function () {{
      var blob = new Blob([JSON.stringify(exported, null, 2) + "\\n"], {{ type: "application/json" }});
      var link = document.createElement("a");
      link.href = URL.createObjectURL(blob);
      link.download = "thresholds.json";
      document.body.appendChild(link);
      link.click();
      document.body.removeChild(link);
    }};
  }}
  draw();
}})();
</script>
</body>
</html>
""".format(preamble=preamble, data=data)
