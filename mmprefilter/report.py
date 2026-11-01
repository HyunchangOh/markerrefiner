"""Text and HTML report. Thresholds are read from the config; the page does not edit them."""

import html
import os

from .structure import row_matches


def _numbers(rows, key):
    values = []
    for row in rows:
        text = (row.get(key) or "").strip()
        if not text:
            continue
        try:
            values.append(float(text))
        except ValueError:
            continue
    return values


def _svg_histogram(values, cutoff, greater, title):
    width, height = 360, 140
    left, right, top, bottom = 36, 8, 22, 28
    plot_w = width - left - right
    plot_h = height - top - bottom
    if not values:
        return "<p>No numeric values for {0}.</p>".format(html.escape(title))
    low = min(min(values), cutoff)
    high = max(max(values), cutoff)
    if high == low:
        high = low + 1
    bins = 12
    span = (high - low) / float(bins)
    counts = [0] * bins
    for value in values:
        index = int((value - low) / span)
        if index == bins:
            index = bins - 1
        counts[index] += 1
    peak = max(counts) or 1
    bars = []
    for index, count in enumerate(counts):
        x = left + index * (plot_w / float(bins))
        bar_h = plot_h * (count / float(peak))
        y = top + plot_h - bar_h
        edge = low + index * span
        passes = edge > cutoff if greater else (edge + span) < cutoff
        # A bin can straddle the cutoff. Colour it if its midpoint passes.
        midpoint = edge + span / 2.0
        passes = midpoint > cutoff if greater else midpoint < cutoff
        colour = "#2f6f4e" if passes else "#c5cdd6"
        bars.append(
            '<rect x="{x:.1f}" y="{y:.1f}" width="{w:.1f}" height="{h:.1f}" fill="{c}"/>'.format(
                x=x + 1,
                y=y,
                w=max(1.0, plot_w / float(bins) - 2),
                h=bar_h,
                c=colour,
            )
        )
    cut_x = left + (cutoff - low) / float(high - low) * plot_w
    label = "{0} {1} {2:g}".format(title, ">" if greater else "<", cutoff)
    return """<figure>
  <figcaption>{label}</figcaption>
  <svg viewBox="0 0 {width} {height}" width="100%" role="img">
    <line x1="{left}" y1="{base}" x2="{right}" y2="{base}" stroke="#1c1c1c"/>
    {bars}
    <line x1="{cut:.1f}" y1="{top}" x2="{cut:.1f}" y2="{base}" stroke="#9a3412" stroke-width="1.5"/>
    <text x="{left}" y="{height_minus}" font-size="11" fill="#333">{low:.0f}</text>
    <text x="{right_label}" y="{height_minus}" font-size="11" fill="#333" text-anchor="end">{high:.0f}</text>
  </svg>
</figure>""".format(
        label=html.escape(label),
        width=width,
        height=height,
        left=left,
        base=top + plot_h,
        right=left + plot_w,
        bars="\n    ".join(bars),
        cut=cut_x,
        top=top,
        height_minus=height - 8,
        low=low,
        high=high,
        right_label=left + plot_w,
    )


METRICS = [
    ("n_exons", "n_exons_gt", True),
    ("n_introns", "n_introns_lt", False),
    ("max_exon_bp", "max_exon_bp_gt", True),
    ("max_intron_bp", "max_intron_bp_lt", False),
    ("avg_exon_bp", "avg_exon_bp_gt", True),
    ("avg_intron_bp", "avg_intron_bp_lt", False),
]


def pass_counts(rows, thresholds):
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
        "Thresholds come from the config file. A blank metric does not fail a row.",
        "Cutoffs: {0}".format(thresholds),
        "",
    ]
    text_lines.extend(qc_lines)
    text_lines.append("")
    sections = []
    for run_name, rows in run_rows:
        total, n_loci, kept, kept_loci = pass_counts(rows, thresholds)
        text_lines.append(
            "{0}: {1} rows, {2} loci; passing: {3} rows, {4} loci".format(
                run_name, total, n_loci, kept, kept_loci
            )
        )
        figures = []
        for column, key, greater in METRICS:
            figures.append(
                _svg_histogram(
                    _numbers(rows, column),
                    float(thresholds[key]),
                    greater,
                    column,
                )
            )
        sections.append(
            "<section><h2>{0}</h2><p>{1} rows ({2} loci). Passing all six cutoffs: {3} rows, {4} loci.</p><div class=\"grid\">{5}</div></section>".format(
                html.escape(run_name),
                total,
                n_loci,
                kept,
                kept_loci,
                "\n".join(figures),
            )
        )
    if comparison_rows:
        text_lines.append("")
        text_lines.append("GFF genomic lengths versus alignment-mask lengths:")
        differ_intron = 0
        coding_agree = 0
        compared = 0
        for row in comparison_rows:
            if not row.get("gff_cds_bp") or not row.get("msa_coding_bp"):
                continue
            compared += 1
            if row.get("coding_bp_match") == "yes":
                coding_agree += 1
            if row.get("gff_max_intron_bp") != row.get("msa_max_intron_bp"):
                differ_intron += 1
        sentence = (
            "Compared {0} loci. Coding length agrees for {1}. "
            "Longest-intron length differs for {2}. "
            "GFF intron lengths are the genomic distances between CDS blocks. "
            "Alignment intron lengths are the N runs on the reference row. "
            "The filter uses the alignment, because that is the sequence opened in Geneious."
        ).format(compared, coding_agree, differ_intron)
        text_lines.append(sentence)
        sections.append("<section><h2>Why the filter uses the alignment</h2><p>{0}</p></section>".format(html.escape(sentence)))

    text_path = os.path.join(output_dir, "report.txt")
    with open(text_path, "w") as handle:
        handle.write("\n".join(text_lines) + "\n")

    html_path = os.path.join(output_dir, "report.html")
    document = """<!DOCTYPE html>
<html lang="en">
<head>
<meta charset="utf-8"/>
<title>mmprefilter report</title>
<style>
  body {{ font-family: Georgia, serif; margin: 2rem auto; max-width: 980px; color: #1c1c1c; line-height: 1.45; }}
  h1 {{ font-size: 1.6rem; }}
  .grid {{ display: grid; grid-template-columns: 1fr 1fr; gap: 0.5rem 1rem; }}
  figure {{ margin: 0; }}
  figcaption {{ font-family: ui-sans-serif, sans-serif; font-size: 0.85rem; }}
  pre {{ font-family: ui-monospace, monospace; font-size: 0.85rem; background: #f4f1ea; padding: 0.8rem; }}
</style>
</head>
<body>
<h1>mmprefilter</h1>
<p>Green bins pass the cutoff drawn in red. The cutoff is the value in the config file.</p>
<pre>{pre}</pre>
{sections}
</body>
</html>
""".format(
        pre=html.escape("\n".join(text_lines)),
        sections="\n".join(sections),
    )
    with open(html_path, "w") as handle:
        handle.write(document)
    return text_path, html_path
