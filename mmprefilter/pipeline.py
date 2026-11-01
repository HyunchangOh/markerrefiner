"""Run the stages named in the config, in order."""

import os

from .config import (
    STAGES,
    ConfigError,
    enabled,
    require,
    resolve_thresholds,
    section,
)
from .qc import (
    read_species_map,
    species_map_from_list,
    write_busco_summary,
    write_longest_isoforms,
    write_species_map,
    write_transcriptome_stats,
)
from .report import write_report
from .structure import (
    MSA_FIELDS,
    alignment_records,
    compare_masked_to_spliced,
    exon_qc_rows,
    filter_run,
    gff_metric_rows,
    iter_alignments,
    load_transcript_map,
    loci_in_alignment_dir,
    msa_metric_rows,
    read_table,
    reference_block_summary,
    write_table,
)


def run_config(cfg, start=None, end=None):
    if start is None:
        start = STAGES[0]
    if end is None:
        end = STAGES[-1]
    if start not in STAGES or end not in STAGES:
        raise ConfigError(
            "Stage must be one of: {0}".format(", ".join(STAGES))
        )
    start_i = STAGES.index(start)
    end_i = STAGES.index(end)
    if end_i < start_i:
        raise ConfigError("--from is after --to.")
    selected = set(STAGES[start_i : end_i + 1])
    output_dir = cfg["output_dir"]
    if not os.path.isdir(output_dir):
        os.makedirs(output_dir)
    suffix = cfg["alignment_suffix"]

    if "transcriptome_qc" in selected and enabled(cfg, "transcriptome_qc"):
        fasta_dir = require(section(cfg, "transcriptome_qc"), "fasta_dir", "transcriptome_qc")
        _must_dir(fasta_dir)
        write_transcriptome_stats(
            fasta_dir, os.path.join(output_dir, "transcriptome_stats.tsv")
        )

    if "busco_summary" in selected and enabled(cfg, "busco_summary"):
        busco_root = require(section(cfg, "busco_summary"), "busco_root", "busco_summary")
        _must_dir(busco_root)
        write_busco_summary(busco_root, os.path.join(output_dir, "busco_summary.tsv"))

    if "longest_isoform" in selected and enabled(cfg, "longest_isoform"):
        fasta_dir = require(section(cfg, "longest_isoform"), "fasta_dir", "longest_isoform")
        _must_dir(fasta_dir)
        write_longest_isoforms(fasta_dir, os.path.join(output_dir, "longest_isoform"))

    species_map_path = os.path.join(output_dir, "species_map.tsv")
    if "species_map" in selected and enabled(cfg, "species_map"):
        list_path = require(section(cfg, "species_map"), "transcriptome_list", "species_map")
        _must_file(list_path)
        write_species_map(species_map_path, species_map_from_list(list_path))

    if "reference_cds_check" in selected and enabled(cfg, "reference_cds_check"):
        block = section(cfg, "reference_cds_check")
        rows = compare_masked_to_spliced(
            require(block, "masked_cds_fasta", "reference_cds_check"),
            require(block, "spliced_cds_fasta", "reference_cds_check"),
            require(block, "alignment_dir", "reference_cds_check"),
            suffix,
        )
        write_table(
            os.path.join(output_dir, "reference_cds_check.tsv"),
            ["locus", "status", "transcript_id", "masked_coding_bp", "spliced_bp"],
            rows,
        )

    if "gff_metrics" in selected and enabled(cfg, "gff_metrics"):
        block = section(cfg, "gff_metrics")
        gff_path = require(block, "gff", "gff_metrics")
        alignment_dir = require(block, "alignment_dir", "gff_metrics")
        _must_file(gff_path)
        _must_dir(alignment_dir)
        map_path = block.get("transcript_map") or os.path.join(
            output_dir, "reference_cds_check.tsv"
        )
        _must_file(map_path)
        transcript_by_locus = load_transcript_map(map_path)
        rows = gff_metric_rows(
            gff_path,
            transcript_by_locus,
            loci_in_alignment_dir(alignment_dir, suffix),
        )
        write_table(
            os.path.join(output_dir, "gff_metrics.tsv"),
            [
                "locus",
                "transcript_id",
                "n_exons",
                "n_introns",
                "max_exon_bp",
                "max_intron_bp",
                "avg_exon_bp",
                "avg_intron_bp",
                "cds_bp",
            ],
            rows,
        )

    if "exon_qc" in selected and enabled(cfg, "exon_qc"):
        for run in _runs(section(cfg, "exon_qc"), "exon_qc"):
            _must_dir(run["alignment_dir"])
            rows = exon_qc_rows(run["alignment_dir"], suffix)
            write_table(
                os.path.join(output_dir, "exon_qc__{0}.tsv".format(run["name"])),
                [
                    "locus",
                    "ingroup_id",
                    "gap_frac_exonic",
                    "ident_vs_ref_exonic",
                    "n_both_ACGT",
                    "n_exonic_cols",
                ],
                rows,
            )

    if "msa_metrics" in selected and enabled(cfg, "msa_metrics"):
        for run in _runs(section(cfg, "msa_metrics"), "msa_metrics"):
            _must_dir(run["alignment_dir"])
            rows, skipped = msa_metric_rows(run["alignment_dir"], suffix)
            write_table(
                os.path.join(output_dir, "msa_metrics__{0}.csv".format(run["name"])),
                MSA_FIELDS,
                rows,
            )
            if skipped:
                with open(
                    os.path.join(output_dir, "msa_skipped__{0}.txt".format(run["name"])),
                    "w",
                ) as handle:
                    handle.write("\n".join(skipped) + "\n")

    if "compare_gff_msa" in selected and enabled(cfg, "compare_gff_msa"):
        gff_table = os.path.join(output_dir, "gff_metrics.tsv")
        _must_file(gff_table)
        alignment_dir = require(
            section(cfg, "gff_metrics"), "alignment_dir", "gff_metrics"
        )
        _must_dir(alignment_dir)
        by_locus = {}
        for locus, path in iter_alignments(alignment_dir, suffix):
            records = alignment_records(path)
            if not records:
                continue
            by_locus[locus] = reference_block_summary(records[-1][1])
        comparison = []
        for row in read_table(gff_table):
            reference = by_locus.get(row["locus"])
            if not reference:
                continue
            gff_cds = row.get("cds_bp") or ""
            msa_coding = reference.get("coding_bp", "")
            comparison.append(
                {
                    "locus": row["locus"],
                    "transcript_id": row.get("transcript_id", ""),
                    "gff_n_exons": row.get("n_exons", ""),
                    "gff_max_exon_bp": row.get("max_exon_bp", ""),
                    "gff_max_intron_bp": row.get("max_intron_bp", ""),
                    "gff_cds_bp": gff_cds,
                    "msa_n_exons": reference.get("n_exons", ""),
                    "msa_max_exon_bp": reference.get("max_exon_bp", ""),
                    "msa_max_intron_bp": reference.get("max_intron_bp", ""),
                    "msa_coding_bp": msa_coding,
                    "coding_bp_match": "yes"
                    if str(gff_cds) != "" and str(gff_cds) == str(msa_coding)
                    else "no",
                }
            )
        write_table(
            os.path.join(output_dir, "gff_vs_msa.tsv"),
            [
                "locus",
                "transcript_id",
                "gff_n_exons",
                "gff_max_exon_bp",
                "gff_max_intron_bp",
                "gff_cds_bp",
                "msa_n_exons",
                "msa_max_exon_bp",
                "msa_max_intron_bp",
                "msa_coding_bp",
                "coding_bp_match",
            ],
            comparison,
        )

    if "filter" in selected and enabled(cfg, "filter"):
        block = section(cfg, "filter")
        thresholds = resolve_thresholds(block, required=False)
        if thresholds is None:
            print(
                "No cutoffs yet. report.html will let you choose them. "
                "Download thresholds.json, set filter.thresholds_file, "
                "and rerun from the filter stage."
            )
        else:
            header_fixes = cfg.get("header_fixes") or []
            if os.path.isfile(species_map_path):
                species_by_code = read_species_map(species_map_path)
            elif block.get("species_map"):
                species_by_code = read_species_map(block["species_map"])
            else:
                species_by_code = {}
            run_names = block.get("runs")
            available = _runs(section(cfg, "msa_metrics"), "msa_metrics")
            if run_names:
                wanted = set(run_names)
                available = [run for run in available if run["name"] in wanted]
            for run in available:
                table = os.path.join(
                    output_dir, "msa_metrics__{0}.csv".format(run["name"])
                )
                _must_file(table)
                written, _kept = filter_run(
                    read_table(table),
                    run["alignment_dir"],
                    suffix,
                    os.path.join(output_dir, "filtered", run["name"]),
                    thresholds,
                    species_by_code,
                    header_fixes,
                )
                print("{0}: wrote {1} alignments".format(run["name"], len(written)))

    if "report" in selected and enabled(cfg, "report"):
        thresholds = resolve_thresholds(section(cfg, "filter"), required=False)
        run_rows = []
        for run in _runs(section(cfg, "msa_metrics"), "msa_metrics"):
            table = os.path.join(output_dir, "msa_metrics__{0}.csv".format(run["name"]))
            if os.path.isfile(table):
                run_rows.append((run["name"], read_table(table)))
        comparison_path = os.path.join(output_dir, "gff_vs_msa.tsv")
        comparison_rows = read_table(comparison_path) if os.path.isfile(comparison_path) else []
        qc_lines = _qc_lines(output_dir)
        write_report(output_dir, thresholds, run_rows, comparison_rows, qc_lines)
        print("Report: {0}".format(os.path.join(output_dir, "report.html")))
    return 0


def _runs(block, label):
    runs = block.get("runs")
    if not isinstance(runs, list) or not runs:
        raise ConfigError("{0}.runs must be a non-empty list.".format(label))
    for run in runs:
        if "name" not in run or "alignment_dir" not in run:
            raise ConfigError(
                "{0}.runs entries need name and alignment_dir.".format(label)
            )
    return runs


def _must_dir(path):
    if not os.path.isdir(path):
        raise ConfigError("Directory not found: {0}".format(path))


def _must_file(path):
    if not os.path.isfile(path):
        raise ConfigError("File not found: {0}".format(path))


def _qc_lines(output_dir):
    lines = []
    stats = os.path.join(output_dir, "transcriptome_stats.tsv")
    if os.path.isfile(stats):
        rows = read_table(stats)
        lines.append("Transcriptome assemblies: {0}".format(len(rows)))
        for row in rows:
            if int(row["n_transcripts"]) != int(row["n_genes"]):
                lines.append(
                    "  {0}: {1} transcripts, {2} Trinity genes".format(
                        row["file"], row["n_transcripts"], row["n_genes"]
                    )
                )
    busco = os.path.join(output_dir, "busco_summary.tsv")
    if os.path.isfile(busco):
        lines.append("BUSCO assemblies summarised: {0}".format(len(read_table(busco))))
        for row in read_table(busco):
            lines.append(
                "  {0}: C:{1}% S:{2}% D:{3}% F:{4}% M:{5}% n:{6}".format(
                    row["sample"],
                    row["Cpct"],
                    row["Spct"],
                    row["Dpct"],
                    row["Fpct"],
                    row["Mpct"],
                    row["n"],
                )
            )
    cds = os.path.join(output_dir, "reference_cds_check.tsv")
    if os.path.isfile(cds):
        rows = read_table(cds)
        n_match = sum(1 for row in rows if row["status"] == "MATCH")
        lines.append(
            "Reference CDS check: {0}/{1} loci MATCH after stripping N masks.".format(
                n_match, len(rows)
            )
        )
    isoform = os.path.join(output_dir, "longest_isoform", "summary.tsv")
    if os.path.isfile(isoform):
        for row in read_table(isoform):
            if row["n_transcripts"] != row["n_genes_kept"]:
                lines.append(
                    "Longest isoform {0}: {1} transcripts -> {2} genes".format(
                        row["file"], row["n_transcripts"], row["n_genes_kept"]
                    )
                )
    return lines
