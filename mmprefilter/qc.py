"""Checks done on Trinity assemblies before MarkerMiner."""

import csv
import os
import re
import statistics

from .fastaio import read_fasta, write_fasta


TRINITY_GENE = re.compile(r"^(TRINITY_DN\d+_c\d+_g\d+)_i\d+$")
BUSCO_LINE = re.compile(
    r"C:(?P<Cpct>[\d.]+)%\[S:(?P<Spct>[\d.]+)%,D:(?P<Dpct>[\d.]+)%\],"
    r"F:(?P<Fpct>[\d.]+)%,M:(?P<Mpct>[\d.]+)%,n:(?P<n>\d+)"
)
BUSCO_FILE = re.compile(
    r"^#\s*Summarized benchmarking in BUSCO notation for file:?\s+(\S+)"
)
COUNT_LINE = {
    "C": re.compile(r"^\s*(\d+)\s+Complete BUSCOs \(C\)"),
    "S": re.compile(r"^\s*(\d+)\s+Complete and single-copy BUSCOs \(S\)"),
    "D": re.compile(r"^\s*(\d+)\s+Complete and duplicated BUSCOs \(D\)"),
    "F": re.compile(r"^\s*(\d+)\s+Fragmented BUSCOs \(F\)"),
    "M": re.compile(r"^\s*(\d+)\s+Missing BUSCOs \(M\)"),
}


def gene_id(seq_id):
    match = TRINITY_GENE.match(seq_id)
    if match:
        return match.group(1)
    return seq_id


def n50(lengths):
    if not lengths:
        return 0
    ordered = sorted(lengths, reverse=True)
    half = sum(ordered) / 2.0
    running = 0
    for length in ordered:
        running += length
        if running >= half:
            return length
    return ordered[-1]


def gc_fraction(sequences):
    gc = 0
    atgc = 0
    for seq in sequences:
        for base in seq.upper():
            if base in "GCgc":
                gc += 1
                atgc += 1
            elif base in "AT":
                atgc += 1
    if atgc == 0:
        return ""
    return round(gc / float(atgc), 4)


def transcriptome_stats(fasta_dir):
    rows = []
    names = sorted(
        name
        for name in os.listdir(fasta_dir)
        if name.lower().endswith((".fasta", ".fa", ".fna"))
    )
    for name in names:
        records = read_fasta(os.path.join(fasta_dir, name))
        lengths = [len(seq) for _, seq in records]
        genes = {gene_id(seq_id) for seq_id, _ in records}
        rows.append(
            {
                "file": name,
                "n_transcripts": len(records),
                "n_genes": len(genes),
                "total_bp": sum(lengths),
                "n50_bp": n50(lengths),
                "median_bp": int(statistics.median(lengths)) if lengths else 0,
                "mean_bp": round(statistics.mean(lengths), 1) if lengths else 0,
                "gc_fraction": gc_fraction([seq for _, seq in records]),
            }
        )
    return rows


def write_transcriptome_stats(fasta_dir, output_tsv):
    rows = transcriptome_stats(fasta_dir)
    _write_dicts(
        output_tsv,
        [
            "file",
            "n_transcripts",
            "n_genes",
            "total_bp",
            "n50_bp",
            "median_bp",
            "mean_bp",
            "gc_fraction",
        ],
        rows,
    )
    return rows


def longest_isoform_records(records):
    """Keep the longest sequence for each Trinity gene id.

    This is the same grouping rule as Trinity's
    get_longest_isoform_seq_per_trinity_gene.pl. On a length tie the first
    sequence in the file is kept.
    """
    best = {}
    order = []
    for seq_id, seq in records:
        gene = gene_id(seq_id)
        if gene not in best:
            order.append(gene)
            best[gene] = (seq_id, seq)
        elif len(seq) > len(best[gene][1]):
            best[gene] = (seq_id, seq)
    return [best[gene] for gene in order]


def write_longest_isoforms(fasta_dir, output_dir):
    if not os.path.isdir(output_dir):
        os.makedirs(output_dir)
    summary_path = os.path.join(output_dir, "summary.tsv")
    rows = []
    names = sorted(
        name
        for name in os.listdir(fasta_dir)
        if name.lower().endswith((".fasta", ".fa", ".fna"))
    )
    for name in names:
        records = read_fasta(os.path.join(fasta_dir, name))
        kept = longest_isoform_records(records)
        stem = name
        for suffix in (".fasta", ".fa", ".fna"):
            if stem.lower().endswith(suffix):
                stem = stem[: -len(suffix)]
                break
        out_name = stem + ".longest_isoform.fasta"
        write_fasta(os.path.join(output_dir, out_name), kept)
        rows.append(
            {
                "file": name,
                "output": out_name,
                "n_transcripts": len(records),
                "n_genes_kept": len(kept),
            }
        )
    _write_dicts(
        summary_path,
        ["file", "output", "n_transcripts", "n_genes_kept"],
        rows,
    )
    return rows


def sample_name_from_busco_input(input_path, run_dir):
    if not input_path:
        return os.path.basename(run_dir)
    base = os.path.basename(input_path)
    for suffix in (".Trinity.fasta", ".fasta", ".fa", ".fna"):
        if base.endswith(suffix):
            base = base[: -len(suffix)]
            break
    parts = base.split("_")
    if len(parts) >= 2:
        return "_".join(parts[:2])
    return base or os.path.basename(run_dir)


def parse_busco_short_summary(path, run_dir):
    with open(path, encoding="utf-8", errors="replace") as handle:
        text = handle.read()
    input_path = ""
    for line in text.splitlines():
        match = BUSCO_FILE.match(line.strip())
        if match:
            input_path = match.group(1)
            break
    scores = None
    for line in text.splitlines():
        match = BUSCO_LINE.search(line)
        if match:
            scores = match.groupdict()
            break
    if not scores:
        return None
    counts = {}
    for line in text.splitlines():
        for key, pattern in COUNT_LINE.items():
            match = pattern.match(line)
            if match:
                counts[key] = match.group(1)
    return {
        "sample": sample_name_from_busco_input(input_path, run_dir),
        "input_file": input_path,
        "results_dir": run_dir,
        "Cpct": scores["Cpct"],
        "Spct": scores["Spct"],
        "Dpct": scores["Dpct"],
        "Fpct": scores["Fpct"],
        "Mpct": scores["Mpct"],
        "C": counts.get("C", ""),
        "S": counts.get("S", ""),
        "D": counts.get("D", ""),
        "F": counts.get("F", ""),
        "M": counts.get("M", ""),
        "n": scores["n"],
    }


def write_busco_summary(busco_root, output_tsv):
    rows = []
    for dirpath, _, filenames in os.walk(busco_root):
        for name in sorted(filenames):
            if not name.startswith("short_summary"):
                continue
            if not name.endswith(".txt"):
                continue
            parsed = parse_busco_short_summary(os.path.join(dirpath, name), dirpath)
            if parsed:
                rows.append(parsed)
    rows.sort(key=lambda row: row["sample"])
    _write_dicts(
        output_tsv,
        [
            "sample",
            "input_file",
            "results_dir",
            "Cpct",
            "Spct",
            "Dpct",
            "Fpct",
            "Mpct",
            "C",
            "S",
            "D",
            "F",
            "M",
            "n",
        ],
        rows,
    )
    return rows


def parse_species_from_filename(path):
    """DAT8-lithospermum_erythrorhizon_localscratch_trinity.Trinity.fasta -> DAT8, species."""
    basename = os.path.basename(path.strip())
    if "-" not in basename:
        return None, None
    code, rest = basename.split("-", 1)
    if "_localscratch" in rest:
        species = rest.split("_localscratch")[0]
    else:
        species = re.sub(r"\.(?:fasta|fa|fna).*$", "", rest, flags=re.I)
        species = re.sub(r"_trinity.*$", "", species, flags=re.I).rstrip("_")
    return code, species


def read_path_list(path):
    base = os.path.dirname(os.path.abspath(path))
    paths = []
    with open(path) as handle:
        for line in handle:
            line = line.strip()
            if not line or line.startswith("#"):
                continue
            if not os.path.isabs(line):
                line = os.path.normpath(os.path.join(base, line))
            paths.append(line)
    return paths


def species_map_from_list(list_path):
    mapping = {}
    for path in read_path_list(list_path):
        code, species = parse_species_from_filename(path)
        if code and species:
            mapping[code] = species
    return mapping


def read_species_map(path):
    mapping = {}
    with open(path, newline="") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        for row in reader:
            code = (row.get("dat_code") or "").strip()
            species = (row.get("species") or "").strip()
            if code and species:
                mapping[code] = species
    return mapping


def write_species_map(path, mapping):
    directory = os.path.dirname(path)
    if directory and not os.path.isdir(directory):
        os.makedirs(directory)
    with open(path, "w") as handle:
        handle.write("dat_code\tspecies\n")
        for code in sorted(mapping):
            handle.write("{0}\t{1}\n".format(code, mapping[code]))


def annotate_seq_id(seq_id, species_by_code, header_fixes):
    if "_" not in seq_id:
        return apply_header_fixes(seq_id, header_fixes)
    code, rest = seq_id.split("_", 1)
    if code in species_by_code:
        seq_id = "{0}_{1}_{2}".format(code, species_by_code[code], rest)
    return apply_header_fixes(seq_id, header_fixes)


def apply_header_fixes(text, header_fixes):
    for wrong, right in header_fixes or []:
        text = text.replace(wrong, right)
    return text


def _write_dicts(path, fieldnames, rows):
    directory = os.path.dirname(path)
    if directory and not os.path.isdir(directory):
        os.makedirs(directory)
    delimiter = "\t" if path.endswith(".tsv") else ","
    with open(path, "w", newline="") as handle:
        writer = csv.DictWriter(
            handle, fieldnames=fieldnames, extrasaction="ignore", delimiter=delimiter
        )
        writer.writeheader()
        for row in rows:
            writer.writerow(row)
