"""Exon and intron structure from a reference CDS, a GFF, and the alignment itself.

MarkerMiner puts the reference sequence last in each
*.mafft.added.reference.align.fna file. In that sequence, ACGT runs are exons
and N runs are intron placeholders. Those placeholders are not genomic intron
lengths. The GFF stage records the genomic lengths so the two can be compared.
The filter uses the alignment.
"""

import csv
import os
import re
from collections import defaultdict

from .fastaio import read_fasta, read_fasta_full_headers, write_fasta
from .qc import annotate_seq_id


ALIGNMENT_TOKEN = re.compile(r"N+|[ACGT]+")


def strip_intron_n(seq):
    return re.sub(r"N+", "", seq.upper())


def parse_spliced_header(header):
    """TAIR-style header -> (locus_id, transcript_id)."""
    token = header.split("|", 1)[0].strip().split()[0]
    if "." in token:
        locus, version = token.rsplit(".", 1)
        if version.isdigit():
            return locus, token
    return token, token


def loci_in_alignment_dir(alignment_dir, suffix):
    stems = []
    for name in sorted(os.listdir(alignment_dir)):
        if name.endswith(suffix):
            stems.append(name[: -len(suffix)])
    return stems


def msa_segments_from_ref(ref_seq):
    """Ordered (kind, [column indexes]) with kind 'coding' or 'intron'.

    Gaps are removed before splitting. ACGT runs are coding. N runs are the
    intron mask. Same tokenisation as the original post-MarkerMiner scripts.
    """
    chars = []
    columns = []
    for index, base in enumerate(ref_seq):
        upper = base.upper()
        if upper in "-.":
            continue
        if upper not in "ACGTN":
            continue
        chars.append(upper)
        columns.append(index)
    if not chars:
        return []
    segments = []
    cursor = 0
    for token in ALIGNMENT_TOKEN.findall("".join(chars)):
        length = len(token)
        token_columns = columns[cursor : cursor + length]
        cursor += length
        kind = "intron" if token[0] == "N" else "coding"
        segments.append((kind, token_columns))
    return segments


def _gap_count(query, columns):
    return sum(1 for index in columns if query[index].upper() in "-.")


def lengths_for_query(query, segments):
    exon_lengths = []
    intron_lengths = []
    for kind, columns in segments:
        if not columns:
            continue
        observed = len(columns) - _gap_count(query, columns)
        if kind == "coding":
            exon_lengths.append(observed)
        else:
            intron_lengths.append(observed)
    return exon_lengths, intron_lengths


def summarise_lengths(exon_lengths, intron_lengths):
    def pack(lengths):
        if not lengths:
            return 0, "", ""
        average = round(sum(lengths) / float(len(lengths)), 2)
        return len(lengths), max(lengths), average

    n_exons, max_exon, avg_exon = pack(exon_lengths)
    n_introns, max_intron, avg_intron = pack(intron_lengths)
    return {
        "n_exons": n_exons,
        "n_introns": n_introns,
        "max_exon_bp": max_exon,
        "max_intron_bp": max_intron,
        "avg_exon_bp": avg_exon,
        "avg_intron_bp": avg_intron,
    }


def _number(value, as_int):
    if value is None:
        return None
    text = str(value).strip()
    if text == "":
        return None
    try:
        number = int(text) if as_int else float(text)
    except ValueError:
        return None
    return number


def row_matches(row, thresholds):
    """AND of the six config cutoffs. A blank cell does not fail the row."""
    checks = [
        ("n_exons", True, "n_exons_gt", True),
        ("n_introns", True, "n_introns_lt", False),
        ("max_exon_bp", True, "max_exon_bp_gt", True),
        ("max_intron_bp", True, "max_intron_bp_lt", False),
        ("avg_exon_bp", False, "avg_exon_bp_gt", True),
        ("avg_intron_bp", False, "avg_intron_bp_lt", False),
    ]
    for column, as_int, key, greater in checks:
        value = _number(row.get(column), as_int)
        if value is None:
            continue
        cutoff = thresholds[key]
        if greater and not (value > cutoff):
            return False
        if not greater and not (value < cutoff):
            return False
    return True


def iter_alignments(alignment_dir, suffix):
    for name in sorted(os.listdir(alignment_dir)):
        if name.endswith(suffix):
            yield name[: -len(suffix)], os.path.join(alignment_dir, name)


def alignment_records(path):
    records = read_fasta(path)
    if len(records) < 2:
        return None
    lengths = {len(seq) for _, seq in records}
    if len(lengths) != 1:
        return None
    return records


def msa_metric_rows(alignment_dir, suffix):
    rows = []
    skipped = []
    for locus, path in iter_alignments(alignment_dir, suffix):
        records = alignment_records(path)
        if not records:
            skipped.append(locus)
            continue
        ref_id, ref_seq = records[-1]
        segments = msa_segments_from_ref(ref_seq)
        if not any(kind == "coding" for kind, _ in segments):
            skipped.append(locus)
            continue
        for seq_id, seq in records:
            if seq_id == ref_id:
                continue
            exons, introns = lengths_for_query(seq, segments)
            if not exons:
                continue
            summary = summarise_lengths(exons, introns)
            summary["locus"] = locus
            summary["ingroup_id"] = seq_id
            rows.append(summary)
    return rows, skipped


def reference_block_summary(ref_seq):
    segments = msa_segments_from_ref(ref_seq)
    exons = [len(columns) for kind, columns in segments if kind == "coding"]
    introns = [len(columns) for kind, columns in segments if kind == "intron"]
    summary = summarise_lengths(exons, introns)
    summary["coding_bp"] = sum(exons)
    return summary


def exon_qc_rows(alignment_dir, suffix):
    """Gap fraction and identity on reference ACGT columns, per ingroup row."""
    rows = []
    for locus, path in iter_alignments(alignment_dir, suffix):
        records = alignment_records(path)
        if not records:
            continue
        ref_id, ref_seq = records[-1]
        coding_columns = [
            index
            for index, base in enumerate(ref_seq)
            if base.upper() in "ACGT"
        ]
        n_exonic = len(coding_columns)
        if n_exonic == 0:
            continue
        for seq_id, seq in records:
            if seq_id == ref_id:
                continue
            gaps = 0
            matches = 0
            comparable = 0
            for index in coding_columns:
                query = seq[index].upper()
                reference = ref_seq[index].upper()
                if query in "-.":
                    gaps += 1
                if reference in "ACGT" and query in "ACGT":
                    comparable += 1
                    if reference == query:
                        matches += 1
            identity = ""
            if comparable:
                identity = round(matches / float(comparable), 6)
            rows.append(
                {
                    "locus": locus,
                    "ingroup_id": seq_id,
                    "gap_frac_exonic": round(gaps / float(n_exonic), 6),
                    "ident_vs_ref_exonic": identity,
                    "n_both_ACGT": comparable,
                    "n_exonic_cols": n_exonic,
                }
            )
    return rows


def compare_masked_to_spliced(masked_fasta, spliced_fasta, alignment_dir, suffix):
    masked = {}
    for header, seq in read_fasta_full_headers(masked_fasta):
        masked[header.split()[0]] = seq.upper().replace(" ", "")
    spliced = {}
    for header, seq in read_fasta_full_headers(spliced_fasta):
        locus, transcript = parse_spliced_header(header)
        spliced[locus] = (transcript, seq.upper().replace(" ", ""))
    rows = []
    for locus in loci_in_alignment_dir(alignment_dir, suffix):
        if locus not in masked or locus not in spliced:
            rows.append(
                {
                    "locus": locus,
                    "status": "MISSING",
                    "transcript_id": spliced.get(locus, ("",))[0] if locus in spliced else "",
                    "masked_coding_bp": len(strip_intron_n(masked[locus])) if locus in masked else "",
                    "spliced_bp": len(spliced[locus][1]) if locus in spliced else "",
                }
            )
            continue
        transcript, spliced_seq = spliced[locus]
        coding = strip_intron_n(masked[locus])
        status = "MATCH" if coding == spliced_seq else "MISMATCH"
        rows.append(
            {
                "locus": locus,
                "status": status,
                "transcript_id": transcript,
                "masked_coding_bp": len(coding),
                "spliced_bp": len(spliced_seq),
            }
        )
    return rows


def parse_gff_attributes(field):
    attributes = {}
    for part in field.strip().split(";"):
        part = part.strip()
        if not part or "=" not in part:
            continue
        key, value = part.split("=", 1)
        attributes[key.strip()] = value.strip()
    return attributes


def cds_blocks_by_transcript(gff_path):
    raw = defaultdict(list)
    with open(gff_path) as handle:
        for line in handle:
            if not line or line[0] == "#":
                continue
            parts = line.rstrip("\n").split("\t")
            if len(parts) < 9 or parts[2] != "CDS":
                continue
            attributes = parse_gff_attributes(parts[8])
            parent = attributes.get("Parent", "").split(",")[0].strip()
            if not parent:
                continue
            raw[parent].append(
                (parts[0], int(parts[3]), int(parts[4]), parts[6])
            )
    index = {}
    for transcript, rows in raw.items():
        strand = rows[0][3]
        if strand == "+":
            rows.sort(key=lambda item: item[1])
        else:
            rows.sort(key=lambda item: item[1], reverse=True)
        blocks = []
        for chrom, start, end, strand_value in rows:
            blocks.append(
                {
                    "chr": chrom,
                    "start": start,
                    "end": end,
                    "length": end - start + 1,
                    "strand": strand_value,
                }
            )
        index[transcript] = blocks
    return index


def genomic_intron_lengths(blocks):
    if len(blocks) < 2:
        return []
    strand = blocks[0]["strand"]
    lengths = []
    for upstream, downstream in zip(blocks, blocks[1:]):
        if upstream["chr"] != downstream["chr"]:
            continue
        if strand == "+":
            length = downstream["start"] - upstream["end"] - 1
        else:
            length = upstream["start"] - downstream["end"] - 1
        if length >= 0:
            lengths.append(length)
    return lengths


def gff_metric_rows(gff_path, transcript_by_locus, loci):
    index = cds_blocks_by_transcript(gff_path)
    rows = []
    for locus in loci:
        transcript = transcript_by_locus.get(locus, "")
        blocks = index.get(transcript) if transcript else None
        if not blocks:
            rows.append(
                {
                    "locus": locus,
                    "transcript_id": transcript,
                    "n_exons": "",
                    "n_introns": "",
                    "max_exon_bp": "",
                    "max_intron_bp": "",
                    "avg_exon_bp": "",
                    "avg_intron_bp": "",
                    "cds_bp": "",
                }
            )
            continue
        exon_lengths = [block["length"] for block in blocks]
        intron_lengths = genomic_intron_lengths(blocks)
        summary = summarise_lengths(exon_lengths, intron_lengths)
        summary["locus"] = locus
        summary["transcript_id"] = transcript
        summary["cds_bp"] = sum(exon_lengths)
        summary["n_introns"] = max(0, len(exon_lengths) - 1)
        rows.append(summary)
    return rows


def load_transcript_map(path):
    mapping = {}
    with open(path, newline="") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        for row in reader:
            if row.get("status") == "MATCH" and row.get("transcript_id"):
                mapping[row["locus"]] = row["transcript_id"]
    return mapping


def write_table(path, fieldnames, rows):
    directory = os.path.dirname(path)
    if directory and not os.path.isdir(directory):
        os.makedirs(directory)
    with open(path, "w", newline="") as handle:
        writer = csv.DictWriter(
            handle, fieldnames=fieldnames, extrasaction="ignore", delimiter="\t"
            if path.endswith(".tsv")
            else ","
        )
        writer.writeheader()
        for row in rows:
            writer.writerow(row)


MSA_FIELDS = [
    "locus",
    "ingroup_id",
    "n_exons",
    "n_introns",
    "max_exon_bp",
    "max_intron_bp",
    "avg_exon_bp",
    "avg_intron_bp",
]


def read_table(path):
    delimiter = "\t" if path.endswith(".tsv") else ","
    with open(path, newline="") as handle:
        return list(csv.DictReader(handle, delimiter=delimiter))


def filter_run(
    rows,
    alignment_dir,
    suffix,
    output_dir,
    thresholds,
    species_by_code,
    header_fixes,
):
    if os.path.isdir(output_dir):
        for name in os.listdir(output_dir):
            os.remove(os.path.join(output_dir, name))
    else:
        os.makedirs(output_dir)
    by_locus = defaultdict(list)
    for row in rows:
        if row_matches(row, thresholds):
            by_locus[row["locus"]].append(row)
    written = []
    kept_rows = []
    for locus in sorted(by_locus):
        keep_ids = {row["ingroup_id"] for row in by_locus[locus]}
        source = os.path.join(alignment_dir, locus + suffix)
        if not os.path.isfile(source):
            continue
        records = read_fasta(source)
        if not records:
            continue
        reference = records[-1]
        trimmed = [item for item in records[:-1] if item[0] in keep_ids]
        trimmed.append(reference)
        if len(trimmed) < 2:
            continue
        renamed = [
            (annotate_seq_id(seq_id, species_by_code, header_fixes), seq)
            for seq_id, seq in trimmed
        ]
        destination = os.path.join(output_dir, locus + suffix)
        write_fasta(destination, renamed)
        written.append(locus)
        kept_rows.extend(by_locus[locus])
    with open(os.path.join(output_dir, "namelist.txt"), "w") as handle:
        for locus in written:
            handle.write(locus + "\n")
    write_table(os.path.join(output_dir, "filtered_rows.csv"), MSA_FIELDS, kept_rows)
    with open(os.path.join(output_dir, "summary.txt"), "w") as handle:
        handle.write("loci_written: {0}\n".format(len(written)))
        handle.write("rows_kept: {0}\n".format(len(kept_rows)))
        handle.write("thresholds: {0}\n".format(thresholds))
    return written, kept_rows
