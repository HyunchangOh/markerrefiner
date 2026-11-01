"""FASTA readers and writers. Sequence ids are the first header token."""

import os


def read_fasta(path):
    """Return [(seq_id, sequence), ...] in file order."""
    records = []
    header = None
    chunks = []
    with open(path) as handle:
        for line in handle:
            line = line.strip()
            if not line:
                continue
            if line.startswith(">"):
                if header is not None:
                    records.append((header.split()[0], "".join(chunks)))
                header = line[1:]
                chunks = []
            else:
                chunks.append(line)
        if header is not None:
            records.append((header.split()[0], "".join(chunks)))
    return records


def read_fasta_full_headers(path):
    """Return [(full header without '>', sequence), ...]."""
    records = []
    header = None
    chunks = []
    with open(path) as handle:
        for line in handle:
            line = line.strip()
            if not line:
                continue
            if line.startswith(">"):
                if header is not None:
                    records.append((header, "".join(chunks)))
                header = line[1:]
                chunks = []
            else:
                chunks.append(line)
        if header is not None:
            records.append((header, "".join(chunks)))
    return records


def write_fasta(path, records, width=80):
    directory = os.path.dirname(path)
    if directory and not os.path.isdir(directory):
        os.makedirs(directory)
    with open(path, "w") as handle:
        for seq_id, seq in records:
            handle.write(">{0}\n".format(seq_id))
            for start in range(0, len(seq), width):
                handle.write(seq[start : start + width] + "\n")
