# mmprefilter

`mmprefilter` takes MarkerMiner alignments and writes the subset to open in Geneious. It also records the transcriptome checks that were done before MarkerMiner.

MarkerMiner is not included. Run it first (Chamala et al., Applications in Plant Sciences, 2015). This repository starts from its `MAFFT_ADD_REF_ALIGN_FASTA` directory and stops at filtered FASTA files. It does not run Geneious and it does not design probes.

Every path, cutoff, and header correction is set in a JSON config file. The program does not have a second set of hidden defaults for those choices.

## What it does

Before MarkerMiner:

1. **Transcriptome statistics.** Sequence count, Trinity gene count, total bases, N50, median, mean, and GC, from the assembly FASTA. This replaces the counts we used from `TrinityStats.pl` for the decision of whether an assembly was usable.
2. **BUSCO summary.** Reads BUSCO `short_summary*.txt` files (transcriptome mode) and writes one table. It does not run BUSCO.
3. **Longest isoform.** Keeps the longest sequence per Trinity gene (`TRINITY_DNx_c0_g1`). This is the same grouping as Trinity's `get_longest_isoform_seq_per_trinity_gene.pl`, so MarkerMiner can be run a second time on one transcript per gene.

After MarkerMiner:

4. **Reference CDS check.** Strips `N` runs from an intron-masked reference CDS and compares what remains to an external spliced CDS (TAIR10 representative CDS in this project). A match says which transcript the reference corresponds to.
5. **GFF exon and intron lengths.** For that transcript, CDS blocks in a GFF give genomic exon lengths and the genomic gaps between them.
6. **Exonic-column QC.** On columns where the alignment reference is `ACGT`, each ingroup sequence gets a gap fraction and an identity to the reference.
7. **Alignment exon and intron lengths.** The last sequence in each MarkerMiner alignment is the reference. After gaps are removed, `ACGT` runs are exons and `N` runs are introns. For every other sequence, the length of a block is the number of those columns minus that sequence's gaps. This is the table the filter uses.
8. **GFF versus alignment.** Coding length and longest intron from the GFF, next to the same numbers from the reference row. Where they differ, the alignment is what Geneious will show, so the filter follows the alignment.
9. **Filter.** A row is kept only if it passes every cutoff in `filter.thresholds`. A blank cell does not fail the row. A locus is written if at least one ingroup row passes; only those rows, plus the reference, are kept. `DAT` codes in the headers are replaced with species names from the transcriptome filenames. `header_fixes` is applied after that (in this project, `echitum` to `echium`).
10. **Report.** `report.txt` and `report.html`. The histograms mark the config cutoffs. They do not edit the cutoffs.

## Requirements

Python 3.8 or newer. No packages to install.

BUSCO and Trinity are needed only if you still have to produce the assemblies, the BUSCO runs, or the MarkerMiner alignments. This program reads those outputs.

## Quick start

From this directory:

```bash
python3 -m mmprefilter run --config example/config.json
```

Open `example/output/report.html`.

The example config uses the cutoffs from the Boraginaceae screen:

| Cutoff | Value |
|---|---|
| `n_exons` | > 0 |
| `n_introns` | < 5 |
| `max_exon_bp` | > 220 |
| `max_intron_bp` | < 121 |
| `avg_exon_bp` | > 152.51259 |
| `avg_intron_bp` | < 102.46 |

On the example alignments that keeps `AT1G01180` and `AT1G01280`, and drops `AT1G01880` (6 introns) and `VV00G01720` (8 introns, longest exon 131 bp). `DAT5_TRINITY_...` is renamed to `DAT5_echium_plantagineum_TRINITY_...`.

To try other cutoffs, edit `filter.thresholds` and rerun from the filter onward. The metric tables are left in place:

```bash
python3 -m mmprefilter run --config example/config.json --from filter
```

`--from` and `--to` choose which stages to execute. They are not analysis settings. Analysis settings are only in the config.

## Running it on a full MarkerMiner output

Copy `example/config.json`. Point `alignment_dir` entries at each `MAFFT_ADD_REF_ALIGN_FASTA` directory, and point the transcriptome, BUSCO, CDS, and GFF entries at your files. Set `enabled` to `false` for any stage you do not have inputs for. The reference CDS check and the GFF stage apply to the Arabidopsis-reference run. The alignment metrics and the filter apply to every run listed under `msa_metrics.runs`.

`config.study.json` is the same cutoffs and the four runs from this project (Arabidopsis or Vitis reference, all isoforms or longest isoform), with the paths left for you to fill in. On the full alignments those cutoffs kept:

| Run | Loci written |
|---|---|
| Arabidopsis, all isoforms | 606 |
| Arabidopsis, longest isoform | 540 |
| Vitis, all isoforms | 409 |
| Vitis, longest isoform | 369 |

A Vitis-to-Arabidopsis BLAST was used while choosing cutoffs and was not applied when those subsets were written. It is not part of this program.

## Config

Paths in the config are relative to the config file unless they are absolute. Paths inside `transcriptome_list` are relative to that list file.

| Key | Role |
|---|---|
| `output_dir` | Where tables, filtered alignments, and the report are written |
| `alignment_suffix` | MarkerMiner alignment filename ending |
| `reference_row` | Must be `last` |
| `header_fixes` | Pairs of strings to replace in sequence names, in order |
| `transcriptome_qc.fasta_dir` | Trinity FASTA files |
| `busco_summary.busco_root` | Directory tree of BUSCO `short_summary*.txt` files |
| `longest_isoform.fasta_dir` | FASTA files to reduce to one sequence per Trinity gene |
| `species_map.transcriptome_list` | One assembly path per line. The filename must look like `DAT8-lithospermum_erythrorhizon...fasta` |
| `reference_cds_check` | Intron-masked CDS, spliced CDS, and the alignment directory whose filenames are the locus ids |
| `gff_metrics` | GFF3 with `CDS` features. Transcript ids come from the CDS check (`MATCH` rows) unless `transcript_map` is set |
| `exon_qc.runs` / `msa_metrics.runs` | Each entry has `name` and `alignment_dir` |
| `filter.runs` | Which metric tables to filter. Omit it to filter every `msa_metrics` run |
| `filter.thresholds` | The six cutoffs. All six keys are required |
| `report.enabled` | Write `report.txt` and `report.html` |

`filter.thresholds` uses a strict inequality (`>` or `<`, not `>=`). A blank `max_intron_bp` or `avg_intron_bp`, as on a single-exon gene, is ignored for that cutoff. `n_introns` of `0` is a real value and is tested.

## Outputs

| File | Contents |
|---|---|
| `transcriptome_stats.tsv` | One row per assembly FASTA |
| `busco_summary.tsv` | One row per BUSCO short summary |
| `longest_isoform/` | One FASTA per input, plus `summary.tsv` |
| `species_map.tsv` | `dat_code` and `species` |
| `reference_cds_check.tsv` | `MATCH`, `MISMATCH`, or `MISSING` per locus |
| `gff_metrics.tsv` | Genomic exon and intron lengths |
| `exon_qc__<run>.tsv` | Exonic gap fraction and identity |
| `msa_metrics__<run>.csv` | One row per locus and ingroup sequence |
| `gff_vs_msa.tsv` | Genomic lengths beside alignment-mask lengths |
| `filtered/<run>/` | Alignments, `namelist.txt`, `filtered_rows.csv`, `summary.txt` |
| `report.html` | Counts and histograms for the current cutoffs |

## Example data

`example/data/PROVENANCE.txt` lists what was copied. The alignments, the three reference CDS records, the GFF CDS lines, and the BUSCO scores are from the Boraginaceae project. The transcriptome files are the first record of each assembly, plus three isoforms of one *Echium vulgare* gene so the longest-isoform step removes something.

## Tests

```bash
python3 -m unittest discover -s tests
```

## Citation

If you use this filter, cite MarkerMiner for the alignments and cite this repository for the subset. Replace the author line in `CITATION.cff` before submission.

BUSCO summaries follow Manni et al., Molecular Biology and Evolution, 2021. Longest-isoform grouping follows Trinity (Grabherr et al., Nature Biotechnology, 2011).

## Licence

MIT. See `LICENSE`.
