# mmprefilter

`mmprefilter` takes MarkerMiner alignments and writes the subset to open in Geneious. It also records the transcriptome checks that were done before MarkerMiner.

MarkerMiner is not included. Run it first (Chamala et al., Applications in Plant Sciences, 2015). This repository starts from its `MAFFT_ADD_REF_ALIGN_FASTA` directory and stops at filtered FASTA files. It does not run Geneious and it does not design probes.

Every path, cutoff, and header correction is set in a JSON config file. The program does not have a second set of hidden defaults for those choices.

## Workflow

The assemblies are already made. MarkerMiner has already been run. The question this program answers is which of those loci are worth opening in Geneious.

Two checks happen before MarkerMiner, so a poor assembly is not mistaken for a poor locus. The rest happens on the MarkerMiner alignments.

### 1. Transcriptome statistics

- **Input.** Trinity assembly FASTA files.
- **Role.** Count transcripts, Trinity genes, total bases, N50, median length, mean length, and GC. These are the figures used to decide that an assembly is intact enough to mine.
- **Output.** `transcriptome_stats.tsv`, one row per assembly.

### 2. BUSCO summary

- **Input.** BUSCO `short_summary*.txt` files from transcriptome-mode runs.
- **Role.** Put complete, single-copy, duplicated, fragmented, and missing scores in one table. A low complete score means the assembly is a bad starting point, before any locus is judged.
- **Output.** `busco_summary.tsv`. BUSCO is not run here.

### 3. Longest isoform

- **Input.** The same Trinity FASTA files.
- **Role.** Keep the longest sequence of each Trinity gene (`TRINITY_DNx_c0_g1`). MarkerMiner can then be run twice: once on every isoform, once on one sequence per gene. The grouping is the same as Trinity's `get_longest_isoform_seq_per_trinity_gene.pl`.
- **Output.** `longest_isoform/*.fasta` and `summary.tsv`. This step prepares an input. It does not read MarkerMiner alignments.

### 4. Reference CDS check

- **Input.** MarkerMiner's intron-masked reference CDS, an external spliced CDS (TAIR10 representative models in this project), and the alignment filenames, which are the locus ids.
- **Role.** Delete the `N` runs from the masked CDS and compare what remains with the spliced CDS. A match confirms that the reference really is that gene, and records the transcript id (`TAIR_tx`) whose exons should be read from the GFF.
- **Output.** `reference_cds_check.tsv`, with `MATCH`, `MISMATCH`, or `MISSING` per locus.

### 5. GFF exon and intron lengths

- **Input.** A GFF3 file of `CDS` features, and the transcript id from step 4.
- **Role.** Measure the gene as it is annotated on the genome: exon lengths are CDS block lengths, intron lengths are the genomic gaps between those blocks. This is the structure you would see on a genome browser, and it is independent of how the alignment was built.
- **Output.** `gff_metrics.tsv`, one row per locus.

### 6. Exonic-column QC

- **Input.** Each `*.mafft.added.reference.align.fna` alignment. The reference is the last sequence.
- **Role.** On columns where that reference is `ACGT`, measure how gappy each ingroup sequence is and how often it matches the reference. This asks whether the alignment is usable, separately from how many exons the gene has.
- **Output.** `exon_qc__<run>.tsv`, one row per locus and ingroup sequence.

### 7. Alignment exon and intron lengths

- **Input.** The same alignments.
- **Role.** Read structure from the reference row itself. Gaps are removed, then `ACGT` runs are exons and `N` runs are introns. For each ingroup sequence, the length of a block is the width of that run minus the gaps in that sequence. These are the lengths present in the file that will be opened in Geneious, so this is the table the filter uses.
- **Output.** `msa_metrics__<run>.csv`, one row per locus and ingroup sequence: `n_exons`, `n_introns`, `max_exon_bp`, `max_intron_bp`, `avg_exon_bp`, `avg_intron_bp`.

### 8. GFF versus alignment

- **Input.** `gff_metrics.tsv` and the reference row of each alignment.
- **Role.** Place genomic CDS length and genomic intron length next to the `ACGT` and `N` lengths from the alignment. Coding lengths should agree if step 4 matched. Intron lengths are a different measurement: the GFF gap is a genomic distance, the alignment intron is an `N` run. The filter follows the alignment.
- **Output.** `gff_vs_msa.tsv`.

### 9. Filter

- **Input.** The step 7 tables, the original alignments, the transcriptome filenames, and the six cutoffs in `filter.thresholds`.
- **Role.** Keep a row only when every cutoff holds. A blank cell does not fail the row. A locus is written when at least one ingroup row passes, and only those rows are kept, plus the reference. `DAT` codes in the headers become species names taken from the assembly filenames. `header_fixes` then corrects known spelling (`echitum` to `echium` in this project).
- **Output.** `filtered/<run>/`, containing one alignment per kept locus, `namelist.txt`, `filtered_rows.csv`, and `summary.txt`.

### 10. Report

- **Input.** The metric tables and the cutoffs in the config.
- **Role.** Show how many rows and loci pass, and draw the six distributions with the cutoff marked. Changing a cutoff means editing the config and rerunning from the filter. The page does not store a new cutoff.
- **Output.** `report.txt` and `report.html`.

Steps 4 and 5 need an external annotation, so they apply to an Arabidopsis-reference MarkerMiner run. Steps 6, 7, 9, and 10 apply to every run listed in the config, including a Vitis reference.

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
