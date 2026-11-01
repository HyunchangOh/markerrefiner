# mmprefilter

`mmprefilter` takes MarkerMiner alignments and writes the subset to open in Geneious. It also records the transcriptome checks that were done before MarkerMiner.

MarkerMiner is not included. Run it first (Chamala et al., Applications in Plant Sciences, 2015). This repository starts from its `MAFFT_ADD_REF_ALIGN_FASTA` directory and stops at filtered FASTA files. It does not run Geneious and it does not design probes.

Every path, cutoff, and header correction is set in a JSON config file. The program does not have a second set of hidden defaults for those choices.

## Workflow

MarkerMiner returns a locus for every reference gene it can find. Most of those loci are a poor use of manual time in Geneious. The usual reason to throw one out, once it is open, is that the gene is the wrong shape for a marker: no exon long enough to use, too many intron breaks, or a long empty stretch, or that the ingroup sequence is barely there. This program measures those problems first and writes only the loci that pass.

The first three stages look at the assemblies, before MarkerMiner. A broken assembly makes every locus from that species look bad. The later stages look at the alignments. Only stage 9 deletes loci. The earlier stages tell you whether the numbers in that deletion are meaningful.

### 1. Transcriptome statistics

- **Purpose.** Decide whether each assembly is complete enough to mine. A tiny or fragmented assembly will miss loci for reasons that have nothing to do with the gene.
- **Input.** Trinity assembly FASTA files.
- **Role.** Count transcripts, Trinity genes, total bases, N50, median length, mean length, and GC.
- **Red flags.** An assembly that is far smaller than the others, an N50 of only a few hundred bases, or a GC fraction that sits away from the rest of the species. That species will contribute gaps later, and those gaps are an assembly failure.
- **Okay to continue.** The assemblies you intend to mine sit in the same range as each other. This stage does not delete files. You set an outlier aside before running MarkerMiner on it.
- **Output.** `transcriptome_stats.tsv`, one row per assembly.

### 2. BUSCO summary

- **Purpose.** Ask whether conserved plant genes are actually present. Transcript counts can look fine when most genes are still missing.
- **Input.** BUSCO `short_summary*.txt` files from transcriptome-mode runs.
- **Role.** Put complete, single-copy, duplicated, fragmented, and missing scores in one table.
- **Red flags.** A low complete percentage and a high missing percentage. In this project *Echium vulgare* and *Mertensia paniculata* were near 40% complete, while *Lithospermum erythrorhizon* was about 95% complete. The low assemblies still went into MarkerMiner. The warning is that those species will be absent from many loci, so a missing sequence there is not evidence that the locus is bad.
- **Okay to continue.** You can name the weak assemblies before you interpret gaps. This program does not apply a BUSCO cutoff.
- **Output.** `busco_summary.tsv`. BUSCO is not run here.

### 3. Longest isoform

- **Purpose.** Give MarkerMiner a second input in which each Trinity gene appears once. Several isoforms of one gene can hand MarkerMiner a short fragment, or a different splice form, of the same locus.
- **Input.** The same Trinity FASTA files.
- **Role.** Keep the longest sequence of each Trinity gene (`TRINITY_DNx_c0_g1`). The grouping is the same as Trinity's `get_longest_isoform_seq_per_trinity_gene.pl`.
- **Red flags.** A gene with several isoforms of very different lengths. The short ones are the sequences most likely to look like a truncated exon later.
- **Okay to continue.** The output has one sequence per gene, and `n_genes_kept` is lower than `n_transcripts` wherever a gene had more than one isoform. MarkerMiner is then run on this FASTA as well as on the full assembly. This stage does not read MarkerMiner alignments.
- **Output.** `longest_isoform/*.fasta` and `summary.tsv`.

### 4. Reference CDS check

- **Purpose.** Confirm that the reference sequence inside MarkerMiner is the gene you think it is, before any exon is measured from an annotation.
- **Input.** MarkerMiner's intron-masked reference CDS, an external spliced CDS (TAIR10 representative models in this project), and the alignment filenames, which are the locus ids.
- **Role.** Delete the `N` runs from the masked CDS and compare what remains with the spliced CDS. A match records the transcript id whose exons should be read from the GFF.
- **Red flags.** `MISMATCH`: the masked reference is not that spliced CDS, so exon boundaries taken from the annotation would describe a different sequence. `MISSING`: the locus id is absent from one of the two files, so there is nothing to check.
- **Okay to continue.** Status is `MATCH`, and the row carries a transcript id. In this project nearly every Arabidopsis locus matched. A mismatch stays out of the GFF measurement.
- **Output.** `reference_cds_check.tsv`.

### 5. GFF exon and intron lengths

- **Purpose.** See the gene as it is annotated on the genome, which is the picture you would use if you aligned the locus to a reference in Geneious and looked at exon and intron sizes there.
- **Input.** A GFF3 file of `CDS` features, and the transcript id from step 4.
- **Role.** Exon length is the length of a CDS block. Intron length is the genomic gap between adjacent blocks.
- **Red flags.** No exon long enough to be a marker (a few hundred bases was the range we had in mind), a long run of tiny exons, a large intron count, or one very long intron. Those are the structures that look awkward in a genome browser and get discarded by eye.
- **Okay to continue.** The locus has a CDS model and the lengths are filled in. This stage records the annotation. It does not drop loci. The drop uses the alignment, because that is the file Geneious will open.
- **Output.** `gff_metrics.tsv`, one row per locus.

### 6. Exonic-column QC

- **Purpose.** Separate a bad alignment from a gene that merely has a difficult structure. A locus can have a normal exon and still be useless if the ingroup sequence is mostly empty, or if it does not match the reference.
- **Input.** Each `*.mafft.added.reference.align.fna` alignment. The reference is the last sequence.
- **Role.** On columns where that reference is `ACGT`, measure the gap fraction of each ingroup sequence and its identity to the reference.
- **Red flags.** A large share of coding columns are gaps, so that species is not really in the locus. Identity to the reference is far below the other sequences in the same alignment, which suggests the wrong transcript was recruited.
- **Okay to continue.** Coding columns are mostly occupied, and identity is in line with the other ingroup sequences at that locus. This stage does not apply a numeric cutoff. It is the table you check when a locus that passed the structure filter still looks wrong.
- **Output.** `exon_qc__<run>.tsv`, one row per locus and ingroup sequence.

### 7. Alignment exon and intron lengths

- **Purpose.** Score the structure that is actually in the alignment. Geneious shows this file, not the genomic coordinates from the GFF.
- **Input.** The same alignments.
- **Role.** Gaps are removed from the reference row. `ACGT` runs are exons and `N` runs are introns. For each ingroup sequence, the length of a block is the width of that run minus the gaps in that sequence.
- **Red flags.** The same biological problems as in stage 5, now measured on the alignment: no real exon left after gaps, many intron breaks, a long `N` run, or exons that are short on average even when one of them is long.
- **Okay to continue.** Every ingroup sequence you might keep has a row, with exon counts filled in. A single-exon gene has `n_introns` of 0 and a blank intron length. Blank means "there is no intron to judge", and the intron cutoffs in the next stage ignore that cell.
- **Output.** `msa_metrics__<run>.csv`, one row per locus and ingroup sequence.

### 8. GFF versus alignment

- **Purpose.** Check that the two ways of measuring the gene describe the same CDS before you trust either one. If they disagree, filtering on the genomic intron length would throw out loci for a gap that is not the gap in the alignment.
- **Input.** `gff_metrics.tsv` and the reference row of each alignment.
- **Role.** Place genomic CDS length and genomic intron length next to the `ACGT` and `N` lengths from the alignment.
- **Red flags.** `coding_bp_match` is `no`: the reference row is not the annotated CDS, so do not use that locus's GFF exons. A large difference in longest-intron length means a genomic cutoff and an alignment cutoff are not the same decision.
- **Okay to continue.** Coding length agrees. You then filter on the alignment columns. Intron lengths may still differ, and that difference is the reason the next stage does not use the GFF.
- **Output.** `gff_vs_msa.tsv`.

### 9. Filter

- **Purpose.** Remove the loci that would be discarded in Geneious for structure, and write the remaining alignments with readable species names.
- **Input.** The stage 7 tables, the original alignments, the transcriptome filenames, and the six cutoffs in `filter.thresholds`.
- **Role.** Keep a row only when every cutoff holds. A locus is written when at least one ingroup row passes, and only those rows are kept, plus the reference. `DAT` codes become species names. `header_fixes` then corrects known spelling (`echitum` to `echium` in this project).
- **Red flags.** Any one of these fails the row. The values are the ones used for the Boraginaceae screen:

  | Red flag | Cutoff that removes it |
  |---|---|
  | No exon left after gaps | `n_exons` not greater than 0 |
  | Too many intron breaks to treat as one marker | `n_introns` of 5 or more |
  | No exon long enough to be useful | `max_exon_bp` of 220 or less |
  | An `N` run longer than the screen allows | `max_intron_bp` of 121 or more |
  | Exons are short even if the longest one passes | `avg_exon_bp` of 152.51259 or less |
  | Intron masks are long on average | `avg_intron_bp` of 102.46 or more |

- **Okay to continue.** At least one ingroup sequence at the locus passes all six, and the locus is listed in `namelist.txt`. That file is what you import into Geneious.
- **Output.** `filtered/<run>/`, with one alignment per kept locus, `namelist.txt`, `filtered_rows.csv`, and `summary.txt`.

### 10. Report

- **Purpose.** See whether those six cutoffs actually separate awkward loci from the ones worth opening, before anyone spends time in Geneious.
- **Input.** The metric tables and the cutoffs in the config.
- **Role.** Count rows and loci that pass, and draw the six distributions with the cutoff marked.
- **Red flags.** Almost every locus passes, so the screen is not removing the structures you meant to avoid. Or almost nothing passes, so the cutoffs are tighter than the genes you have.
- **Okay to continue.** The kept set is small enough to curate and still contains the loci you want. Editing a cutoff means changing `filter.thresholds` and rerunning from stage 9. The page does not store a new cutoff.
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
