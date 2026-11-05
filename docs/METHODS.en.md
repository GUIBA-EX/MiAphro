# Methods and reproducibility

[Catalogue](../README.en.md) · [中文](METHODS.md)

## Design and naming

Each gene and target unit uses one frozen reference per named species. Units are nested and species overlap between panels. Direction-adjusted MAFFT alignments guide window searches. Ordinary/D counterparts share all four endpoints; degeneracy is ≤8 per oligo and limited to its 3′ half. Selected spans are displayed as S 100–200 bp, M 201–300 bp and L 301–480 bp. The original Aphroditiformia and early order-wide searches used narrower ranges; the source scripts preserve these original search spaces.

Two concentration scenarios (50 nM per variant and 200 nM per direction) require minimum predicted Tm≥45°C, F/R midpoint difference≤5°C and structure Tm<55°C. The latest expanded searches use 18–30 nt oligos and 1 nt start increments, with staged coverage- and exact-oriented candidate selection. Shared coordinates do not imply that an ordinary oligo is necessarily a member of its D counterpart's expansion.

## Matching models

Target coverage requires ≤2 mismatches per end and compatibility at each terminal three bases, product length 100–500 bp and internal sequence≥50 bp. Exact matches have zero mismatches at both ends; IUPAC-compatible bases match. All reference species remain in the denominator, including fragments that do not span a window. Background scans permit internal gaps of 0–1000 bp. Strict WC/GT allows ≤2 mismatches per end and compatible terminal three bases; loose WC allows ≤3 mismatches, loose GT≤2, without terminal constraints. GT additionally allows primer-oriented G:A and T:C compatibility. These are matching sensitivity models, not measured amplification probabilities.

## Background panels and selection

- Original Aphroditiformia study: a six-group polychaete background, with each strict WC/GT count≤10 and each loose count≤50. This is not an evaluation against the later expanded marine background.
- Archived Phyllodocida order-wide results: the expanded marine background, with each strict WC/GT count≤10 and each loose count≤50. The original selection and its corresponding background are retained.
- Three additional suborders and three families: marine references plus Phyllodocida outside the target unit. Strict WC≤1%, GT≤3%; loose WC≤5%, GT≤10%; each strict close-relative count≤max(2, floor(3%×close-relative N)).
- Latest Syllidae 16S and Nereididae 28S: the current profile does not reduce ordinary/D coverage or exact-match counts relative to the preceding selection. Ranking then favors D coverage, D exact matches, ordinary coverage, ordinary exact matches and fewer GT background hits. The archived conditional alternatives retained strict and close-relative gates but removed loose-model hard vetoes; their original counts remain available.

The later marine background uses existing WoRMS marine-genus screening and deduplicated references, rather than complete species- and locality-level verification. Each unit retains its own panel and denominator. Reference-species coverage is not coverage of every described species in that taxon, and percentages are not common-panel comparisons.

## Files and reproduction

`panels/` contains current results, FASTA, species outcomes and frozen references. `conditional/` is an archive pointer with empty current exports. `archive/` preserves all original exports, the previous integrated PDF and per-set decisions. `studies/Aphroditiformia/` preserves the original study, base-matching figures, WoRMS coverage and preprint. `analysis/` contains compact source scripts, frozen inputs and candidate score tables; large per-candidate hit matrices, caches and raw downloaded GenBank text are regenerated rather than duplicated. `provenance.csv` maps unified names to source rows; `analysis/source_manifest.json` records imported-file hashes.

Install `requirements.txt`. Run `python3 scripts/build_catalog.py`, `python3 scripts/build_docs.py`, then `python3 scripts/check_release.py`. Full verification with `python3 scripts/recheck_catalog.py` requires a C++17 compiler and compares new C++ scans with an independent Python implementation record by record; a new run writes `audit/current_rescan.json` and preserves the original full scan. `python3 scripts/plot_catalog.py` regenerates the overview using matplotlib.

Complete search regeneration additionally needs MAFFT and primer3-py. Keep the source directory layout, compile the corresponding `raw_scan.cpp`, and run the design scripts on the frozen inputs. Regenerate full candidate hit matrices before invoking legacy finish scripts that depend on them. Unified-catalogue verification does not need the omitted intermediate matrices.

The current 48 pairs are unchanged subsets of the 102 independently checked original pairs. Full original scan records remain intact; release checks verify retained identity and selection decisions. The integrated preprint reports the curated catalogue; the original Aphroditiformia preprint remains a historical study.

## Catalogue curation on 5 November 2026

Retain standard ordinary/D sets with D coverage ≥50% and at least one exact reference match, under the original full background profile. Withdraw Phyllodocida and Nereidiformia broad-target designs and archive all relaxed-background alternatives. The result is six target units and 24 ordinary/D sets (48 pairs; 96 oligos); 54 pairs are archived. Ordinary primers remain same-coordinate comparators, not necessarily designs reaching the D coverage criterion.

The 50% and nonzero-exact rules are operational editorial criteria, not experimentally established efficacy thresholds. Exact matching remains low for some retained designs and is shown alongside coverage. Different lengths address different fragment-length requirements; a length is not removed solely because another length has higher coverage. Every decision is recorded in `archive/selection_decisions.csv`. No sequences or original match statistics are changed.
