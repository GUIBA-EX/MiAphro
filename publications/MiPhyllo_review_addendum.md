# MiPhyllo: computational review addendum

This addendum accompanies the 5 November 2026 curated preprint and does not replace its primary catalog or claim experimental validation.

## Methods

All 48 retained pairs were rescanned against their own frozen target and background panels using independent C++ and Python implementations. Homologous projected terminal 1/3/5 compatibility was summarized per end and record. Primer3-py 2.4.0 was used to expand all IUPAC variants and calculate hairpin, homodimer, within-direction heterodimer and forward/reverse heterodimer temperatures under 50 nM per variant and 200 nM per direction. We recorded the difference between minimum primer Tm and maximum modeled structure Tm. Clean recovered internal inserts were compared by exact identity against inserts recovered by the same assay.

## Results

All 48 pairs passed the basic structure gates, while eight had a minimum margin below 5°C in at least one scenario, three with negative margins. Seven of 24 D pairs lost reference records under the terminal-triplet condition relative to an otherwise equivalent unconstrained-terminal paired scan. Among recovered inserts, for example, 37 of 56 Aphroditiformia 28S MD records were unique against the other recovered species labels. These denominators do not include undetected competitors.

A finite single-substitution screen of the eight low-margin pairs held terminal triplets fixed. One Nereididae 16S LD candidate, F `TTAATTGAAGGCYRGTATGA` and R `CAATGTGATCTCTCTARY`, yielded a 5.67°C minimum margin versus −2.08°C for the original. It matched the same 58/61 target records and the same 11 exactly; strict WC/GT background counts were 30/58 for both, while relaxed WC/GT counts were 166/285 and 166/282, respectively. An independent implementation confirmed target accession sets and background counts. This is an exploratory computational candidate, not a replacement in the released catalog.

## Limitations and data

The retrospective search uses the design libraries; there is no independently sampled test set or wet-lab confirmation. Reference insert uniqueness is conditional on recovery, one reference per species and clean inserts; it is not full-library species resolution or a barcode-gap analysis. Some homologous windows cannot distinguish indels from incompleteness. The one-substitution search is finite and cannot rule out better shifted-window or clade-specific alternatives. Complete tables, source code, settings and terminal 1/3/5 audits are in `audit/review_20261107/` and `scripts/review_*.py`.
