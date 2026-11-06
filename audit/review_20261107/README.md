# 2026-11-07 computational review

Frozen input catalog, primary references and sequences are unchanged. `all_pair_structures.csv` contains 48 pairs in two concentration scenarios, all concrete variants, worst hairpin/homodimer/within-direction/F-R heterodimer and participating variant sequences. Primer3-py 2.4.0; 50 mM monovalent, 1.5 mM Mg, 0.6 mM dNTP, no adapters. A negative margin indicates a model warning, not an experimental failure.

`terminal_windows_records.csv` and `terminal_windows_summary.csv` report homologous projected terminal 1/3/5 compatibility for both ends and the same record. Ambiguous bases and alignment-window indel or incompleteness are reported separately from called-base mismatches; indel and incomplete windows are not distinguished from each other. Homologous-site compatibility is not identical to an unrestricted paired scan.

`recovered_insert_identity.csv` and `recovered_inserts.csv` report identity uniqueness against other species labels recovered by the same assay. Records with no clean insert and records with multiple distinct clean inserts are excluded from the single-insert denominator. Undetected competing taxa and within-species variation are not assessed. This is not full-library taxonomic resolution or external validation.

Reproduce from repository root: `python scripts/recheck_catalog.py`, then `python scripts/review_terminal_windows.py`, `python scripts/review_structures.py` and `python scripts/review_insert_identity.py`. The rescan checks C++ against independent Python on each target and background record using its own panel and writes `audit/current_rescan.json`; it also writes temporary scanner output to `.recheck/`. Terminal-three compatibility and losses are in `optimization/terminal_2026-11-05/current_terminal/`.

No new replacement primer, independent test set, expanded non-target sampling or experiment is claimed. Historical search and original primer evidence are preserved.

## Low-margin candidate search

`low_margin_search_candidates.csv` and `low_margin_search_summary.csv` record a finite single-base substitution search of the eight pairs with modeled margins below 5 °C. The final three 3′ bases of each primer and the other primer were held fixed; candidates had to reach a ≥5 °C minimum margin in both scenarios, pass original chemistry and background gates, and lose neither target coverage count nor exact-match count. The search did not explore shifted endpoints, multi-position changes, all clades or all possible primer sequences.

One exploratory Nereididae 16S LD candidate passed these gates: F `TTAATTGAAGGCYRGTATGA`; R `CAATGTGATCTCTCTARY`. Its minimum margin was 5.67 °C versus −2.08 °C for the original. Both cover the same 58 of 61 targets and exactly match the same 11; strict WC/GT background counts remain 30/58 and loose WC/GT change from 166/285 to 166/282. These are frozen-panel computational comparisons. An independent Python implementation checked all target and background counts, near-relative strict hits, and accession gains and losses (`low_margin_independent_verification.json`). It is not a validated replacement; the official catalog remains unchanged. Seven other pairs had no qualifying candidate in this finite search.

Reproduce with `python scripts/review_low_margin_search.py` followed by `python scripts/review_candidate_verify.py`; this is an additional retrospective search and not an independent biological test.
