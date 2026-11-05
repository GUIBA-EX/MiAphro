# MiPhyllo · Joint terminal degeneracy optimization

[Current catalogue](../../README.en.md) · [Comparison data](comparison.json) · [Updated sequences](selected_results.csv) · [Previous PDF](pre_optimization_preprint.pdf)

Nine D pairs were evaluated. Four were replaced; five were retained. Ordinary sequences, lengths and all four endpoints are unchanged. Each oligo has degeneracy ≤8, restricted to its 3′ half. Replacement requires no loss in target coverage or exact matching and improvement in at least one, plus the original thermodynamic, background and close-relative gates.

Joint redistribution resolves previous ambiguities and reallocates them simultaneously. Seeds include canonical expansions and frequent target-window suffixes. Up to three 3′-half positions are jointly assigned masks; terminal-triplet masks are also enumerated within budget. Coverage- and exact-oriented shortlists are combined and re-scanned on raw sequences. This is an expanded staged search, not a global optimum proof.

|D pair|F/R candidates searched|Coverage|Exact|Strict WC/GT|Loose WC/GT|Terminal exclusions|Decision|
|---|---:|---:|---:|---|---|---:|---|
|Aphroditiformia-28S-MD|72,703/100,464|58→58/68|53→53/68|1/2→1/2|17/2→17/2|1→1|retain|
|Glyceriformia-16S-SD|63,957/43,835|14→15/20|1→1/20|14/33→13/33|113/104→107/105|1→0|replace|
|Nereididae-16S-LD|119,190/39,861|58→58/61|11→11/61|30/58→30/58|166/285→166/285|1→1|retain|
|Nereididae-28S-SD|96,492/80,430|29→30/52|8→19/52|0/0→0/0|0/0→0/0|7→5|replace|
|Phyllodociformia-16S-LD|107,108/99,790|40→41/59|3→3/59|0/1→0/0|76/29→80/31|2→2|replace|
|Phyllodociformia-16S-SD|4,160/56,040|49→50/59|11→11/59|0/8→0/12|81/11→82/15|1→0|replace|
|Phyllodociformia-28S-MD|11,689/86,668|52→52/63|5→5/63|21/26→21/26|95/42→95/42|5→5|retain|
|Phyllodociformia-28S-SD|54,053/15,417|51→51/63|4→4/63|7/23→7/23|53/53→53/53|5→5|retain|
|Syllidae-16S-SD|74,215/82,665|110→110/180|7→7/180|1/2→1/2|15/2→15/2|16→16|retain|

Syllidae 16S SD: 74,215 forward and 82,665 reverse allocations were generated. No shortlisted pair improved coverage/exact matching without a loss in the other metric. Its sequence is retained; terminal mismatch exclusions remain 16 species.

“Terminal exclusions” counts reference species gained by removing only terminal compatibility from the same scan; it is not a measure of experimental amplification. Every species counted as covered retains compatible terminal three bases at both ends.

[Current terminal audit](current_terminal/summary.csv) · [Species-level terminal audit](current_terminal/species.csv) · [Updated species outcomes](selected_species.csv)

## Reproduction

Install biopython and primer3-py. `joint_search.py` regenerates expanded candidates using the frozen `baseline_catalog.csv`, the preserved initial `candidates.json` lists and repository reference panels. Remove only the generated `joint_pairs.json` files in a temporary copy to rerun search. `background.py` scans target/background panels with the compiled repository scanner, applies full gates and independently checks selected pairs in Python. `audit_terminal_current.py` audits the published current catalogue. Candidate-generation and screening details are retained per pair.

No experimental results are added. Background counts are not amplification probabilities; modest increases in some accepted candidates are shown explicitly.
