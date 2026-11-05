![MiPhyllo](figures/MiPhyllo_overview.png)

# MiPhyllo: taxon-stratified 16S and 28S PCR primer candidates for environmental DNA metabarcoding of Phyllodocida

**Fei Xia**

[中文](README.md) | English

[Preprint PDF](publications/MiPhyllo_preprint.pdf) · [Manuscript](publications/MiPhyllo_preprint.md)

## Overview

The curated catalogue contains **48 primer pairs (96 oligos), arranged in 24 same-coordinate ordinary/D sets across six target units**: three suborders and three families. The original search considered eight units. Each retained unit keeps its own species-deduplicated reference and background panels. Phyllodocida remains the project scope; no order-wide primer is retained.

[Current FASTA](primers.fasta) · [Catalogue](catalog.csv) · [Species outcomes](species.csv) · [Methods](docs/METHODS.en.md) · [Archive and selection](archive/README.en.md)

The banner illustrates taxonomic context. Nereidiformia is represented by three family panels, not a retained suborder-wide assay.

## Retained designs

|Target|16S lengths|28S lengths|16S references|28S references|Results|
|---|---|---|---:|---:|---|
|Aphroditiformia|S|S/M/L|94|68|[Aphroditiformia](panels/Aphroditiformia/README.en.md)|
|Glyceriformia|S/M|S/M/L|20|12|[Glyceriformia](panels/Glyceriformia/README.en.md)|
|Phyllodociformia|S/L|S/M|59|63|[Phyllodociformia](panels/Phyllodociformia/README.en.md)|
|Syllidae|S|—|180|28|[Syllidae](panels/Syllidae/README.en.md)|
|Nereididae|S/M/L|S|61|52|[Nereididae](panels/Nereididae/README.en.md)|
|Hesionidae|S/M/L|S/M/L|53|45|[Hesionidae](panels/Hesionidae/README.en.md)|

Reference counts describe the available libraries; a count does not imply a retained primer for that gene. Each listed length includes its ordinary and D pair.

## Selection and naming

Retain a standard ordinary/D set when D covers at least 50% of its reference species and exactly matches at least one, under its original full background selection profile. Phyllodocida and Nereidiformia broad-target designs and all relaxed-background alternatives are archived. This is an editorial catalogue rule, not an experimentally established efficacy threshold. Ordinary primers are retained as same-coordinate comparators; the D criterion does not imply that ordinary coverage is also ≥50%.

`MiPhyllo-{taxon}-{gene}-{S|M|L}[D]_{F|R}`. Ordinary primers contain only A/C/G/T; D marks degeneracy in at least one oligo. Each oligo has degeneracy ≤8, restricted to its 3′ half. Ordinary/D counterparts share all four binding endpoints. S/M/L denote length classes, not mandatory sets.

## Results

![D coverage and exact matching](figures/catalog_coverage.svg)

Filled markers show D coverage; open markers show exact matching. Grey connectors show the gap between them. Missing lengths were removed or had no original candidate. Percentages use each unit’s own full reference denominator and are not common-panel comparisons.

- **Aphroditiformia**: Only S is retained for 16S; S/M/L are retained for 28S. The 16S SD covers 62/94 species but exactly matches only one.
- **Glyceriformia**: 16S retains S/M and 28S retains S/M/L; both 16S D pairs exactly match only one species. Lacydoniidae is absent from the reference panel.
- **Phyllodociformia**: 16S retains S/L and 28S retains S/M. The references represent Phyllodocidae only.
- **Syllidae**: Only 16S S/SD is retained; 28S and the other lengths are archived.
- **Nereididae**: 16S retains S/M/L; 28S retains S only. The 16S LD has the highest coverage but more background matches than SD.
- **Hesionidae**: Both genes retain S/M/L. The 28S SD covers 40/45 species with eight exact matches.

## Terminal optimization

Joint redistribution screened nine D pairs. Four replacements improve coverage without losing exact matches; five pairs, including Syllidae 16S SD, retain their sequences. All replacements retain coordinates, per-oligo degeneracy ≤8 and the original full background gates.

[Full comparison and methods](optimization/terminal_2026-11-05/README.en.md)

|D pair|Covered before → after|Exact before → after|Strict WC/GT before → after|Loose WC/GT before → after|Decision|
|---|---:|---:|---|---|---|
|Aphroditiformia-28S-MD|58 → 58 / 68|53 → 53 / 68|1/2 → 1/2|17/2 → 17/2|retain|
|Glyceriformia-16S-SD|14 → 15 / 20|1 → 1 / 20|14/33 → 13/33|113/104 → 107/105|replace|
|Nereididae-16S-LD|58 → 58 / 61|11 → 11 / 61|30/58 → 30/58|166/285 → 166/285|retain|
|Nereididae-28S-SD|29 → 30 / 52|8 → 19 / 52|0/0 → 0/0|0/0 → 0/0|replace|
|Phyllodociformia-16S-LD|40 → 41 / 59|3 → 3 / 59|0/1 → 0/0|76/29 → 80/31|replace|
|Phyllodociformia-16S-SD|49 → 50 / 59|11 → 11 / 59|0/8 → 0/12|81/11 → 82/15|replace|
|Phyllodociformia-28S-MD|52 → 52 / 63|5 → 5 / 63|21/26 → 21/26|95/42 → 95/42|retain|
|Phyllodociformia-28S-SD|51 → 51 / 63|4 → 4 / 63|7/23 → 7/23|53/53 → 53/53|retain|
|Syllidae-16S-SD|110 → 110 / 180|7 → 7 / 180|1/2 → 1/2|15/2 → 15/2|retain|

![Before and after](figures/terminal_optimization.svg)

## Reproducibility

[Selection decisions](archive/selection_decisions.csv) · [Curation audit](audit/curation.json) · [Independent scan](audit/catalog_rescan.json) · [Provenance](provenance.csv) · [Checksums](SHA256SUMS)

The catalogue retains 48 pairs: 44 are unchanged from the original independently scanned collection; four D pairs were replaced after joint degeneracy redistribution and independent target/background verification. Ordinary primers and all binding endpoints are unchanged. Original results and previous preprints remain archived. No experimental validation is added.

```bash
python3 -m pip install -r requirements.txt
python3 scripts/build_catalog.py
python3 scripts/build_docs.py
python3 scripts/plot_catalog.py  # requires matplotlib
python3 publications/render_integrated.py  # requires pandoc and Chrome
python3 scripts/update_checksums.py
python3 scripts/check_release.py
```

## Historical study

[Aphroditiformia study](studies/Aphroditiformia/README.en.md) · [Historical PDF](studies/Aphroditiformia/publications/MiAphro_preprint.pdf)

The historical study preserves the original complete search and figures. Use the current catalogue above for the curated selection.

## Citation and copyright

Xia, F. (2026). *MiPhyllo: taxon-stratified 16S and 28S PCR primer candidates for environmental DNA metabarcoding of Phyllodocida*. [MiPhyllo](https://github.com/GUIBA-EX/MiPhyllo).

[Citation metadata](CITATION.cff) · [Copyright and permissions](COPYRIGHT.md)
