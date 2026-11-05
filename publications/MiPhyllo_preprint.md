MiPhyllo: taxon-stratified 16S and 28S PCR primer candidates for
environmental DNA metabarcoding of Phyllodocida

<div class="section page">

# MiPhyllo: taxon-stratified 16S and 28S PCR primer candidates for environmental DNA metabarcoding of Phyllodocida

Fei Xia

Computational preprint · 5 November 2026

## Abstract

Taxonomic breadth and primer specificity can conflict when designing
environmental DNA assays for diverse marine annelids. We integrated
mitochondrial 16S and nuclear 28S designs across eight nested target
units within Phyllodocida: the order, four suborders and three families.
Frozen, species-deduplicated panels were searched for short, medium and
long windows. Ordinary and degenerate counterparts shared binding
endpoints, with degeneracy limited to eight variants per oligo. The
catalogue comprises 94 current pairs and eight conditional alternatives.
All 102 pairs were re-scanned using independent C++ and Python
implementations with agreement in target and background outcomes.
Family-level designs improved focal coverage: Nereididae 16S LD covered
58/61 references, while Hesionidae 28S SD covered 40/45. Order-wide
coverage remained low. Latest Nereididae 28S SD covered 29/52 references
with eight exact matches. A conditional Syllidae 16S S2D increased
coverage from 110 to 126/180 and exact matches from seven to 48, but
exhibited substantially greater loose-model background matching. These
results support taxon-stratified candidate selection rather than a
single universal assay. They represent computational matching, not
experimentally measured amplification performance.

## 1. Introduction

Phyllodocida comprises divergent polychaete lineages for which a single
short rRNA assay may not reconcile coverage and specificity. Scale-worm
phylogenetic reference data provide an established starting point \[1\],
while environmental DNA has been applied to marine metazoan surveys
\[2\]. Broad-target primer designs such as MiFish illustrate the value
of conserved binding regions surrounding variable sequence \[3\]; this
study applies that general principle to annelid reference panels,
without transferring fish-primer sequences or assuming equivalent
performance.

We compare designs at the order, suborder and family levels, retain
multiple amplicon lengths, and explicitly separate strict and loose
background matching. MiPhyllo denotes the integrated Phyllodocida
collection. The repository and current sequence identifiers use
MiPhyllo; the original Aphroditiformia preprint remains a separate case
study.

</div>

<div class="section page">

## 2. Materials and methods

### 2.1 Frozen references and taxonomic units

Each gene panel retained one reference per named species, prioritizing
sequence completeness during source curation. Target units were
Phyllodocida; Aphroditiformia, Glyceriformia, Nereidiformia and
Phyllodociformia; and Syllidae, Nereididae and Hesionidae. Panels
overlap taxonomically and must not be summed as independent sampling.
Direction-adjusted MAFFT alignments \[4\] guided homologous windows.
Reference accessions, family assignments, alignments and frozen
backgrounds are distributed with each unit.

### 2.2 Candidate generation

Current selected spans fall within S 100–200 bp, M 201–300 bp and L
301–480 bp. Earlier searches used narrower ranges, retained in source
scripts. Latest expanded searches evaluated 18–30 nt oligos at 1 nt
start increments, combining frequent haplotypes and consensus seeds,
coverage-oriented and exact-oriented branches, spatially distributed
windows and local refinement. Ordinary/D pairs share four binding
endpoints; each oligo has degeneracy ≤8, restricted to its 3′ half. An
ordinary sequence need not be nested within the D expansion.

Two concentration scenarios—50 nM per variant and 200 nM per
direction—required minimum predicted Tm ≥45°C, F/R midpoint difference
≤5°C and predicted secondary-structure Tm \<55°C. Search procedures were
staged heuristics, not exhaustive enumeration of all possible sequences.

### 2.3 Matching and selection

Target coverage required ≤2 mismatches per end, compatibility at each
terminal three bases, total product length 100–500 bp and internal
length ≥50 bp. Exact matching required zero mismatches at both ends,
including IUPAC-compatible bases. All reference species remained in the
denominator, including fragments lacking a complete window. Background
scans allowed an internal gap of 0–1000 bp. Strict WC/GT used ≤2
mismatches and terminal compatibility; loose WC allowed ≤3 mismatches
and loose GT ≤2 without terminal constraints. GT additionally permits
primer-oriented G:A and T:C compatibility and is a sensitivity model.

The original Aphroditiformia study used a six-group polychaete
background; Phyllodocida used an expanded marine background. Both
retained absolute limits of ten strict and fifty loose hits per model.
Later suborder/family designs required strict WC ≤1%, GT ≤3%, loose WC
≤5% and GT ≤10%, plus each strict close-relative count ≤max(2, floor(3%
× close-relative N)). Marine screening was based on existing genus-level
WoRMS assignments, not comprehensive species/locality verification.

Latest balanced selections did not reduce ordinary/D coverage or exact
matching relative to preceding selections. Conditional alternatives
retain strict and close-relative gates but remove loose-model hard
vetoes. They are separate candidates, identified by 2 (e.g., S2/S2D),
not two-pair mixtures.

### 2.4 Independent verification

All 102 selected pairs were freshly scanned on their own frozen panels,
comparing C++ and independent Python outcomes record by record. Sequence
exports, coordinates, degeneracy, species totals, bilingual tables and
checksums were also checked. No wet-lab validation, common-panel
comparison or independent holdout experiment is reported here.

</div>

<div class="section page">

## 3. Results

### 3.1 Reference panels and coverage

| Target unit      | 16S target | 16S background | 28S target | 28S background |
|------------------|------------|----------------|------------|----------------|
| Aphroditiformia  | 94         | 381            | 68         | 214            |
| Phyllodocida     | 515        | 3340           | 302        | 2085           |
| Glyceriformia    | 20         | 3835           | 12         | 2375           |
| Nereidiformia    | 315        | 3540           | 133        | 2254           |
| Phyllodociformia | 59         | 3796           | 63         | 2324           |
| Syllidae         | 180        | 3675           | 28         | 2359           |
| Nereididae       | 61         | 3794           | 52         | 2335           |
| Hesionidae       | 53         | 3802           | 45         | 2342           |

**Table 1.** Frozen reference species counts. Nested panels overlap;
background definitions differ between design units.

![](../figures/catalog_coverage.png)

**Figure 1.** Current D coverage within each reference panel. S, M and L
denote selected length classes. Missing Syllidae 28S L indicates no
passing pair. Percentages are not common-panel performance comparisons.

Phyllodocida designs covered at most 71/515 16S and 53/302 28S
references. Nereidiformia remained limited (62/315 and 45/133).
Phyllodociformia reached 49/59 and 52/63, but its references represented
Phyllodocidae only. Glyceriformia 16S SD/MD covered 14/20, with only one
exact match each. Detailed ordinary/D outcomes and sequences for all
units appear in Appendix A.

</div>

<div class="section page">

## 3.2 Family-level optimization and conditional alternatives

Nereididae 16S LD covered 58/61 references, but its strict WC/GT
background counts (30/58) exceeded those of SD (0/12), which covered
57/61. Latest 28S SD increased coverage from 28 to 29/52 while retaining
eight exact matches; all four strict/loose background counts were zero.
The 28S MD and LD retained low current coverage (11/52 and 4/52).

Syllidae current 16S SD, MD and LD covered 110, 104 and 87 of 180
references, respectively. No ordinary/D replacement improved or retained
all four target metrics under the original complete background profile.
Conditional S2D covered 126/180 with 48 exact matches, compared with 110
and seven for SD. However, strict WC/GT hits were 33/42 of 3,675
background references and loose hits were 1,007/601. Thus the apparent
gain depends on how mismatch tolerance is treated. The conditional set
is not merged into the current FASTA.

## 4. Discussion

The data favor stratified target design: broad order-level candidates
show poor reference coverage, whereas several narrower panels produce
higher within-panel coverage. This pattern does not establish
superiority in a common sample or reference panel. Some gains reflect
restricted target scope, different reference availability and a
different near-relative background.

Degeneracy can increase compatible target diversity without increasing
the number of primer pairs, but exact matching, background compatibility
and the concentration of each variant remain distinct considerations. A
high exact-match count does not rescue a low-coverage design, and low
strict-background counts do not imply low loose-model matching. The
conditional Syllidae designs make this tradeoff explicit.

Reference limitations include incomplete fragments, uneven taxonomic
representation and inherited taxonomy. Lacydoniidae is absent from
Glyceriformia references, Microphthalmidae from Nereidiformia
references, and Phyllodociformia includes only Phyllodocidae. Syllidae
28S remains poorly covered. Neither all described species nor the full
marine environmental community is represented. Primer-site matching also
does not measure taxonomic resolution of the intervening amplicon.

## 5. Conclusion and data availability

MiPhyllo provides a traceable, taxon-stratified collection of 16S/28S
candidates and explicit alternatives for eight Phyllodocida target
units. Candidate use should consider target scope, exact matching and
both background models. Frozen references, 94 current pairs, eight
conditional pairs, per-species outcomes, source algorithms and
independent re-scan records are available at
[github.com/GUIBA-EX/MiPhyllo](https://github.com/GUIBA-EX/MiPhyllo).
Full canonical identifiers follow
MiPhyllo-{target}-{gene}-{length}\[2\]\[D\]\_{F\|R}; these identifiers
are consistent with the repository and manuscript name.

## References

\[1\] Norlinder E et al. (2012). Phylogeny of scale-worms
(Aphroditiformia, Annelida), assessed from 18SrRNA, 28SrRNA, 16SrRNA,
mitochondrial cytochrome c oxidase subunit I (COI), and morphology.
*Molecular Phylogenetics and Evolution*. PMID: 22789762.

\[2\] Environmental DNA metabarcoding for monitoring metazoan
biodiversity in Antarctic nearshore ecosystems (2021). *PeerJ* 9:e12458.
doi:10.7717/peerj.12458.

\[3\] Miya M et al. (2015). MiFish, a set of universal PCR primers for
metabarcoding environmental DNA from fishes: detection of more than 230
subtropical marine species. *Royal Society Open Science* 2:150088.
doi:10.1098/rsos.150088.

\[4\] Katoh K, Standley DM (2013). MAFFT multiple sequence alignment
software version 7: improvements in performance and usability.
*Molecular Biology and Evolution* 30:772–780. doi:10.1093/molbev/mst010.

</div>

<div class="section page">

## Appendix A1. Aphroditiformia

Canonical prefix: MiPhyllo-Aphroditiformia-. All sequences are 5′→3′.
Covered and exact counts use the full target reference denominator.
Strict and loose background values are WC/GT.

| Gene/pair | bp  | Covered | Exact | Strict bg | Loose bg |
|-----------|-----|---------|-------|-----------|----------|
| 16S-S     | 138 | 58/94   | 1/94  | 0/0       | 23/0     |
| 16S-SD    | 138 | 62/94   | 1/94  | 0/0       | 25/0     |
| 16S-M     | 210 | 28/94   | 24/94 | 0/0       | 0/0      |
| 16S-MD    | 210 | 44/94   | 24/94 | 2/9       | 12/9     |
| 16S-L     | 329 | 28/94   | 24/94 | 0/0       | 0/0      |
| 16S-LD    | 329 | 41/94   | 24/94 | 0/7       | 7/7      |
| 28S-S     | 148 | 54/68   | 29/68 | 0/0       | 16/2     |
| 28S-SD    | 148 | 59/68   | 47/68 | 2/2       | 19/2     |
| 28S-M     | 250 | 56/68   | 50/68 | 0/0       | 16/0     |
| 28S-MD    | 250 | 58/68   | 53/68 | 1/2       | 17/2     |
| 28S-L     | 368 | 58/68   | 22/68 | 0/0       | 15/0     |
| 28S-LD    | 368 | 59/68   | 33/68 | 0/1       | 18/1     |

### Primer sequences

| Gene/pair | Forward                | Reverse                |
|-----------|------------------------|------------------------|
| 16S-S     | ATAAACTACCTTAGGGATAACA | GAGATTTAACAGTTGAACAAAC |
| 16S-SD    | ATAAACTACCTYRGGGATAACA | GAGATTTAACAGTTGAACRAAC |
| 16S-M     | AAGACAAGAAGACCCTATAG   | ACCATCCTATCAATAAGAGCTC |
| 16S-MD    | AAGACAAGAAGACCCYRTWG   | ACCATCCTATCAATWWGARCTC |
| 16S-L     | GTAGCATAATCATTTGCCTT   | ACCATCCTATCAATAAGAGC   |
| 16S-LD    | GTAGCATAATMAYTTGCCTY   | ACCATCCTATCAATWWGARC   |
| 28S-S     | CTCAGAGGCAAACTGATGGGAC | GCTCTCCTCGTGGAAAGGGC   |
| 28S-SD    | CTCAGAGGCAAACYGATGGGAC | GCTCTCCTCGYRGAAAGKGC   |
| 28S-M     | AATCCCCTGTCGCTGCGGCA   | TCCCTCACGGTACTTGTCCGCT |
| 28S-MD    | AATCCCCTGTCGYTGCGGCA   | TCCCTCACGGTACTTGTCCKCT |
| 28S-L     | CTCAGAGGCAAACTGATGGGAC | CGGAACACACTGTAAACAGTCC |
| 28S-LD    | CTCAGAGGCAAACYGATGGGAC | CGGAACACACTGTRAACAGTCC |

Binding coordinates, accession identifiers and per-species matches:
panels/Aphroditiformia/results.csv and species.csv in the repository.
Ordinary/D pairs have identical endpoints within each design.

</div>

<div class="section page">

## Appendix A2. Phyllodocida

Canonical prefix: MiPhyllo-Phyllodocida-. All sequences are 5′→3′.
Covered and exact counts use the full target reference denominator.
Strict and loose background values are WC/GT.

| Gene/pair | bp  | Covered | Exact  | Strict bg | Loose bg |
|-----------|-----|---------|--------|-----------|----------|
| 16S-S     | 178 | 36/515  | 13/515 | 0/3       | 0/3      |
| 16S-SD    | 178 | 40/515  | 16/515 | 0/6       | 0/6      |
| 16S-M     | 256 | 69/515  | 10/515 | 2/5       | 34/7     |
| 16S-MD    | 256 | 71/515  | 10/515 | 3/6       | 36/8     |
| 16S-L     | 380 | 69/515  | 14/515 | 3/4       | 30/5     |
| 16S-LD    | 380 | 70/515  | 14/515 | 3/4       | 33/5     |
| 28S-S     | 136 | 39/302  | 16/302 | 0/0       | 0/0      |
| 28S-SD    | 136 | 44/302  | 24/302 | 0/0       | 0/0      |
| 28S-M     | 262 | 45/302  | 15/302 | 1/2       | 4/3      |
| 28S-MD    | 262 | 53/302  | 22/302 | 4/6       | 15/9     |
| 28S-L     | 320 | 46/302  | 26/302 | 0/0       | 0/0      |
| 28S-LD    | 320 | 51/302  | 38/302 | 0/0       | 4/0      |

### Primer sequences

| Gene/pair | Forward                | Reverse                |
|-----------|------------------------|------------------------|
| 16S-S     | AGGACACACAAGAATAAAGCTG | TGTCTTTTGGTTGTCCCAAC   |
| 16S-SD    | AGGACACACAAGRATWAARCTG | TGTCTTTTGGTYGYCCCARC   |
| 16S-M     | AAGGACAAGAAGACCCTATAGA | ATAGTTACACTAAGCCAACATC |
| 16S-MD    | AAGGACAAGAAGACCCTRTNGA | ATAGTTACACTRAGCCAACATC |
| 16S-L     | CAAAGGTAGCATAATCATTTGC | ATAGTTACACTAAGCCAACATC |
| 16S-LD    | CAAAGGTAGCAYAATMAYTTGC | ATAGTTACACTRAGCCAACATC |
| 28S-S     | GTGTTATAGTCCCTGGGTTC   | GGAACACACTGTAAACAGTCCA |
| 28S-SD    | GTGTTATAGTCCCYGGKKTC   | GGAACACACTGTRRACAGTCCW |
| 28S-M     | AGAGGCAAACGGATAGGACC   | GCCACCCGAGGCTATAACAC   |
| 28S-MD    | AGAGGCAAACKGAYRGGACC   | GCCACCCGAGTCTATAMYRC   |
| 28S-L     | GGAGGGCGCCCTTTCCACGAGG | CGATCGATTCGCCACGGACCCT |
| 28S-LD    | GGAGGGCGCCCTTTCCRYGAGS | CGATCGATTCGCCASNGACCCT |

Binding coordinates, accession identifiers and per-species matches:
panels/Phyllodocida/results.csv and species.csv in the repository.
Ordinary/D pairs have identical endpoints within each design.

</div>

<div class="section page">

## Appendix A3. Glyceriformia

Canonical prefix: MiPhyllo-Glyceriformia-. All sequences are 5′→3′.
Covered and exact counts use the full target reference denominator.
Strict and loose background values are WC/GT.

| Gene/pair | bp  | Covered | Exact | Strict bg | Loose bg |
|-----------|-----|---------|-------|-----------|----------|
| 16S-S     | 110 | 8/20    | 0/20  | 2/3       | 60/27    |
| 16S-SD    | 110 | 14/20   | 1/20  | 14/33     | 113/104  |
| 16S-M     | 252 | 7/20    | 1/20  | 1/1       | 40/17    |
| 16S-MD    | 252 | 14/20   | 1/20  | 10/13     | 81/35    |
| 16S-L     | 316 | 5/20    | 1/20  | 0/0       | 2/0      |
| 16S-LD    | 316 | 8/20    | 2/20  | 0/7       | 10/8     |
| 28S-S     | 106 | 7/12    | 4/12  | 0/2       | 1/28     |
| 28S-SD    | 106 | 7/12    | 6/12  | 0/3       | 7/31     |
| 28S-M     | 236 | 5/12    | 4/12  | 0/4       | 90/35    |
| 28S-MD    | 236 | 7/12    | 4/12  | 0/4       | 93/35    |
| 28S-L     | 324 | 5/12    | 4/12  | 0/1       | 49/26    |
| 28S-LD    | 324 | 7/12    | 4/12  | 0/1       | 49/26    |

### Primer sequences

| Gene/pair | Forward                | Reverse                |
|-----------|------------------------|------------------------|
| 16S-S     | TAAATTTTGTTGGGGCGACA   | ATCCCTAAGGTAGCTATTTCCT |
| 16S-SD    | TAAATTTTGTTGGGGCGACA   | ATCCCTATGGTAGCTRWTTYCT |
| 16S-M     | CACGAAAATAAAGCTGTCTC   | CCCTAAGGTAGCTATTTCCT   |
| 16S-MD    | CACGAAAATAAWDCTGTCTC   | CCCTAAGGTAGCTRTTTYCT   |
| 16S-L     | AACTAGACTTCAAGTGAA     | AACCTTTCTAGGCTTCTG     |
| 16S-LD    | AACTAGACCYYARGTGAA     | AACCTTTCTWRGCTTCTR     |
| 28S-S     | ACGGCGAGTGAAGCGGGAAGAG | CTCGGACGCCAGAACGCTGC   |
| 28S-SD    | ACGGCGAGTGAAGCGGGAWGAG | CTCGGACGCCASAMCRCTGC   |
| 28S-M     | GCGGACAAGTACCGTGAG     | GAAAACCGGTCGTGGCTA     |
| 28S-MD    | GCGGACAAGTACCGTGAG     | GAAAACCGGTCGTGGCYA     |
| 28S-L     | CGTCCTTGGAGTCGGGTTGT   | CGGAAAACCGGTCGTGGCTA   |
| 28S-LD    | CGTCCTTGGAGTCGGGTTGT   | CGGAAAACCGGTCGTGGCYA   |

Binding coordinates, accession identifiers and per-species matches:
panels/Glyceriformia/results.csv and species.csv in the repository.
Ordinary/D pairs have identical endpoints within each design.

</div>

<div class="section page">

## Appendix A4. Nereidiformia

Canonical prefix: MiPhyllo-Nereidiformia-. All sequences are 5′→3′.
Covered and exact counts use the full target reference denominator.
Strict and loose background values are WC/GT.

| Gene/pair | bp  | Covered | Exact  | Strict bg | Loose bg |
|-----------|-----|---------|--------|-----------|----------|
| 16S-S     | 162 | 45/315  | 3/315  | 0/0       | 0/0      |
| 16S-SD    | 162 | 50/315  | 15/315 | 0/0       | 0/0      |
| 16S-M     | 212 | 55/315  | 8/315  | 2/4       | 36/7     |
| 16S-MD    | 212 | 57/315  | 8/315  | 3/7       | 46/12    |
| 16S-L     | 314 | 48/315  | 3/315  | 1/7       | 30/8     |
| 16S-LD    | 314 | 62/315  | 0/315  | 15/40     | 74/53    |
| 28S-S     | 194 | 23/133  | 4/133  | 0/0       | 0/1      |
| 28S-SD    | 194 | 31/133  | 0/133  | 0/6       | 53/55    |
| 28S-M     | 252 | 28/133  | 6/133  | 1/4       | 17/6     |
| 28S-MD    | 252 | 36/133  | 7/133  | 3/6       | 51/9     |
| 28S-L     | 324 | 27/133  | 6/133  | 0/1       | 22/6     |
| 28S-LD    | 324 | 45/133  | 6/133  | 11/26     | 100/32   |

### Primer sequences

| Gene/pair | Forward              | Reverse                |
|-----------|----------------------|------------------------|
| 16S-S     | TGGGACAACCTAAAGACA   | AATAGTTACACTAAGCCAACAT |
| 16S-SD    | TGGGACAACCWWAAGMCA   | AATGGTTACACTAAGCCAACAT |
| 16S-M     | AAGGACAAGAAGACCCTA   | AATAGTTACACTAAGCCAAC   |
| 16S-MD    | AAGGACAAGAAGACCCTA   | AATAGTTACAYWRAGCCAAC   |
| 16S-L     | GCCTTTTAATTGGAGGCTAG | AATGGTTACACTAAGCCA     |
| 16S-LD    | GCCTTTTAATTRGAGGMTWG | AATAGGTACVYTAAGCCA     |
| 28S-S     | GGGACCGCAGAGTCTGCC   | CCACCATCCCGCCGGCCT     |
| 28S-SD    | GGGACCGCAAAGTCTRSY   | CTACCATCCCSTSGGMCT     |
| 28S-M     | GACGAGCGCACTTTCCGC   | GAAGCACTGTGAACAGTCCT   |
| 28S-MD    | GACGAGCGCACTTTCYGB   | AAAGCACTGTRRACAGTCCW   |
| 28S-L     | GACGCGTGCACTTTCCGG   | GTGCCGACAGATTCGCCA     |
| 28S-LD    | GACGCGTGCACTTTCCRB   | GTGCCGACARAWTCGCYA     |

Binding coordinates, accession identifiers and per-species matches:
panels/Nereidiformia/results.csv and species.csv in the repository.
Ordinary/D pairs have identical endpoints within each design.

</div>

<div class="section page">

## Appendix A5. Phyllodociformia

Canonical prefix: MiPhyllo-Phyllodociformia-. All sequences are 5′→3′.
Covered and exact counts use the full target reference denominator.
Strict and loose background values are WC/GT.

| Gene/pair | bp  | Covered | Exact | Strict bg | Loose bg |
|-----------|-----|---------|-------|-----------|----------|
| 16S-S     | 118 | 31/59   | 2/59  | 0/0       | 10/3     |
| 16S-SD    | 118 | 49/59   | 11/59 | 0/8       | 81/11    |
| 16S-M     | 258 | 33/59   | 0/59  | 12/49     | 47/57    |
| 16S-MD    | 258 | 44/59   | 0/59  | 15/93     | 76/101   |
| 16S-L     | 308 | 27/59   | 2/59  | 0/0       | 9/4      |
| 16S-LD    | 308 | 40/59   | 3/59  | 0/1       | 76/29    |
| 28S-S     | 198 | 49/63   | 4/63  | 5/18      | 35/30    |
| 28S-SD    | 198 | 51/63   | 4/63  | 7/23      | 53/53    |
| 28S-M     | 298 | 12/63   | 6/63  | 0/0       | 44/12    |
| 28S-MD    | 298 | 52/63   | 5/63  | 21/26     | 95/42    |
| 28S-L     | 302 | 12/63   | 0/63  | 0/0       | 36/11    |
| 28S-LD    | 302 | 50/63   | 0/63  | 14/20     | 82/35    |

### Primer sequences

| Gene/pair | Forward                  | Reverse                  |
|-----------|--------------------------|--------------------------|
| 16S-S     | ATAGTTACCATAGGGATAACAGGC | ACCAACCCTTATAGACTGCTAC   |
| 16S-SD    | ATAGTTACCATAGGGATAACAGGC | ACCAACCCTTATARRCTGCTRC   |
| 16S-M     | GACGAAAAGACCCTATTGAG     | CTATAGGATATCCTAAGCCAACAT |
| 16S-MD    | GACGAAAAGACCCTATYGAK     | CTATAGGATATCYTAAGCCAACAT |
| 16S-L     | ACCTAAATTTTAATGAAAGACA   | TTGAACAAACCAACCCTT       |
| 16S-LD    | ACCTAAATTTTAATRWARGACA   | TTGAACAAACYWACCCTY       |
| 28S-S     | GAAAGCCACGACCGGTTC       | CACTGTGAACAGTCCTCGCCAG   |
| 28S-SD    | GAAAGCCACGACCGGTTC       | CACTGTGAACAGTCCTCVCCWG   |
| 28S-M     | GTCCTTGGAGTCGGGTTGTTTG   | AAAGTGCGCCCACCGACCAC     |
| 28S-MD    | GTCCTTGGAGTCGGGTTGTTTG   | AAAGTGCACCCACCGACNRC     |
| 28S-L     | GCGTCCTCGGAGTCGGGT       | GGAAAGTGCGCCCACCGACCAC   |
| 28S-LD    | GCGTCCTCGGAGTCGGGT       | GGAAAGTGCACCCACCGACNRC   |

Binding coordinates, accession identifiers and per-species matches:
panels/Phyllodociformia/results.csv and species.csv in the repository.
Ordinary/D pairs have identical endpoints within each design.

</div>

<div class="section page">

## Appendix A6. Syllidae

Canonical prefix: MiPhyllo-Syllidae-. All sequences are 5′→3′. Covered
and exact counts use the full target reference denominator. Strict and
loose background values are WC/GT.

| Gene/pair | bp  | Covered | Exact | Strict bg | Loose bg |
|-----------|-----|---------|-------|-----------|----------|
| 16S-S     | 142 | 55/180  | 2/180 | 0/2       | 4/2      |
| 16S-SD    | 142 | 110/180 | 7/180 | 1/2       | 15/2     |
| 16S-M     | 282 | 41/180  | 3/180 | 2/4       | 53/4     |
| 16S-MD    | 282 | 104/180 | 0/180 | 6/23      | 162/24   |
| 16S-L     | 308 | 39/180  | 3/180 | 0/0       | 17/4     |
| 16S-LD    | 308 | 87/180  | 7/180 | 2/11      | 169/41   |
| 28S-S     | 118 | 2/28    | 2/28  | 0/0       | 0/0      |
| 28S-SD    | 118 | 2/28    | 2/28  | 0/0       | 0/0      |
| 28S-M     | 202 | 1/28    | 1/28  | 0/0       | 0/0      |
| 28S-MD    | 202 | 1/28    | 1/28  | 0/0       | 0/0      |

### Primer sequences

| Gene/pair | Forward              | Reverse                      |
|-----------|----------------------|------------------------------|
| 16S-S     | AAAGACAGAAAGACCCTA   | CCCTACAGTAGCTTACTTAT         |
| 16S-SD    | AAAGACAGAAAGWCCCTH   | CCCTACAGTARCTTACTTRW         |
| 16S-M     | GCGAAGGTAGCATAATCATT | AAAATATTAGCCTGTTATCCCTACAGTA |
| 16S-MD    | GCGAAGGTAGCAYAATCAYK | AAAGTATTAGCCTGTTAYCCCTAVAGTA |
| 16S-L     | TTCAACGGCCGCGGTATC   | AAAGTATTAGCCTGTTATCCCTACAG   |
| 16S-LD    | TTCAACGGCCGCGGTABY   | AAAGTATTAGCCTGTTATCCCTANAG   |
| 28S-S     | CTGTGGCGTGTGGGACGG   | GGACCATAGGGAGGAGGCCGGA       |
| 28S-SD    | CTGTGGCGTGTGGGACGG   | GGACCATAGGGAGGAGRSCSGA       |
| 28S-M     | TGCGCCGGTGGCAGCATCTG | ACTTTCCCTCACGGTACTTGTTCGCTAT |
| 28S-MD    | TGCGCCGGTGGCAGCATCTG | ACTTTCCCTCACGGTACTTGTYCRCTAT |

Binding coordinates, accession identifiers and per-species matches:
panels/Syllidae/results.csv and species.csv in the repository.
Ordinary/D pairs have identical endpoints within each design.

</div>

<div class="section page">

## Appendix A7. Nereididae

Canonical prefix: MiPhyllo-Nereididae-. All sequences are 5′→3′. Covered
and exact counts use the full target reference denominator. Strict and
loose background values are WC/GT.

| Gene/pair | bp  | Covered | Exact | Strict bg | Loose bg |
|-----------|-----|---------|-------|-----------|----------|
| 16S-S     | 194 | 30/61   | 7/61  | 0/0       | 16/6     |
| 16S-SD    | 194 | 57/61   | 7/61  | 0/12      | 51/25    |
| 16S-M     | 212 | 52/61   | 8/61  | 5/9       | 48/12    |
| 16S-MD    | 212 | 52/61   | 8/61  | 8/14      | 66/19    |
| 16S-L     | 412 | 38/61   | 6/61  | 0/40      | 66/106   |
| 16S-LD    | 412 | 58/61   | 11/61 | 30/58     | 166/285  |
| 28S-S     | 150 | 19/52   | 1/52  | 0/0       | 0/0      |
| 28S-SD    | 150 | 29/52   | 8/52  | 0/0       | 0/0      |
| 28S-M     | 202 | 1/52    | 1/52  | 0/0       | 0/0      |
| 28S-MD    | 202 | 11/52   | 1/52  | 0/0       | 0/0      |
| 28S-L     | 385 | 2/52    | 2/52  | 0/4       | 28/32    |
| 28S-LD    | 385 | 4/52    | 0/52  | 14/28     | 57/51    |

### Primer sequences

| Gene/pair | Forward                    | Reverse                      |
|-----------|----------------------------|------------------------------|
| 16S-S     | AATTTGGTTGGGACAACC         | CTGTCAATGTGGTCTCTCTAGT       |
| 16S-SD    | AATTTGGTTGGGACRAYY         | CTGTCAATGTGRTCTCTCTARY       |
| 16S-M     | AAGGACAAGAAGACCCTA         | AATAGTTACACTAAGCCAAC         |
| 16S-MD    | AAGGACAAGAAGACCCTA         | AATAGTTACAYWRAGCCAAC         |
| 16S-L     | TTAATTGAAGGCTAGTATGA       | CAATGTGATCTCTCTAGT           |
| 16S-LD    | TTAATTGAAGGYYRGTATGA       | CAATGTGATCTCTCTARY           |
| 28S-S     | GACGCGTGCACTTTCCGG         | CAGAGGGCCTGCGACTCG           |
| 28S-SD    | GACGCGTGCAMTTTCCRG         | CAGAGGGCCTGCGACYHG           |
| 28S-M     | CTGTAGCGTTAGGGACGTCGATCGAG | CGGACTCGTGCCAGTATTTAGCCTTAGA |
| 28S-MD    | CTGTAGCGTTAGGGACGTCGRTCGAV | CGGACTCGTGCCAGTRTTTAGCSTTRGA |
| 28S-L     | GGTGAGCACGGAAGTCTCG        | CCGATTCCAGGGAGCTCGT          |
| 28S-LD    | GGTGGGCACRGAAGYCTYG        | CCGATTCCAGGGWRCTCGY          |

Binding coordinates, accession identifiers and per-species matches:
panels/Nereididae/results.csv and species.csv in the repository.
Ordinary/D pairs have identical endpoints within each design.

</div>

<div class="section page">

## Appendix A8. Hesionidae

Canonical prefix: MiPhyllo-Hesionidae-. All sequences are 5′→3′. Covered
and exact counts use the full target reference denominator. Strict and
loose background values are WC/GT.

| Gene/pair | bp  | Covered | Exact | Strict bg | Loose bg |
|-----------|-----|---------|-------|-----------|----------|
| 16S-S     | 188 | 28/53   | 2/53  | 1/7       | 1/7      |
| 16S-SD    | 188 | 37/53   | 1/53  | 1/43      | 32/44    |
| 16S-M     | 258 | 24/53   | 1/53  | 22/43     | 48/47    |
| 16S-MD    | 258 | 30/53   | 3/53  | 22/44     | 68/48    |
| 16S-L     | 362 | 28/53   | 2/53  | 2/8       | 7/27     |
| 16S-LD    | 362 | 37/53   | 3/53  | 2/60      | 50/69    |
| 28S-S     | 196 | 31/45   | 3/45  | 5/14      | 22/17    |
| 28S-SD    | 196 | 40/45   | 8/45  | 10/44     | 78/70    |
| 28S-M     | 252 | 28/45   | 6/45  | 1/4       | 17/7     |
| 28S-MD    | 252 | 35/45   | 7/45  | 2/8       | 47/11    |
| 28S-L     | 416 | 38/45   | 1/45  | 22/29     | 77/33    |
| 28S-LD    | 416 | 38/45   | 1/45  | 22/29     | 77/33    |

### Primer sequences

| Gene/pair | Forward                      | Reverse                  |
|-----------|------------------------------|--------------------------|
| 16S-S     | TGAACTAAGGTTTTAACGGC         | AAACTCTATAGGGTCTTCTTGTCT |
| 16S-SD    | TGAATTACGGYTTHAACGGC         | AAACTCAATAGGGTCTTCTTGTCY |
| 16S-M     | ACGGTTTTAACGGCCGCGGT         | TCCCCAGTCGCCCCAACTAA     |
| 16S-MD    | ACGGTTTTAACGGCYGCGGT         | TCCCCAGTCGCCCCAACYAA     |
| 16S-L     | TGAACTAAGGTTTTAACG           | ATTAGCCTGTTATCCCTA       |
| 16S-LD    | TGAATTACGGYTTHAACG           | ATTAGCCTGTTAYCCCTR       |
| 28S-S     | CACTTTCCGCGGGGAGAG           | TGTGAACAGTCCTCCCCAGT     |
| 28S-SD    | CACTTTCCGYSGRGAGAG           | TGTGAACAGTCCTYCCMRGT     |
| 28S-M     | GACGAGCGCACTTTCCGC           | GAAGCACTGTGAACAGTCCT     |
| 28S-MD    | GACGAGCGCACTTTCCRB           | GAAGCACTGTRRRCAGTCCT     |
| 28S-L     | GAGTCGGGTTGTTTGAGAATGCAGCCCA | CGAGCCATAGGACGGGCACGGC   |
| 28S-LD    | GAGTCGGGTTGTTTGRGAATGCAGCCCA | CGAGCCATAGGACGGGCACGGC   |

Binding coordinates, accession identifiers and per-species matches:
panels/Hesionidae/results.csv and species.csv in the repository.
Ordinary/D pairs have identical endpoints within each design.

</div>

<div class="section page">

## Appendix B. Conditional alternatives

Suffix 2 identifies a distinct alternative, not a mixture. Strict and
close-relative gates are retained; loose-model counts are diagnostic
rather than selection vetoes.

| Target/gene/pair   | Covered | Exact  | Strict bg | Loose bg |
|--------------------|---------|--------|-----------|----------|
| Syllidae 16S S2    | 111/180 | 5/180  | 5/11      | 547/199  |
| Syllidae 16S S2D   | 126/180 | 48/180 | 33/42     | 1007/601 |
| Syllidae 16S M2    | 56/180  | 31/180 | 31/34     | 592/88   |
| Syllidae 16S M2D   | 126/180 | 1/180  | 34/40     | 974/526  |
| Syllidae 16S L2    | 44/180  | 14/180 | 3/5       | 993/806  |
| Syllidae 16S L2D   | 91/180  | 50/180 | 12/15     | 2023/987 |
| Nereididae 28S M2  | 15/52   | 11/52  | 0/12      | 3/83     |
| Nereididae 28S M2D | 19/52   | 1/52   | 6/12      | 248/91   |

### Sequences (5′→3′)

| Pair           | Forward                 | Reverse                     |
|----------------|-------------------------|-----------------------------|
| Syllidae S2    | AAAGACAGAAAGACCCTA      | TATTAGCCTGTTATCCCTATAGT     |
| Syllidae S2D   | AAAGACAAGAMGWCCCTR      | TATTAGCCTGTTAYCCCTAHAGT     |
| Syllidae M2    | AGGTAGCATAATCATTTGC     | TATTAGCCTGTTATCCCTACAG      |
| Syllidae M2D   | AGGTAGCACAATAAYTWGY     | TATTAGCCTGTTATCCCTANAG      |
| Syllidae L2    | AAGGTAGCATAATCACTTGCCT  | CCAACATCGAGGTGCCAATCA       |
| Syllidae L2D   | AAGGTAGCATAATMAYKTGCCT  | CCAACATCGAGGTGCYAADCA       |
| Nereididae M2  | AGTAACGGCGAGTGAAGCGGGAA | CAAACAACCCGACTCTTGGGACAGCCG |
| Nereididae M2D | AGTAACGGCGAGTGAAGCGGGAR | CAAACAACCCGACTCYWAGGACAGCCR |

Full names and coordinates are in conditional/results.csv; all eight
pairs were independently re-scanned. Copyright © 2026 Fei Xia.
Repository copyright and third-party data terms apply.

</div>
