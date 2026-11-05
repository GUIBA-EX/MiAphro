# MiAphro: in silico design and evaluation of mitochondrial and nuclear rRNA primers for scale-worm (Aphroditiformia) eDNA metabarcoding

**Fei Xia**

Computational preprint · 5 November 2026

## Abstract

We evaluated in silico 12 ordinary or degenerate primer pairs in three windows per marker for marine scale worms (Aphroditiformia). At 148 bp, 28S-SD covered 59/68 species with 47 paired exact matches. At 250 bp, 28S-MD covered 58/68 species with 53 exact matches. The short 16S degenerate pair recovered Eulepethidae missed by 28S but exactly matched only one species across the 16S panel. The 28S short and medium D pairs each matched two or fewer background species under either strict pairing model. All selected assays missed Aphroditidae in the frozen panels.

**Keywords:** Aphroditiformia, environmental DNA, mitochondrial 16S, nuclear 28S, primer degeneracy

## Introduction

Multigene phylogenetic studies provide specimen sequences across scale-worm families [1]. Nearshore Antarctic seawater surveys have also detected Polynoidae using environmental DNA [2]. We compared how primer length and degeneracy affect family coverage, exact template matching and predicted amplification in other polychaetes. We asked whether adding degeneracy at fixed coordinates increased paired exact matches, and whether either marker recovered families missed by the other.

## Materials and methods

### Data acquisition, alignment and candidate generation

For each marker, we selected one annotated GenBank rRNA sequence per named target species, prioritizing unambiguous sequence length. After orienting and aligning sequences with MAFFT, we searched reference spans of 135–200 bp (S), 201–280 bp (M) and 281–400 bp (L). Each ordinary and D pair shares all four reference endpoints. Degenerate positions occur in the 3′ half with at most eight variants per oligo. In both concentration scenarios, minimum predicted Tm was ≥45°C, the difference between forward and reverse Tm midpoints was ≤5°C, and predicted structure Tm was <55°C.

### Matching models and selection criteria

We searched both strands for 100–500 bp products with ≥50 internal bases. Coverage allows at most two mismatches per end with compatible terminal triplets. Exact matches require zero mismatches at both ends. Denominators include all panel species, including partial references that may not span both sites. Background families were Nereididae, Phyllodocidae, Glyceridae, Nephtyidae, Hesionidae and Syllidae. Strict background counts are given for Watson–Crick (WC) and G–T-permissive pairing, each capped at ten records. Relaxed scans omit the terminal constraint and allow three mismatches per end for WC and two for GT, each capped at 50 records. Eligible windows were ranked by summed coverage and exact counts across the ordinary and D pairs.

### Reference panels and thermodynamic evaluation

The 16S target panel included 69 Polynoidae among 94 species; 28S included 41 among 68. Background panels comprised 381 and 214 species, respectively. Coordinates are 1-based and inclusive on the ungapped reference fragments MG905038.1 (16S) and DQ790064.1 (28S). Binding-site arrows point in synthesis orientation. Candidate oligos were 20 or 22 bases long. D variants were assessed at 50 nM per variant and 200 nM per direction divided among variants.

### Independent implementation and consensus checks

We screened 3,793 distinct marker/window combinations across the six size classes, including first-round windows. A second implementation checked raw paired hits on both strands. Family counts were tallied from species outcomes. For the consensus audit, we counted bases on complete, indel-free binding sites and re-evaluated substitutions when the most frequent template base was absent from the oligo. Reported background gains in this audit count record–model combinations.

## Results

### Mitochondrial 16S coverage and exact matching

The 138 bp SD pair covered 62 of 94 species, but only one species matched both ends exactly. The medium D pair covered 44 species with 24 exact matches. It gained two covered species over the first-round 16S D pair at the same exact-match count.

![Binding sites](../figures/16S_positions.svg)

**Figure 1.** Primer binding sites and inclusive reference spans. Each ordinary/D pair shares both sites.



### Nuclear 28S coverage and exact matching

The 148 bp SD and 368 bp LD pairs each covered 59 of 68 species. The short pair had 47 exact matches, compared with 33 for LD. The 250 bp MD pair yielded the highest exact count (53), with one fewer covered species than SD.

![Binding sites](../figures/28S_positions.svg)

**Figure 2.** Primer binding sites and inclusive reference spans. Each ordinary/D pair shares both sites.



### Consensus audit and primer–template comparisons

We checked the 12 current pairs and four first-round pairs. Among complete, indel-free templates for each end, we flagged sites where the most common base was absent from the oligo and more frequent than any included base. We tested single and combined replacements. No current 28S pair had such a site. For 16S-SD, reverse C10T would raise coverage from 62 to 69 but lose five previously covered species and add 79 background record–model hits. Other tested 16S substitutions also gained background hits. We retained the published sequences.

![16S SD](../figures/template_matching/MiAphro-16S-SD_preview.svg)

**Figure 3.** 16S-SD binding sites, with one representative per family. Blue indicates compatibility, orange mismatch and grey a gap or unknown base. The terminal three bases are boxed. The hit column summarizes the complete-sequence paired scan.

![28S MD](../figures/template_matching/MiAphro-28S-MD_preview.svg)

**Figure 4.** 28S-MD binding sites, with one representative per family. Blue indicates compatibility, orange mismatch and grey a gap or unknown base. The terminal three bases are boxed. The hit column summarizes the complete-sequence paired scan.

### Family-level coverage

![Family coverage](../figures/family_coverage.svg)

**Figure 5.** The heatmap uses WoRMS-verified family assignments. Each cell shows coverage percentage and covered/total reference species; grey NA indicates no references. Genera recorded under Pholoidae and Pisionidae are assigned to Sigalionidae. Original GenBank labels remain in the source tables.

### Reference taxonomic sampling

![Reference taxonomic sampling](../figures/reference_taxonomy.svg)

**Figure 6.** Reference sampling is summarized as unique genera and species for 16S, 28S and their union. Genus names and family assignments were checked against WoRMS; species retain source-record names without species-level synonym consolidation. Counts describe sampled taxa, not proportions of all accepted WoRMS genera.

[WoRMS](https://www.marinespecies.org/aphia.php?p=taxdetails&id=927) · 2026-11-05 · [Taxonomic mapping](../taxonomy/references.csv) · [Genus counts and names](../taxonomy/reference_coverage.csv)

## Discussion

28S-SD offers a short product and broad predicted coverage. 28S-MD provides the most exact matches. The long 28S-LD pair matches as many species as SD but yields only 33 exact matches. No selected pair detected Aphroditidae on these panels. All 28S pairs missed Eulepethidae, whereas 16S-SD covered all three representatives. The short 16S pair's broad mismatch-tolerant coverage coincides with low exact matching.

All four Aphroditidae 16S records span all three windows. All three 28S records in that family span S and M, but only two span L. Both Eulepethidae 28S records span all windows. Thus one Aphroditidae record lacks the full 28S L interval; fragment length alone cannot explain the other misses in these families.

## Conclusions

The short and medium 28S D pairs provide the strongest balance of predicted coverage and exact matching in the curated panels. The 16S short pair adds Eulepethidae coverage. Further design should prioritize the shared Aphroditidae gap.

## Data availability and reproducibility

All primer sequences, target alignments, background references, per-species predictions and review scripts are available in the [MiAphro repository](https://github.com/GUIBA-EX/MiAphro). The 12 S/M/L pairs are the current panel. Four first-round pairs are archived separately.

## References

1. Norlinder et al. (2012). Phylogeny of scale-worms (Aphroditiformia, Annelida), assessed from 18SrRNA, 28SrRNA, 16SrRNA, mitochondrial cytochrome c oxidase subunit I (COI), and morphology. [PubMed](https://pubmed.ncbi.nlm.nih.gov/22789762/).
2. Environmental DNA metabarcoding for monitoring metazoan biodiversity in Antarctic nearshore ecosystems (2021). *PeerJ* 9, e12458. [DOI](https://doi.org/10.7717/peerj.12458).

## Table 1. Mitochondrial 16S paired matching

Species-level target outcomes and background record counts. Target denominator: 94; background: 381.

| Pair | Span (bp) | Covered | Exact | Strict WC / GT | Relaxed WC / GT |
|---|---|---|---|---|---|
| S | 138 | 58 | 1 | 0 / 0 | 23 / 0 |
| SD | 138 | 62 | 1 | 0 / 0 | 25 / 0 |
| M | 210 | 28 | 24 | 0 / 0 | 0 / 0 |
| MD | 210 | 44 | 24 | 2 / 9 | 12 / 9 |
| L | 329 | 28 | 24 | 0 / 0 | 0 / 0 |
| LD | 329 | 41 | 24 | 0 / 7 | 7 / 7 |

## Table 2. Nuclear 28S paired matching

Species-level target outcomes and background record counts. Target denominator: 68; background: 214.

| Pair | Span (bp) | Covered | Exact | Strict WC / GT | Relaxed WC / GT |
|---|---|---|---|---|---|
| S | 148 | 54 | 29 | 0 / 0 | 16 / 2 |
| SD | 148 | 59 | 47 | 2 / 2 | 19 / 2 |
| M | 250 | 56 | 50 | 0 / 0 | 16 / 0 |
| MD | 250 | 58 | 53 | 1 / 2 | 17 / 2 |
| L | 368 | 58 | 22 | 0 / 0 | 15 / 0 |
| LD | 368 | 59 | 33 | 0 / 1 | 18 / 1 |

## Table 3. Primer sequences

**Table 3.** Both oligos are shown in synthesis orientation (5′–3′).

| Pair | Forward | Reverse |
|---|---|---|
| MiAphro-16S-S | ATAAACTACCTTAGGGATAACA | GAGATTTAACAGTTGAACAAAC |
| MiAphro-16S-SD | ATAAACTACCTYRGGGATAACA | GAGATTTAACAGTTGAACRAAC |
| MiAphro-16S-M | AAGACAAGAAGACCCTATAG | ACCATCCTATCAATAAGAGCTC |
| MiAphro-16S-MD | AAGACAAGAAGACCCYRTWG | ACCATCCTATCAATWWGARCTC |
| MiAphro-16S-L | GTAGCATAATCATTTGCCTT | ACCATCCTATCAATAAGAGC |
| MiAphro-16S-LD | GTAGCATAATMAYTTGCCTY | ACCATCCTATCAATWWGARC |
| MiAphro-28S-S | CTCAGAGGCAAACTGATGGGAC | GCTCTCCTCGTGGAAAGGGC |
| MiAphro-28S-SD | CTCAGAGGCAAACYGATGGGAC | GCTCTCCTCGYRGAAAGKGC |
| MiAphro-28S-M | AATCCCCTGTCGCTGCGGCA | TCCCTCACGGTACTTGTCCGCT |
| MiAphro-28S-MD | AATCCCCTGTCGYTGCGGCA | TCCCTCACGGTACTTGTCCKCT |
| MiAphro-28S-L | CTCAGAGGCAAACTGATGGGAC | CGGAACACACTGTAAACAGTCC |
| MiAphro-28S-LD | CTCAGAGGCAAACYGATGGGAC | CGGAACACACTGTRAACAGTCC |
