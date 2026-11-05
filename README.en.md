# MiAphro: in silico design and evaluation of mitochondrial and nuclear rRNA primers for scale-worm (Aphroditiformia) eDNA metabarcoding

**Fei Xia**

[中文](README.md) | English

## Start here

[Current 12-pair FASTA](SML/primers.fasta) · [Current results](SML/results.csv) · [Four first-round baseline pairs](baseline/results.csv)

[Preprint PDF](publications/MiAphro_preprint.pdf) · [Manuscript](publications/MiAphro_preprint.md)

## Abstract

We evaluated in silico 12 ordinary or degenerate primer pairs in three windows per marker for marine scale worms (Aphroditiformia). At 148 bp, 28S-SD covered 59/68 species with 47 paired exact matches. At 250 bp, 28S-MD covered 58/68 species with 53 exact matches. The short 16S degenerate pair recovered Eulepethidae missed by 28S but exactly matched only one species across the 16S panel.

## Introduction

Multigene phylogenetic studies provide specimen sequences across scale-worm families [1]. Nearshore Antarctic seawater surveys have also detected Polynoidae using environmental DNA [2]. We compared how primer length and degeneracy affect family coverage, exact template matching and predicted amplification in other polychaetes.

## Materials and methods

For each marker, we selected one annotated GenBank rRNA sequence per named target species, prioritizing unambiguous sequence length. After orienting and aligning sequences with MAFFT, we searched reference spans of 135–200 bp (S), 201–280 bp (M) and 281–400 bp (L). Each ordinary and D pair shares all four reference endpoints. Degenerate positions occur in the 3′ half with at most eight variants per oligo. In both concentration scenarios, minimum predicted Tm was ≥45°C, the difference between forward and reverse Tm midpoints was ≤5°C, and predicted structure Tm was <55°C.

We searched both strands for 100–500 bp products with ≥50 internal bases. Coverage allows at most two mismatches per end with compatible terminal triplets. Exact matches require zero mismatches at both ends. Denominators include all panel species, including partial references that may not span both sites. Background families were Nereididae, Phyllodocidae, Glyceridae, Nephtyidae, Hesionidae and Syllidae. Strict background counts are given for Watson–Crick (WC) and G–T-permissive pairing, each capped at ten records. Relaxed scans omit the terminal constraint and allow three mismatches per end for WC and two for GT, each capped at 50 records. Eligible windows were ranked by summed coverage and exact counts across the ordinary and D pairs.

## Results

### 16S

Target panel: 94 species. Background panel: 381 species.

![16S binding sites](figures/16S_positions.svg)

| Pair | bp | F (5′→3′) | R (5′→3′) | Covered | Exact | Strict WC/GT | Alignment |
|---|---:|---|---|---:|---:|---:|---|
| S | 138 | `ATAAACTACCTTAGGGATAACA` | `GAGATTTAACAGTTGAACAAAC` | 58/94 | 1/94 | 0/0 | [S](figures/template_matching/MiAphro-16S-S.svg) |
| SD | 138 | `ATAAACTACCTYRGGGATAACA` | `GAGATTTAACAGTTGAACRAAC` | 62/94 | 1/94 | 0/0 | [SD](figures/template_matching/MiAphro-16S-SD.svg) |
| M | 210 | `AAGACAAGAAGACCCTATAG` | `ACCATCCTATCAATAAGAGCTC` | 28/94 | 24/94 | 0/0 | [M](figures/template_matching/MiAphro-16S-M.svg) |
| MD | 210 | `AAGACAAGAAGACCCYRTWG` | `ACCATCCTATCAATWWGARCTC` | 44/94 | 24/94 | 2/9 | [MD](figures/template_matching/MiAphro-16S-MD.svg) |
| L | 329 | `GTAGCATAATCATTTGCCTT` | `ACCATCCTATCAATAAGAGC` | 28/94 | 24/94 | 0/0 | [L](figures/template_matching/MiAphro-16S-L.svg) |
| LD | 329 | `GTAGCATAATMAYTTGCCTY` | `ACCATCCTATCAATWWGARC` | 41/94 | 24/94 | 0/7 | [LD](figures/template_matching/MiAphro-16S-LD.svg) |

Blue indicates compatible bases, orange mismatches and grey gaps or unknown bases. Boxes enclose terminal triplets. The excerpt shows one species per family; the table links show every target species.

![16S SD](figures/template_matching/MiAphro-16S-SD_preview.svg)

### 28S

Target panel: 68 species. Background panel: 214 species.

![28S binding sites](figures/28S_positions.svg)

| Pair | bp | F (5′→3′) | R (5′→3′) | Covered | Exact | Strict WC/GT | Alignment |
|---|---:|---|---|---:|---:|---:|---|
| S | 148 | `CTCAGAGGCAAACTGATGGGAC` | `GCTCTCCTCGTGGAAAGGGC` | 54/68 | 29/68 | 0/0 | [S](figures/template_matching/MiAphro-28S-S.svg) |
| SD | 148 | `CTCAGAGGCAAACYGATGGGAC` | `GCTCTCCTCGYRGAAAGKGC` | 59/68 | 47/68 | 2/2 | [SD](figures/template_matching/MiAphro-28S-SD.svg) |
| M | 250 | `AATCCCCTGTCGCTGCGGCA` | `TCCCTCACGGTACTTGTCCGCT` | 56/68 | 50/68 | 0/0 | [M](figures/template_matching/MiAphro-28S-M.svg) |
| MD | 250 | `AATCCCCTGTCGYTGCGGCA` | `TCCCTCACGGTACTTGTCCKCT` | 58/68 | 53/68 | 1/2 | [MD](figures/template_matching/MiAphro-28S-MD.svg) |
| L | 368 | `CTCAGAGGCAAACTGATGGGAC` | `CGGAACACACTGTAAACAGTCC` | 58/68 | 22/68 | 0/0 | [L](figures/template_matching/MiAphro-28S-L.svg) |
| LD | 368 | `CTCAGAGGCAAACYGATGGGAC` | `CGGAACACACTGTRAACAGTCC` | 59/68 | 33/68 | 0/1 | [LD](figures/template_matching/MiAphro-28S-LD.svg) |

Blue indicates compatible bases, orange mismatches and grey gaps or unknown bases. Boxes enclose terminal triplets. The excerpt shows one species per family; the table links show every target species.

![28S MD](figures/template_matching/MiAphro-28S-MD_preview.svg)

### Family-level coverage

![Family-level species coverage](figures/family_coverage.svg)

Bars show species coverage within each family; labels give covered/total species and percentages. Panels use marker-specific references and include partial sequences. Ordinary and D pairs are shown separately.

## Discussion

28S-SD offers a short product and broad predicted coverage. 28S-MD provides the most exact matches. The long 28S-LD pair matches as many species as SD but yields only 33 exact matches. No selected pair detected Aphroditidae on these panels. All 28S pairs missed Eulepethidae, whereas 16S-SD covered all three representatives. The short 16S pair's broad mismatch-tolerant coverage coincides with low exact matching.

All four Aphroditidae 16S records span all three windows. All three 28S records in that family span S and M, but only two span L. Both Eulepethidae 28S records span all windows. Thus one Aphroditidae record lacks the full 28S L interval; fragment length alone cannot explain the other misses in these families.

## Data and reproducibility

We checked the 12 current pairs and four first-round pairs. Among complete, indel-free templates for each end, we flagged sites where the most common base was absent from the oligo and more frequent than any included base. We tested single and combined replacements. No current 28S pair had such a site. For 16S-SD, reverse C10T would raise coverage from 62 to 69 but lose five previously covered species and add 79 background record–model hits. Other tested 16S substitutions also gained background hits. We retained the published sequences.

[FASTA](SML/primers.fasta) · [SML results](SML/results.csv) · [Family coverage](SML/families.csv) · [Species records](SML/species.csv) · [First-round baseline](baseline/results.csv) · [Pair audit](audit/pair_review.csv) · [Position counts](audit/position_support.csv) · [Consensus trials](audit/consensus_trials.csv)

With Biopython and primer3-py installed, run `python3 scripts/finish_sml.py` to verify the frozen per-record scans and `python3 scripts/review.py` to regenerate base matching and substitution tests.

Run `python3 scripts/check_release.py` to check the documentation, reference panels, FASTA exports and matching figures together.

## Citation

Xia, F. (2026). *MiAphro: in silico design and evaluation of mitochondrial and nuclear rRNA primers for scale-worm (Aphroditiformia) eDNA metabarcoding*. [MiAphro](https://github.com/GUIBA-EX/MiAphro).

## References

1. Norlinder et al. (2012). Phylogeny of scale-worms (Aphroditiformia, Annelida), assessed from 18SrRNA, 28SrRNA, 16SrRNA, mitochondrial cytochrome c oxidase subunit I (COI), and morphology. [PubMed](https://pubmed.ncbi.nlm.nih.gov/22789762/).
2. Environmental DNA metabarcoding for monitoring metazoan biodiversity in Antarctic nearshore ecosystems (2021). *PeerJ* 9, e12458. [DOI](https://doi.org/10.7717/peerj.12458).

## Rebuild the exhaustive scanner

From the repository root, run `c++ -std=c++17 -O2 scripts/raw_scan.cpp -o scripts/scan` before repeating window discovery with MAFFT, Biopython and primer3-py. The per-record verification scripts use the archived scan tables directly.
