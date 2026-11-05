# MiAphro: in silico design and evaluation of mitochondrial and nuclear rRNA primers for scale-worm (Aphroditiformia) eDNA metabarcoding

**Fei Xia**

[English](README.en.md) | 中文

## 使用入口

[现行 12 对引物 FASTA](SML/primers.fasta) · [现行结果](SML/results.csv) · [首轮 4 对基线](baseline/results.csv)

[预印本 PDF](publications/MiAphro_preprint.pdf) · [论文原稿](publications/MiAphro_preprint.md)

## 摘要

本研究通过序列比对和计算扩增，面向海洋鳞虫类（Aphroditiformia），比较线粒体 16S 和核 28S 各三个扩增窗口，共 12 对普通或简并引物。28S-SD 在 148 bp 的参考跨度下覆盖 59/68 种，47 种双端完全匹配。250 bp 的 28S-MD 覆盖 58/68 种，双端完全匹配增至 53 种。16S-SD 补充了 28S 漏检的 Eulepethidae，但在整个 16S 面板中只有 1 种双端完全匹配。

## 引言

鳞虫的多基因系统研究为跨科引物设计提供了标本参考序列 [1]。南极近岸海水的环境 DNA 调查已检出 Polynoidae [2]。本项目比较各长度窗口能覆盖哪些科、两端能否完全匹配，以及能否限制其他多毛类的预测扩增。

## 材料与方法

每个基因从 GenBank 具名目标物种中各保留一条带对应 rRNA 注释的序列，以明确碱基数最多者优先。MAFFT 定向比对后，在参考片段上分别搜索 S（135–200 bp）、M（201–280 bp）和 L（281–400 bp）。普通与 D 的 F、R 四个结合端点完全相同。每端最多 8 种变体，简并位点限于 3′ 半段。两种浓度情景下预测最低 Tm ≥45°C，正反向 Tm 中点之差 ≤5°C，预测结构 Tm <55°C。

在原始序列双链搜索 100–500 bp、内部至少 50 bp 的产物。覆盖要求每端至多两个错配且各自 3′ 末三位兼容；完全匹配要求两端零错配。分母为全部目标参考物种，部分序列可能没有跨越整个窗口。背景包括 Nereididae、Phyllodocidae、Glyceridae、Nephtyidae、Hesionidae 和 Syllidae。表中严格背景分别报告 Watson–Crick（WC）及允许 G–T 配对的记录数，两者各须 ≤10。宽松检查取消末端限制，WC 每端允许三个错配、GT 两个，各须 ≤50。在同档通过门槛的窗口内，按普通和 D 的覆盖与完全匹配数之和排序。

## 结果

### 16S

目标面板 94 种，背景面板 381 种。

![16S binding sites](figures/16S_positions.svg)

| 引物 | bp | F（5′→3′） | R（5′→3′） | 覆盖 | 完全匹配 | 严格 WC/GT | 匹配图 |
|---|---:|---|---|---:|---:|---:|---|
| S | 138 | `ATAAACTACCTTAGGGATAACA` | `GAGATTTAACAGTTGAACAAAC` | 58/94 | 1/94 | 0/0 | [S](figures/template_matching/MiAphro-16S-S.svg) |
| SD | 138 | `ATAAACTACCTYRGGGATAACA` | `GAGATTTAACAGTTGAACRAAC` | 62/94 | 1/94 | 0/0 | [SD](figures/template_matching/MiAphro-16S-SD.svg) |
| M | 210 | `AAGACAAGAAGACCCTATAG` | `ACCATCCTATCAATAAGAGCTC` | 28/94 | 24/94 | 0/0 | [M](figures/template_matching/MiAphro-16S-M.svg) |
| MD | 210 | `AAGACAAGAAGACCCYRTWG` | `ACCATCCTATCAATWWGARCTC` | 44/94 | 24/94 | 2/9 | [MD](figures/template_matching/MiAphro-16S-MD.svg) |
| L | 329 | `GTAGCATAATCATTTGCCTT` | `ACCATCCTATCAATAAGAGC` | 28/94 | 24/94 | 0/0 | [L](figures/template_matching/MiAphro-16S-L.svg) |
| LD | 329 | `GTAGCATAATMAYTTGCCTY` | `ACCATCCTATCAATWWGARC` | 41/94 | 24/94 | 0/7 | [LD](figures/template_matching/MiAphro-16S-LD.svg) |

蓝色表示兼容，橙色为错配，灰色为缺口或未知碱基；框出 3′ 末三位。示例每科取一条，表内链接展示全部目标物种。

![16S SD](figures/template_matching/MiAphro-16S-SD_preview.svg)

### 28S

目标面板 68 种，背景面板 214 种。

![28S binding sites](figures/28S_positions.svg)

| 引物 | bp | F（5′→3′） | R（5′→3′） | 覆盖 | 完全匹配 | 严格 WC/GT | 匹配图 |
|---|---:|---|---|---:|---:|---:|---|
| S | 148 | `CTCAGAGGCAAACTGATGGGAC` | `GCTCTCCTCGTGGAAAGGGC` | 54/68 | 29/68 | 0/0 | [S](figures/template_matching/MiAphro-28S-S.svg) |
| SD | 148 | `CTCAGAGGCAAACYGATGGGAC` | `GCTCTCCTCGYRGAAAGKGC` | 59/68 | 47/68 | 2/2 | [SD](figures/template_matching/MiAphro-28S-SD.svg) |
| M | 250 | `AATCCCCTGTCGCTGCGGCA` | `TCCCTCACGGTACTTGTCCGCT` | 56/68 | 50/68 | 0/0 | [M](figures/template_matching/MiAphro-28S-M.svg) |
| MD | 250 | `AATCCCCTGTCGYTGCGGCA` | `TCCCTCACGGTACTTGTCCKCT` | 58/68 | 53/68 | 1/2 | [MD](figures/template_matching/MiAphro-28S-MD.svg) |
| L | 368 | `CTCAGAGGCAAACTGATGGGAC` | `CGGAACACACTGTAAACAGTCC` | 58/68 | 22/68 | 0/0 | [L](figures/template_matching/MiAphro-28S-L.svg) |
| LD | 368 | `CTCAGAGGCAAACYGATGGGAC` | `CGGAACACACTGTRAACAGTCC` | 59/68 | 33/68 | 0/1 | [LD](figures/template_matching/MiAphro-28S-LD.svg) |

蓝色表示兼容，橙色为错配，灰色为缺口或未知碱基；框出 3′ 末三位。示例每科取一条，表内链接展示全部目标物种。

![28S MD](figures/template_matching/MiAphro-28S-MD_preview.svg)

### 各科覆盖

![Family coverage](figures/family_coverage.svg)

热图按 WoRMS 核验后的科级归属汇总，每格显示覆盖率及覆盖种数/参考种数；灰色 NA 表示缺少参考序列。Pholoidae 和 Pisionidae 的参考属归入 Sigalionidae。原始 GenBank 分类保存在源数据中。

### 参考库分类覆盖

![Reference taxonomic sampling](figures/reference_taxonomy.svg)

参考库图分别统计 16S、28S 及两者合并后的去重属数和物种数。属名和科归属经 WoRMS 核验；物种保留源记录名称，未进行物种级异名归并。这里展示已采样数量，不代表占 WoRMS 全部有效属的比例。

[WoRMS](https://www.marinespecies.org/aphia.php?p=taxdetails&id=927) · 2026-11-05 · [Taxonomic mapping](taxonomy/references.csv) · [Genus counts and names](taxonomy/reference_coverage.csv)

## 讨论

28S-SD 适合优先考虑短片段和覆盖范围，28S-MD 则得到最多的双端完全匹配。28S-LD 与短版覆盖数相同，但完全匹配降至 33 种。在当前面板中，Aphroditidae 在两个基因的全部入选方案中均未命中。28S 全部漏掉 Eulepethidae，16S-SD 则覆盖其中 3/3 种。16S 短版的高容错覆盖和低完全匹配须同时解读。

进一步核对参考片段两端：4 条 Aphroditidae 的 16S 记录均跨越三个窗口；其 3 条 28S 记录均跨越 S 和 M，但仅 2 条跨越 L。Eulepethidae 的 2 条 28S 记录均跨越三个窗口。因此，28S 长档有一条 Aphroditidae 记录的未命中需结合片段缺失解读，其余这些未命中不能仅用片段过短解释。

## 数据与复现

复核涵盖现行 12 对和首轮 4 对。按每端完整、无插缺的模板统计逐位碱基；若最常见碱基未被引物包含且频数高于任一已包含碱基，则重算单点及组合替换。现行 28S 未见这一类偏差。16S-SD 的 R10 C→T 可将覆盖由 62 提高至 69 种，但丢失原先覆盖的 5 种，并新增 79 个背景「记录×模型」命中。其他 16S 修改也会增加背景命中，本轮保留原序列。

[FASTA](SML/primers.fasta) · [SML results](SML/results.csv) · [Family coverage](SML/families.csv) · [Species records](SML/species.csv) · [First-round baseline](baseline/results.csv) · [Pair audit](audit/pair_review.csv) · [Position counts](audit/position_support.csv) · [Consensus trials](audit/consensus_trials.csv)

使用 Biopython 与 primer3-py，运行 `python3 scripts/finish_sml.py` 复核冻结的逐记录搜索结果，再运行 `python3 scripts/review.py` 重算逐位匹配与替换试验。

运行 `python3 scripts/check_release.py` 可检查文档、面板、FASTA 及匹配图的一致性。

## 引用 / Citation

Xia, F. (2026). *MiAphro: in silico design and evaluation of mitochondrial and nuclear rRNA primers for scale-worm (Aphroditiformia) eDNA metabarcoding*. [MiAphro](https://github.com/GUIBA-EX/MiAphro).

## 参考文献

1. Norlinder et al. (2012). Phylogeny of scale-worms (Aphroditiformia, Annelida), assessed from 18SrRNA, 28SrRNA, 16SrRNA, mitochondrial cytochrome c oxidase subunit I (COI), and morphology. [PubMed](https://pubmed.ncbi.nlm.nih.gov/22789762/).
2. Environmental DNA metabarcoding for monitoring metazoan biodiversity in Antarctic nearshore ecosystems (2021). *PeerJ* 9, e12458. [DOI](https://doi.org/10.7717/peerj.12458).

## 重建搜索程序

重新搜索窗口前，在仓库根目录运行 `c++ -std=c++17 -O2 scripts/raw_scan.cpp -o scripts/scan`，并安装 MAFFT、Biopython 和 primer3-py。逐记录复核可直接使用随仓库提供的扫描表。
