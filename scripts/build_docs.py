from pathlib import Path
import csv,json
from collections import Counter
R=Path(__file__).resolve().parents[1]
read=lambda p:list(csv.DictReader(p.open()))
rows=read(R/'catalog.csv');alt=read(R/'conditional/results.csv');panels=read(R/'panels.csv');taxa=list(dict.fromkeys(x['target_taxon'] for x in panels))
notes={
'Aphroditiformia':('保留原鳞虫研究、碱基匹配图、WoRMS属级覆盖与预印本。','The original scale-worm study, base-matching figures, WoRMS genus coverage and preprint are retained.'),
'Phyllodocida':('目级覆盖仍低，保留为广谱搜索结果。','Order-wide coverage remains low; these candidates document the broad-target search.'),
'Glyceriformia':('参考未包含Lacydoniidae；16S覆盖与完全匹配之间存在明显差距。','Lacydoniidae is absent from the reference panel; 16S coverage substantially exceeds exact matching.'),
'Nereidiformia':('亚目级覆盖仍低；另提供三个科的独立设计。','Suborder-wide coverage remains low; three constituent families have separate designs.'),
'Phyllodociformia':('当前参考仅代表Phyllodocidae。','The current references represent Phyllodocidae only.'),
'Syllidae':('优先关注16S。28S的S/M覆盖很低，L无通过候选。','Prioritize examination of 16S. The 28S S/M candidates have very low coverage; no L pair passed.'),
'Nereididae':('16S覆盖较高；28S以S为主要候选，M/L覆盖仍有限。','16S has relatively high panel coverage; S is the principal 28S candidate, while M/L remain limited.'),
'Hesionidae':('28S S覆盖40/45种；完整长度比较见本单元结果。','28S S covers 40/45 reference species; all length comparisons are retained.')}
for en in [False,True]:
 suffix='README.en.md' if en else 'README.md';other='README.md' if en else 'README.en.md'
 title='MiPhyllo: taxon-stratified 16S and 28S PCR primer candidates for environmental DNA metabarcoding of Phyllodocida'
 lines=['# '+title,'','**Fei Xia**','', '[Integrated preprint / 整合版预印本](publications/MiPhyllo_preprint.pdf) · [Manuscript / 原稿](publications/MiPhyllo_preprint.md)', '', 'MiPhyllo names the expanded collection; the MiAphro repository URL and sequence identifiers are retained. / 整合版命名MiPhyllo，保留MiAphro仓库地址与序列编号。', '',f'[中文]({other}) | English' if en else f'[English]({other}) | 中文','', '## Overview' if en else '## 项目概览','',
 'MiPhyllo integrates eight nested target units within Phyllodocida: the order itself, four suborders, and three families. It contains 94 current primer pairs (188 oligos) and 8 distinct conditional alternatives (16 oligos). Each unit retains its own reference panel, background panel and selection criteria; this is a collection of computational designs, not one universal assay for all polychaetes.' if en else 'MiPhyllo整合Phyllodocida内8个相互嵌套的设计单元：目级、4个亚目及3个科。现行目录包含94对引物（188条寡核苷酸），另有8对不同的条件备选（16条寡核苷酸）。各设计单元保留自己的参考库、背景库和筛选条件；这些是分级设计的计算候选，并非一套覆盖全部多毛类的通用实验方案。','',
 '[Current FASTA](primers.fasta) · [Catalogue](catalog.csv) · [Reference panels](panels.csv) · [Conditional alternatives](conditional/README.en.md) · [Methods](docs/METHODS.en.md)' if en else '[现行序列](primers.fasta) · [完整目录](catalog.csv) · [参考库清单](panels.csv) · [条件备选](conditional/README.md) · [方法与复现](docs/METHODS.md)','',
 '## Target units' if en else '## 设计单元','',
 '|Target|Rank|16S references|28S references|Results|' if en else '|目标类群|层级|16S参考物种|28S参考物种|结果入口|','|---|---|---:|---:|---|']
 for taxon in taxa:
  ps={x['gene']:x for x in panels if x['target_taxon']==taxon};rank=ps['16S']['rank'];rank=rank if en else {'order':'目','suborder':'亚目','family':'科'}[rank]
  lines.append(f"|{taxon}|{rank}|{ps['16S']['target_species']}|{ps['28S']['target_species']}|[{taxon}](panels/{taxon}/{suffix})|")
 lines+=['','```mermaid','graph TD','  P[Phyllodocida] --> A[Aphroditiformia]','  P --> G[Glyceriformia]','  P --> N[Nereidiformia]','  P --> F[Phyllodociformia]','  N --> S[Syllidae]','  N --> R[Nereididae]','  N --> H[Hesionidae]','```','',
 'The diagram shows included design units, not a complete taxonomic classification. Nested panels contain overlapping species and must not be summed.' if en else '示意图仅列出本项目的设计单元，不是完整分类系统。上下级面板包含重叠物种，不能直接相加。','',
 '## Names and profiles' if en else '## 命名与方案','',
 '`MiAphro-{taxon}-{gene}-{length}[2][D]_{F|R}`. S/M/L denote short, medium and long designs. D marks degeneracy in at least one oligo of the pair; each oligo has degeneracy ≤8. The suffix 2 identifies a distinct conditional alternative, not a two-pair mixture. Examples: `MiAphro-Syllidae-16S-S_F`, `MiAphro-Syllidae-16S-SD_R`, and `MiAphro-Syllidae-16S-S2D_R`. Ordinary/D counterparts share all four endpoints within a profile. Different profiles or taxa may use different windows and references.' if en else '`MiAphro-{目标类群}-{基因}-{长度}[2][D]_{F|R}`。S/M/L表示短、中、长设计；D表示该对引物至少一条含简并，每条简并度≤8。数字2表示另一套条件候选，不表示两对引物混合。例：`MiAphro-Syllidae-16S-S_F`、`MiAphro-Syllidae-16S-SD_R`、`MiAphro-Syllidae-16S-S2D_R`。同一方案内普通/D的四个结合端点一致；不同方案、类群可能使用不同窗口和参考序列。','',
 '## Results' if en else '## 结果概览','',
 '![Within-panel D coverage](figures/catalog_coverage.svg)','',
 'The figure shows D coverage within each frozen reference panel. Panel sizes and background composition differ; these are not common-panel comparisons.' if en else '图示为各自冻结参考库内的D方案覆盖率。不同设计单元的分母和背景组成不同，不能解释为共同面板上的直接性能比较。','']
 for taxon in taxa:lines.append('- **'+taxon+'**: '+notes[taxon][int(en)])
 lines+=['',
 'The latest Nereididae 28S SD covers 29/52 species with 8 exact matches. A conditional Syllidae 16S S2D covers 126/180 with 48 exact matches: strict background WC/GT hits are 33/42 of 3,675, but loose-model hits are 1,007/601. Conditional candidates remain separate from the current FASTA.' if en else '最新Nereididae 28S SD覆盖29/52种，完全匹配8种。条件备选Syllidae 16S S2D覆盖126/180种、完全匹配48种；在3675种背景中，严格WC/GT命中33/42，宽松模型命中1007/601。条件备选与现行FASTA分开保存。','',
 '## Reproducibility and evidence' if en else '## 数据、复现与证据','',
 '[Species-level results](species.csv) · [Naming provenance](provenance.csv) · [Independent re-scan](audit/catalog_rescan.json) · [Checksums](SHA256SUMS)' if en else '[逐物种结果](species.csv) · [来源与命名对应](provenance.csv) · [独立重扫检查](audit/catalog_rescan.json) · [文件校验值](SHA256SUMS)','',
 '```bash','python3 -m pip install -r requirements.txt','python3 scripts/build_catalog.py','python3 scripts/build_docs.py','python3 scripts/check_release.py','# Optional full independent re-scan / 可选：完整独立重扫','python3 scripts/recheck_catalog.py','```','',
 'Frozen references and alignments are under `panels/`; compact source analyses are under `analysis/`. Large intermediate hit matrices and caches are regenerated rather than duplicated. Selected results were re-scanned with independent C++ and Python implementations; this integration adds no wet-lab validation.' if en else '冻结参考序列和比对位于`panels/`，精简后的源分析位于`analysis/`。大型中间命中矩阵与缓存可重建，不重复收录。入选结果已使用独立C++和Python实现重新扫描；本次整合未添加实验验证。','',
 '## Aphroditiformia study and preprint' if en else '## Aphroditiformia原研究与预印本','',
 '[Study and bibliography](studies/Aphroditiformia/'+suffix+') · [Preprint PDF](studies/Aphroditiformia/publications/MiAphro_preprint.pdf) · [Manuscript](studies/Aphroditiformia/publications/MiAphro_preprint.md)','',
 'The retained preprint reports Aphroditiformia only. The additional order-, suborder- and family-level designs are documented in this catalogue and their individual pages.' if en else '保留的预印本仅报告Aphroditiformia。新增目级、亚目级和科级设计由当前目录及各单元页面报告，不将其归入原预印本结论。','',
 '## Citation and copyright' if en else '## 引用与著作权','',f'Xia, F. (2026). *{title}*. [MiAphro](https://github.com/GUIBA-EX/MiAphro).','', '[Citation metadata](CITATION.cff) · [Copyright and permissions](COPYRIGHT.md)','']
 (R/suffix).write_text('\n'.join(lines))
 for taxon in taxa:
  rs=[x for x in rows if x['target_taxon']==taxon];lines=[f'# MiAphro · {taxon}','', '[Catalogue](../../'+suffix+') | [FASTA](primers.fasta) | [CSV](results.csv) | [Species](species.csv)','',notes[taxon][int(en)],'']
  for gene in ['16S','28S']:
   p=next(x for x in panels if x['target_taxon']==taxon and x['gene']==gene);metadata=read(R/'panels'/taxon/'references'/gene/'target.csv');families=Counter(x['family'] for x in metadata)
   lines+=['## '+gene,'',f"Target / background reference species: **{p['target_species']} / {p['background_species']}**" if en else f"目标 / 背景参考物种：**{p['target_species']} / {p['background_species']}**",'', 'Families: '+', '.join(f'{k} ({v})' for k,v in sorted(families.items())),'',f'[Target FASTA](references/{gene}/target.fa) · [Metadata](references/{gene}/target.csv) · [Alignment](references/{gene}/aligned.fa) · [Background FASTA](references/{gene}/background.fa) · [Background metadata](references/{gene}/background.csv)','',
   '|Pair|bp|F (5′→3′)|R (5′→3′)|Covered|Exact|Strict WC/GT|Loose WC/GT|' if en else '|引物对|bp|F（5′→3′）|R（5′→3′）|覆盖|完全匹配|严格WC/GT|宽松WC/GT|','|---|---:|---|---|---:|---:|---:|---:|']
   for x in rs:
    if x['gene']==gene:lines.append(f"|{x['name']}|{x['span']}|`{x['F']}`|`{x['R']}`|{x['covered']}/{x['target_species']}|{x['exact']}/{x['target_species']}|{x['WC']}/{x['GT']}|{x['WCloose']}/{x['GTloose']}|")
   lines+=['','|Pair|Reference|F start–end|R start–end|Deg F/R|','|---|---|---:|---:|---:|']
   for x in rs:
    if x['gene']==gene:lines.append(f"|{x['name']}|{x['reference']}|{x['F_start']}–{x['F_end']}|{x['R_start']}–{x['R_end']}|{x['deg_F']}/{x['deg_R']}|")
  lines+=['','[Methods and background profiles](../../docs/METHODS'+('.en' if en else '')+'.md)']
  if taxon=='Aphroditiformia':lines+=['','[Original base-matching figures and reference-coverage analysis](../../studies/Aphroditiformia/'+suffix+')']
  if any(x['target_taxon']==taxon for x in alt):lines+=['','[Conditional alternatives](../../conditional/'+suffix+')']
  (R/'panels'/taxon/suffix).write_text('\n'.join(lines)+'\n')
 # Conditional catalogue.
 lines=['# Conditional alternatives' if en else '# 条件备选','', '[Catalogue](../'+suffix+') · [FASTA](primers.fasta) · [CSV](results.csv) · [Species](species.csv)','',
 'These four ordinary/D sets improve or retain all four target metrics relative to their standard counterparts. They retain strict WC≤1%, GT≤3% and the close-relative gate, but treat loose-model hits as diagnostics rather than a hard veto. High loose-model counts must be considered alongside the apparent strict-model specificity. Identical alternatives were excluded; the suffix 2 is an alternative identifier.' if en else '这四组普通/D方案相对现行对应方案，覆盖与完全匹配四项均不退步。它们保留严格WC≤1%、GT≤3%及近缘非目标阈值，但不再用宽松模型命中作硬性否决。严格模型中的低命中须与宽松模型中的较高命中一起解读。与现行序列完全相同的备选已去重；数字2仅标识另一套候选。','',
 '|Pair|F (5′→3′)|R (5′→3′)|Covered|Exact|Strict WC/GT|Loose WC/GT|Background N|','|---|---|---|---:|---:|---:|---:|---:|']
 for x in alt:lines.append(f"|{x['name']}|`{x['F']}`|`{x['R']}`|{x['covered']}/{x['target_species']}|{x['exact']}/{x['target_species']}|{x['WC']}/{x['GT']}|{x['WCloose']}/{x['GTloose']}|{x['background_species']}|")
 (R/'conditional'/suffix).write_text('\n'.join(lines)+'\n')
(R/'docs').mkdir(exist_ok=True)
cn='''# 方法与复现

[目录](../README.md) · [English](METHODS.en.md)

## 设计与命名

每个基因、每个设计单元使用冻结的具名物种参考，一物种一条序列；不同设计单元存在嵌套和重复。MAFFT定向比对后搜索参考窗口，普通/D四个端点一致；每条简并度≤8，简并限于3′半段。现行入选跨度按S 100–200 bp、M 201–300 bp、L 301–480 bp展示。原Aphroditiformia和早期目级搜索的范围较窄，具体搜索脚本和已选坐标保留在源分析中，不以新的分档范围倒推旧搜索空间。

预测最低Tm≥45°C、F/R的Tm中点差≤5°C、二级结构Tm<55°C，在50 nM/变体与200 nM/方向两种情景下检查。新扩展搜索使用18–30 nt长度和1 nt起点步长，以覆盖与完全匹配进行分阶段筛选。普通/D共享坐标并不要求普通序列必须包含于D的展开序列中。

## 匹配定义

目标覆盖：双端各≤2错配、各自3′末三位兼容，产物100–500 bp且内部至少50 bp。完全匹配：双端零错配，IUPAC兼容视为匹配。分母包含不能跨越整个窗口的参考片段。背景扫描允许内部间隔0–1000 bp；严格WC/GT各≤2错配并要求末三位兼容，宽松WC允许≤3错配、宽松GT允许≤2错配且不约束末端。GT模型在引物方向额外允许G:A与T:C兼容，是敏感性分析，不是实测扩增概率。

## 背景与选择条件

- Aphroditiformia原研究：六类群多毛类背景；严格WC和GT各≤10、宽松各≤50；不能将其称为后续扩展海洋背景评估。
- Phyllodocida目级现行结果：扩展海洋背景，严格WC和GT各≤10、宽松各≤50；目录保留原选择与对应背景。
- 三个新增亚目和三个科：海洋背景加目标之外的Phyllodocida。严格WC≤1%、GT≤3%，宽松WC≤5%、GT≤10%；近缘非目标严格WC/GT各≤max(2, floor(3%×近缘参考数))。
- 最新Syllidae 16S与Nereididae 28S：现行方案要求普通/D的覆盖与完全匹配四项不低于前轮，再按D覆盖、D完全匹配、普通覆盖、普通完全匹配、较低GT背景命中排序。第二套条件备选保留严格及近缘阈值，取消宽松模型硬阈值，另列结果。

后续海洋背景以既有WoRMS海洋属筛查和去重参考整理为基础，并非逐物种、逐采样地点的完整海洋属性核验。各面板的背景组成与参考分母独立保留，不混算跨设计单元的覆盖率。覆盖率表示参考物种中的比例，不是所有已知物种的分类覆盖率。

## 数据布局

- `panels/`：按目标类群保存现行结果、FASTA、逐物种结果及参考库。
- `conditional/`：四组独立条件备选；不与现行FASTA混合。
- `studies/Aphroditiformia/`：原研究、图表、WoRMS参考覆盖、原始题目预印本与其可复现代码。
- `analysis/`：后续源分析的脚本、冻结输入、候选评分和最终筛选结果。大型逐候选命中矩阵、缓存、下载原始GenBank文本不重复收录；选定序列可直接重扫。
- `provenance.csv`：统一名称与源结果行的对应；`analysis/source_manifest.json`记录导入文件校验值。

## 运行

安装`requirements.txt`。`python3 scripts/build_catalog.py`从冻结结果构建统一目录，`python3 scripts/build_docs.py`同步中英文文档，`python3 scripts/check_release.py`检查名称、FASTA、面板、指标、坐标、链接与校验值。`python3 scripts/recheck_catalog.py`编译C++扫描器，逐记录对照独立Python实现，复核全部现行及条件候选；需要C++17编译器。绘图可运行`python3 scripts/plot_catalog.py`，需要matplotlib。

完整重新搜索还需要MAFFT和primer3-py。源分析保留其相对目录；先编译相应目录中的`raw_scan.cpp`，再从冻结输入运行设计脚本。完整候选命中矩阵应先重建，再运行依赖这些矩阵的旧汇总脚本。当前统一目录的复核不依赖被省略的中间矩阵。

本次整合保留既有序列和统计，仅选用已经独立复核的最新版方案并统一名称。未新增实验结果。原预印本内容和PDF未更改，其结论仅适用于Aphroditiformia。
'''
en='''# Methods and reproducibility

[Catalogue](../README.en.md) · [中文](METHODS.md)

## Design and naming

Each gene and target unit uses one frozen reference per named species. Units are nested and species overlap between panels. Direction-adjusted MAFFT alignments guide window searches. Ordinary/D counterparts share all four endpoints; degeneracy is ≤8 per oligo and limited to its 3′ half. Selected spans are displayed as S 100–200 bp, M 201–300 bp and L 301–480 bp. The original Aphroditiformia and early order-wide searches used narrower ranges; the source scripts preserve these original search spaces.

Two concentration scenarios (50 nM per variant and 200 nM per direction) require minimum predicted Tm≥45°C, F/R midpoint difference≤5°C and structure Tm<55°C. The latest expanded searches use 18–30 nt oligos and 1 nt start increments, with staged coverage- and exact-oriented candidate selection. Shared coordinates do not imply that an ordinary oligo is necessarily a member of its D counterpart's expansion.

## Matching models

Target coverage requires ≤2 mismatches per end and compatibility at each terminal three bases, product length 100–500 bp and internal sequence≥50 bp. Exact matches have zero mismatches at both ends; IUPAC-compatible bases match. All reference species remain in the denominator, including fragments that do not span a window. Background scans permit internal gaps of 0–1000 bp. Strict WC/GT allows ≤2 mismatches per end and compatible terminal three bases; loose WC allows ≤3 mismatches, loose GT≤2, without terminal constraints. GT additionally allows primer-oriented G:A and T:C compatibility. These are matching sensitivity models, not measured amplification probabilities.

## Background panels and selection

- Original Aphroditiformia study: a six-group polychaete background, with each strict WC/GT count≤10 and each loose count≤50. This is not an evaluation against the later expanded marine background.
- Current Phyllodocida order-wide results: the expanded marine background, with each strict WC/GT count≤10 and each loose count≤50. The original selection and its corresponding background are retained.
- Three additional suborders and three families: marine references plus Phyllodocida outside the target unit. Strict WC≤1%, GT≤3%; loose WC≤5%, GT≤10%; each strict close-relative count≤max(2, floor(3%×close-relative N)).
- Latest Syllidae 16S and Nereididae 28S: the current profile does not reduce ordinary/D coverage or exact-match counts relative to the preceding selection. Ranking then favors D coverage, D exact matches, ordinary coverage, ordinary exact matches and fewer GT background hits. Conditional alternatives retain strict and close-relative gates but remove loose-model hard vetoes; their loose counts remain explicit.

The later marine background uses existing WoRMS marine-genus screening and deduplicated references, rather than complete species- and locality-level verification. Each unit retains its own panel and denominator. Reference-species coverage is not coverage of every described species in that taxon, and percentages are not common-panel comparisons.

## Files and reproduction

`panels/` contains current results, FASTA, species outcomes and frozen references. `conditional/` holds four alternative ordinary/D sets. `studies/Aphroditiformia/` preserves the original study, base-matching figures, WoRMS coverage and preprint. `analysis/` contains compact source scripts, frozen inputs and candidate score tables; large per-candidate hit matrices, caches and raw downloaded GenBank text are regenerated rather than duplicated. `provenance.csv` maps unified names to source rows; `analysis/source_manifest.json` records imported-file hashes.

Install `requirements.txt`. Run `python3 scripts/build_catalog.py`, `python3 scripts/build_docs.py`, then `python3 scripts/check_release.py`. Full verification with `python3 scripts/recheck_catalog.py` requires a C++17 compiler and compares new C++ scans with an independent Python implementation record by record. `python3 scripts/plot_catalog.py` regenerates the overview using matplotlib.

Complete search regeneration additionally needs MAFFT and primer3-py. Keep the source directory layout, compile the corresponding `raw_scan.cpp`, and run the design scripts on the frozen inputs. Regenerate full candidate hit matrices before invoking legacy finish scripts that depend on them. Unified-catalogue verification does not need the omitted intermediate matrices.

Integration selects already independently verified results and unifies names; it adds no experimental results. The retained preprint and PDF are unchanged and report Aphroditiformia only.
'''
(R/'docs/METHODS.md').write_text(cn);(R/'docs/METHODS.en.md').write_text(en)
print('Built bilingual catalogue, 8 panel pages and conditional documentation')
