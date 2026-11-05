"""Generate bilingual release documentation from the curated catalogue."""
from pathlib import Path
import csv,json
R=Path(__file__).resolve().parents[1]
read=lambda p:list(csv.DictReader(p.open()))
rows=read(R/'catalog.csv');panels=read(R/'panels.csv');taxa=list(dict.fromkeys(x['target_taxon'] for x in panels));decisions=read(R/'archive/selection_decisions.csv')
TITLE='MiPhyllo: taxon-stratified 16S and 28S PCR primer candidates for environmental DNA metabarcoding of Phyllodocida'
notes={
'Aphroditiformia':('16S仅保留S；28S保留S/M/L。16S SD覆盖62/94种，完全匹配仅1种。','Only S is retained for 16S; S/M/L are retained for 28S. The 16S SD covers 62/94 species but exactly matches only one.'),
'Glyceriformia':('16S保留S/M，28S保留S/M/L；16S两个D均仅1种完全匹配。参考未包含Lacydoniidae。','16S retains S/M and 28S retains S/M/L; both 16S D pairs exactly match only one species. Lacydoniidae is absent from the reference panel.'),
'Phyllodociformia':('16S保留S/L，28S保留S/M。参考仅代表Phyllodocidae。','16S retains S/L and 28S retains S/M. The references represent Phyllodocidae only.'),
'Syllidae':('仅保留16S S/SD；28S和其余长度归档。','Only 16S S/SD is retained; 28S and the other lengths are archived.'),
'Nereididae':('16S保留S/M/L，28S仅保留S。16S LD覆盖最高，但背景命中高于SD。','16S retains S/M/L; 28S retains S only. The 16S LD has the highest coverage but more background matches than SD.'),
'Hesionidae':('16S与28S均保留S/M/L。28S SD覆盖40/45种，完全匹配8种。','Both genes retain S/M/L. The 28S SD covers 40/45 species with eight exact matches.')}
for en in [False,True]:
 f='README.en.md' if en else 'README.md';other='README.md' if en else 'README.en.md';lang=int(en)
 rule=('Retain a standard ordinary/D set when D covers at least 50% of its reference species and exactly matches at least one, under its original full background selection profile. Phyllodocida and Nereidiformia broad-target designs and all relaxed-background alternatives are archived. This is an editorial catalogue rule, not an experimentally established efficacy threshold. Ordinary primers are retained as same-coordinate comparators; the D criterion does not imply that ordinary coverage is also ≥50%.' if en else '保留D覆盖率≥50%、至少有1种完全匹配且通过原完整背景筛选的普通/D组合。Phyllodocida、Nereidiformia大范围方案及放宽背景限制的条件备选归档。这是目录精简规则，不是实验有效性阈值。普通引物保留为同坐标对照，D达到标准不代表普通引物也达到50%覆盖率。')
 lines=['![MiPhyllo](figures/MiPhyllo_overview.png)','', '# '+TITLE,'','**Fei Xia**','',f'[English]({other}) | 中文' if not en else f'[中文]({other}) | English','', '[Preprint PDF](publications/MiPhyllo_preprint.pdf) · [Manuscript](publications/MiPhyllo_preprint.md)','', '## Overview' if en else '## 项目概览','',
 f'The curated catalogue contains **{len(rows)} primer pairs ({len(rows)*2} oligos), arranged in {len(rows)//2} same-coordinate ordinary/D sets across six target units**: three suborders and three families. The original search considered eight units. Each retained unit keeps its own species-deduplicated reference and background panels. Phyllodocida remains the project scope; no order-wide primer is retained.' if en else f'精简目录保留**{len(rows)}对引物（{len(rows)*2}条寡核苷酸），组成{len(rows)//2}组同坐标普通/D方案，覆盖6个设计单元**：3个亚目和3个科。原始搜索涉及8个单元。各单元独立保存一物种一序列的参考库与背景库。Phyllodocida表示项目研究范围，现行目录不再保留目级引物。','',
 '[Current FASTA](primers.fasta) · [Catalogue](catalog.csv) · [Species outcomes](species.csv) · [Methods](docs/METHODS.en.md) · [Archive and selection](archive/README.en.md)' if en else '[现行序列](primers.fasta) · [现行目录](catalog.csv) · [逐物种结果](species.csv) · [方法](docs/METHODS.md) · [归档与取舍](archive/README.md)','',
 'The banner illustrates taxonomic context. Nereidiformia is represented by three family panels, not a retained suborder-wide assay.' if en else '顶部简介图展示分类背景；Nereidiformia由三个科级面板代表，不表示保留了亚目通用引物。','',
 '## Retained designs' if en else '## 保留方案','',
 '|Target|16S lengths|28S lengths|16S references|28S references|Results|' if en else '|目标类群|16S长度|28S长度|16S参考物种|28S参考物种|结果入口|','|---|---|---|---:|---:|---|']
 for t in taxa:
  ps={x['gene']:x for x in panels if x['target_taxon']==t};bands=lambda g:'/'.join(x['band'] for x in rows if x['target_taxon']==t and x['gene']==g and x['variant']=='D') or '—'
  lines.append(f"|{t}|{bands('16S')}|{bands('28S')}|{ps['16S']['target_species']}|{ps['28S']['target_species']}|[{t}](panels/{t}/{f})|")
 lines+=['', 'Reference counts describe the available libraries; a count does not imply a retained primer for that gene. Each listed length includes its ordinary and D pair.' if en else '参考物种数表示可用参考库规模，不表示该基因一定保留引物。每个保留长度均含普通与D两对。','', '## Selection and naming' if en else '## 筛选与命名','',rule,'',
 '`MiPhyllo-{taxon}-{gene}-{S|M|L}[D]_{F|R}`. Ordinary primers contain only A/C/G/T; D marks degeneracy in at least one oligo. Each oligo has degeneracy ≤8, restricted to its 3′ half. Ordinary/D counterparts share all four binding endpoints. S/M/L denote length classes, not mandatory sets.' if en else '`MiPhyllo-{类群}-{基因}-{S|M|L}[D]_{F|R}`。普通引物仅含A/C/G/T；D表示至少一条含简并。每条简并度≤8，简并限于3′半段。普通/D四个结合端点相同；S/M/L是长度分档，不要求凑齐。','',
 '## Results' if en else '## 结果概览','', '![D coverage and exact matching](figures/catalog_coverage.svg)','',
 'Filled markers show D coverage; open markers show exact matching. Grey connectors show the gap between them. Missing lengths were removed or had no original candidate. Percentages use each unit’s own full reference denominator and are not common-panel comparisons.' if en else '实心点为D覆盖率，空心点为完全匹配率，灰线表示二者差距。缺失长度表示已移除或原本无候选。各自完整参考库作为分母，不作共同面板性能比较。','']
 for t in taxa:lines.append('- **'+t+'**: '+notes[t][lang])
 lines+=['','## Reproducibility' if en else '## 复现与证据','',
 '[Selection decisions](archive/selection_decisions.csv) · [Curation audit](audit/curation.json) · [Independent scan](audit/catalog_rescan.json) · [Provenance](provenance.csv) · [Checksums](SHA256SUMS)','',
 'The 48 retained pairs are sequence-identical subsets of the 102 independently re-scanned original pairs. This revision changes selection and presentation, not sequences or matching outcomes. Original results and the previous integrated PDF remain in the archive. No experimental validation is added.' if en else '48对保留引物均为原102对独立复核引物的原序列子集。本次变更筛选与展示，不修改序列和匹配结果。原始结果与上一版整合PDF保留在归档中，未新增实验验证。','',
 '```bash','python3 -m pip install -r requirements.txt','python3 scripts/build_catalog.py','python3 scripts/build_docs.py','python3 scripts/plot_catalog.py  # requires matplotlib','python3 publications/render_integrated.py  # requires pandoc and Chrome','python3 scripts/update_checksums.py','python3 scripts/check_release.py','```','',
 '## Historical study' if en else '## 原始研究','',f'[Aphroditiformia study](studies/Aphroditiformia/{f}) · [Historical PDF](studies/Aphroditiformia/publications/MiAphro_preprint.pdf)','',
 'The historical study preserves the original complete search and figures. Use the current catalogue above for the curated selection.' if en else '原研究保留完整搜索结果与图表；精简后的保留方案以本页现行目录为准。','',
 '## Citation and copyright' if en else '## 引用与著作权','',f'Xia, F. (2026). *{TITLE}*. [MiPhyllo](https://github.com/GUIBA-EX/MiPhyllo).','', '[Citation metadata](CITATION.cff) · [Copyright and permissions](COPYRIGHT.md)','']
 (R/f).write_text('\n'.join(lines))
 for p in (R/'panels').iterdir():
  if not p.is_dir():continue
  t=p.name;rs=[x for x in rows if x['target_taxon']==t]
  ls=[f'# MiPhyllo · {t}','',f'[Catalogue](../../{f}) · [Archive](../../archive/{f})','']
  if not rs:
   ls += ['This broad-target design is archived and has no current primer pairs. Reference inputs remain available here; complete original results are in the archive.' if en else '该大范围设计已归档，不保留现行引物。参考输入仍保留于本目录；完整原始结果见归档。','']
  else:
   ls += [notes[t][lang],'','[FASTA](primers.fasta) · [CSV](results.csv) · [Species](species.csv)','']
  for g in ['16S','28S']:
   gs=[x for x in rs if x['gene']==g];ls+=['## '+g,'',f'[Target FASTA](references/{g}/target.fa) · [Metadata](references/{g}/target.csv) · [Alignment](references/{g}/aligned.fa) · [Background](references/{g}/background.fa)','']
   if not gs:ls+=['No retained candidate.' if en else '无保留候选。',''];continue
   ls+=['|Pair|bp|F (5′→3′)|R (5′→3′)|Covered|Exact|Strict WC/GT|Loose WC/GT|','|---|---:|---|---|---:|---:|---:|---:|']
   for x in gs:ls.append(f"|{x['name']}|{x['span']}|`{x['F']}`|`{x['R']}`|{x['covered']}/{x['target_species']}|{x['exact']}/{x['target_species']}|{x['WC']}/{x['GT']}|{x['WCloose']}/{x['GTloose']}|")
   ls+=['','|Pair|Reference|F start–end|R start–end|Deg F/R|','|---|---|---:|---:|---:|']
   for x in gs:ls.append(f"|{x['name']}|{x['reference']}|{x['F_start']}–{x['F_end']}|{x['R_start']}–{x['R_end']}|{x['deg_F']}/{x['deg_R']}|")
   ls+=['',f"Background reference species: {gs[0]['background_species']}",'']
  ls += [f'[Methods](../../docs/METHODS'+('.en' if en else '')+'.md)','']
  (p/f).write_text('\n'.join(ls))
 # Explicit archive boundary; original outputs are never overwritten by current exports.
 ls=['# Archive and selection' if en else '# 归档与取舍','',f'[Current catalogue](../{f})','',rule,'',
 'All 94 original standard pairs and eight conditional pairs are preserved. The current catalogue retains 48 pairs; 54 pairs are archived. Removal changes recommendation status, not the recorded results.' if en else '完整保留原94对标准候选与8对条件备选。现行目录保留48对，另54对归档。移除仅改变推荐状态，不删除原始结果。','',
 '[Original catalogue](2026-11-05-unfiltered/catalog.csv) · [Original FASTA](2026-11-05-unfiltered/primers.fasta) · [Species outcomes](2026-11-05-unfiltered/species.csv) · [Conditional sequences](2026-11-05-unfiltered/conditional/primers.fasta) · [Previous PDF](2026-11-05-unfiltered/publications/MiPhyllo_preprint.pdf)','',
 '[Removed candidates](removed_candidates.csv) · [Removed FASTA](removed_primers.fasta) · [All decisions](selection_decisions.csv)','',
 '|Set|D covered|D exact|Decision|Reason|','|---|---:|---:|---|---|']
 for d in decisions:ls.append(f"|{d['set']}|{d['covered']}/{d['target_species']}|{d['exact']}/{d['target_species']}|{d['status']}|{d['reason']}|")
 (R/'archive'/f).write_text('\n'.join(ls)+'\n')
 (R/'conditional'/f).write_text(('# Archived conditional alternatives' if en else '# 已归档条件备选')+f'\n\n[Archive](../archive/{f})\n\n'+('No conditional alternatives remain in the current catalogue.' if en else '现行目录不再保留条件备选。')+'\n')
(R/'publications/README.md').write_text('# MiPhyllo curated preprint / 精简版预印本\n\n[PDF](MiPhyllo_preprint.pdf) · [Manuscript](MiPhyllo_preprint.md) · [HTML](MiPhyllo_preprint.html)\n\nFei Xia · 6 target units · 24 ordinary/D sets · 48 primer pairs · 96 oligos.\n\n[Previous integrated PDF / 原整合版PDF](../archive/2026-11-05-unfiltered/publications/MiPhyllo_preprint.pdf) · [Historical Aphroditiformia PDF / 原鳞虫PDF](MiAphro_preprint.pdf)\n')
print('Built bilingual curated documents and archive decisions')
