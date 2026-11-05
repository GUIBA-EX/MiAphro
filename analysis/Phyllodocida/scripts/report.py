from pathlib import Path
import csv,json,gzip,io,collections,hashlib
from Bio import SeqIO
from chemistry import deg,thermal,terminal
root=Path(__file__).resolve().parents[1]
def read(p):return list(csv.DictReader(p.open()))
summary=[]
for gene in ['16S','28S']:
 for role in ['target','background']:
  d=root/gene;rs=list(SeqIO.parse(d/(role+'.fa'),'fasta'));meta=read(d/(role+'.csv'));assert len(rs)==len(meta)==len({r['species'] for r in meta})
  downloaded=sum(1 for f in d.glob(role+'_*.gb.gz') for r in SeqIO.parse(io.StringIO(gzip.decompress(f.read_bytes()).decode()),'genbank'))
  expected=int(json.loads((d/(role+'.json')).read_text())['esearchresult']['count']);assert downloaded==expected,(gene,role,downloaded,expected)
  summary.append(dict(gene=gene,role=role,downloaded=downloaded if role=='target' else 'multiple sources',species=len(rs),genera=len({r['species'].split()[0] for r in meta}),families=len({r['family'] for r in meta})))
 for role in ['target','background']:
  pass
 assert not ({r['species'] for r in read(d/'target.csv')}&{r['species'] for r in read(d/'background.csv')})
results=read(root/'SML/results.csv');audit=[]
for r in results:
 f,rev=r['F'],r['R'];ok,tm=thermal(f,rev);assert ok and terminal(f) and terminal(rev)
 for q in [f,rev]:assert deg(q)<=8 and set(q[:len(q)//2])<=set('ACGT')
 assert min(min(t['F_tm_min'],t['R_tm_min']) for t in tm)>=45
 name=r['name'];base=name[:-1] if name.endswith('D') else name
 for other in results:
  if other['name'] in [base,base+'D']:assert all(r[k]==other[k] for k in ['F_start','F_end','R_start','R_end'])
 audit.append(dict(name=name,F_degeneracy=deg(f),R_degeneracy=deg(rev),thermal_pass=ok,**{k:r[k] for k in ['F_start','F_end','R_start','R_end']}))
for fn,rows in [('reference_summary.csv',summary),('primer_audit.csv',audit)]:
 with (root/fn).open('w') as h:w=csv.DictWriter(h,rows[0]);w.writeheader();w.writerows(rows)
text=['# Phyllodocida 16S / 28S 引物设计','', '2026-11-05 · 目级候选筛选','', '**本轮未筛选出满足广泛目级覆盖与少量非目标命中的通用方案。** 下列12对为当前搜索范围内通过非目标门槛的对照候选，D覆盖仅7.8%–17.5%，集中于少数科；不作为目级通用推荐。', '', '## 参考库','', '| 基因 | 面板 | 下载记录 | 保留物种 | 属 | 科 |','|---|---|---:|---:|---:|---:|']
for r in summary:text.append('| '+' | '.join(str(r[k]) for k in ['gene','role','downloaded','species','genera','families'])+' |')
text+=['','目标为 NCBI Phyllodocida 下带 16S/28S 标题或基因注释的 150–6000 bp 记录；提取注释 rRNA，每个已命名物种保留有效 ACGT 最长的序列，歧义比例≤3%。背景按海水/潮间带用途整理，经 WoRMS 属级海洋属性核验，排除目标目；每物种保留一条参考。查询式及原始记录保存在各基因目录。','', '## 低非目标命中的对照候选','', '| 引物 | 参考跨度 bp | 覆盖种数 | 完全匹配 | 背景 WC/GT | F 简并度 | R 简并度 |','|---|---:|---:|---:|---:|---:|---:|']
for r in results:text.append(f"| {r['name']} | {r['span']} | {r['covered']}/{r['target_species']} | {r['exact']}/{r['target_species']} | {r['WC']}/{r['GT']} | {deg(r['F'])} | {deg(r['R'])} |")
text+=['','## 序列（5′→3′）','','| 引物 | F | R |','|---|---|---|']
for r in results:text.append(f"| {r['name']} | `{r['F']}` | `{r['R']}` |")
text+=['','## 方法与选择规则','','MAFFT 定向比对。S 135–200 bp，M 201–280 bp，L 281–400 bp；长度为未加 gap 的参考片段坐标跨度。普通与 D 四个起止坐标一致，简并限于 3′ 半段，每条≤8。每个位点从前三个高频单倍型构建普通和简并候选；分散窗口搜索后，每档最多复核500个窗口。热力学使用50 nM/具体变体和200 nM/方向两个情景，最低 Tm≥45°C，双端中点 Tm 差≤5°C，预测结构 Tm<55°C。','','覆盖模型：每端≤2错配，末端3碱基兼容，100–500 bp产物且内部≥50 bp；完全匹配为双端零错配。分母包含所有参考物种。背景严格 WC、GT 各≤10种；宽松 WC（每端≤3错配）及 GT（每端≤2错配）各≤50种，不设末端约束，内部间距0–1000 bp；严格命中另按背景门类控制，每组不超过 max(2, floor(1%×组内物种数))。通过门槛后按普通+D覆盖及完全匹配总数排序；这是本次搜索范围内的优选。C++全窗口扫描与独立Python逐记录复核。','','## 结果文件','','- [每科覆盖](SML/families.csv)','- [每物种覆盖](SML/species.csv)','- [全部候选指标与坐标](SML/results.csv)','- [引物 FASTA](SML/primers.fasta)','- [审计](primer_audit.csv)','- [不同排序目标的备选](SML/alternative_objectives.csv)','','这是计算候选结果。尚未计算窗口物种分辨率或独立留出表现，背景的属级海洋属性不代表每条序列的采集地点或具体潮位带。']
(root/'README.md').write_text('\n'.join(text)+'\n')
import html
(root/'report.html').write_text('<meta charset="utf-8"><title>Phyllodocida primers</title><style>body{max-width:1100px;margin:40px auto;font:16px/1.6 sans-serif}pre{white-space:pre-wrap}</style><pre>'+html.escape('\n'.join(text))+'</pre>')
print('PASS: reference acquisition, uniqueness, coordinates, degeneracy, thermal checks; report written.')
