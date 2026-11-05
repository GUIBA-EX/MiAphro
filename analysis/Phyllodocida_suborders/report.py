from pathlib import Path
import csv,collections,json,subprocess,sys
from Bio import SeqIO
r=Path(__file__).resolve().parent
read=lambda p:list(csv.DictReader(p.open()))
lines=['# 三个 Phyllodocida 亚目的 16S / 28S 引物候选','', '作者：Fei Xia · 2026-11-05','', '目标按 WoRMS 亚目及科级异名归属拆分。普通/D 同坐标，每条简并度≤8，简并仅在3′半段。S 100–200 bp，M 201–300 bp，L 301–480 bp。每条参考序列对应一个命名物种。','', '背景为海洋优先参考库与全部其他 Phyllodocida 参考。严格 WC≤1%、GT≤3%；宽松WC≤5%、GT≤10%；近缘非目标WC/GT各不超过 max(2, floor(3%×其他Phyllodocida物种数))。这些是计算筛选阈值，不是实验扩增概率。','', '## 参考库','','| 亚目 | 基因 | 目标物种 | 背景物种 | 已覆盖科 |','|---|---|---:|---:|---|']
for x in read(r/'reference_summary.csv'):lines.append('| '+' | '.join(x[k] for k in ['suborder','gene','target_species','background_species','families'])+' |')
lines+=['','Glyceriformia 缺少 Lacydoniidae；Nereidiformia 缺少 Microphthalmidae；Phyllodociformia 当前仅有 Phyllodocidae。每科一物种一票，参考物种多的科对优化目标贡献较大。','', '## D候选概览','','| 亚目 | 基因/档 | 跨度bp | 覆盖 | 完全匹配 | 背景WC/GT | 近缘WC/GT |','|---|---|---:|---:|---:|---:|---:|']
allrows=[];groups=[];spans=[]
for name in ['Glyceriformia','Nereidiformia','Phyllodociformia']:
 root=r/name;sys.path.insert(0,str(root/'scripts'));from chemistry import deg,thermal
 rows=read(root/'Expanded/results.csv');by={x['name']:x for x in rows}
 for x in rows:
  gene,band=x['name'].split('-')[1:];band=band.removesuffix('D');base=f'{name}-{gene}-{band}'
  al={z.id.removeprefix('_R_'):str(z.seq) for z in SeqIO.parse(root/gene/'aligned.fa','fasta')};cols=[i for i,b in enumerate(al[x['reference']]) if b!='-'];lo=cols[int(x['F_start'])-1];hi=cols[int(x['R_end'])-1];spans.append(dict(name=x['name'],total=len(al),spanning=sum(bool(z[:lo+1].strip('-')) and bool(z[hi:].strip('-')) for z in al.values())))
  assert all(x[k]==by[base][k] for k in ['F_start','F_end','R_start','R_end'])
  assert all(deg(q)<=8 and set(q[:len(q)//2])<=set('ACGT') for q in [x['F'],x['R']]);assert thermal(x['F'],x['R'])[0]
  md={z['accession']:z for z in read(root/gene/'background.csv')};den=collections.Counter(z['group'] for z in md.values());ct=collections.Counter()
  with (root/'Expanded'/gene/band/'independent_hits.tsv').open() as h:
   for z in csv.DictReader(h,delimiter='\t'):
    if z['name']==x['id'] and z['role']=='background':ct[md[z['accession']]['group'],z['model']]+=int(z['le2last3'])
  for group,n in den.items():groups.append(dict(name=x['name'],group=group,species=n,WC=ct[group,'WC'],GT=ct[group,'GT']))
  if x['name'].endswith('D'):
   near=f"{ct['other Phyllodocida','WC']}/{ct['other Phyllodocida','GT']}"
   lines.append(f"|{name}|{gene}/{band}|{x['span']}|{x['covered']}/{x['target_species']} ({100*int(x['covered'])/int(x['target_species']):.1f}%)|{x['exact']}/{x['target_species']}|{x['WC']}/{x['GT']}|{near}|")
  allrows.append(x)
 lines+=[]
lines+=['','## 全部序列（5′→3′）','','| 名称 | F | R |','|---|---|---|']
for x in allrows:lines.append(f"|{x['name']}|`{x['F']}`|`{x['R']}`|")
lines+=['','## 各科覆盖','','| 名称 | 科 | 覆盖 | 完全匹配 |','|---|---|---:|---:|']
for name in ['Glyceriformia','Nereidiformia','Phyllodociformia']:
 for x in read(r/name/'Expanded/families.csv'):
  if x['name'].endswith('D'):lines.append(f"|{x['name']}|{x['family']}|{x['covered']}/{x['species']}|{x['exact']}|")
lines+=['','## 方法与复核','','每位点搜索18/20/22/24 nt，使用前6个高频单倍型及共识序列扩展候选；每档最多600个同坐标窗口进入扫描。目标模型：每端≤2错配且最后3碱基兼容，100–500 bp产物、内部≥50 bp；完全匹配要求双端零错配。热力学双浓度情景最低Tm≥45°C，双端Tm中点差≤5°C，结构Tm<55°C。参考窗口相对当前参考片段，不代表所有物种的实际产物长度。分母包含不完整序列。','','候选按D覆盖、普通覆盖、合计完全匹配排序。匹配结果经独立Python和C++逐记录复核。属级海洋属性不等同于具体采样地点或物种生境核验；背景未含藻类与微生物。','', '[序列FASTA](primers.fasta) · [指标及坐标](results.csv) · [分组非目标](background_groups.csv) · [参考片段跨越窗口的数量](reference_span_support.csv)','', 'WoRMS 分类快照、源比对和各物种结果分别位于 taxonomy/ 和各亚目目录。']
for fn,rs in [('results.csv',allrows),('background_groups.csv',groups),('reference_span_support.csv',spans)]:
 with (r/fn).open('w') as h:w=csv.DictWriter(h,rs[0],lineterminator='\n');w.writeheader();w.writerows(rs)
(r/'primers.fasta').write_text(''.join(f">{x['name']}_{end}\n{x[end]}\n" for x in allrows for end in ['F','R']))
(r/'README.md').write_text('\n'.join(lines)+'\n');subprocess.run(['pandoc','--standalone','--metadata','title=Suborder primer candidates',str(r/'README.md'),'-o',str(r/'report.html')],check=True)
(r/'validation.json').write_text(json.dumps(dict(status='PASS',pairs=len(allrows),independent_scan='per-record C++/Python checks in each finish_expanded.py',constraints='same endpoints; max degeneracy 8 per oligo; 3prime half only'),indent=2))
print('PASS',len(allrows),'pairs; report completed.')
