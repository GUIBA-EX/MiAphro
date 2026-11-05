from pathlib import Path
import csv,collections,json,subprocess
from chemistry import deg,thermal
r=Path(__file__).resolve().parents[1];d=r/'Expanded'
read=lambda p:list(csv.DictReader(p.open()))
rows=read(d/'results.csv');old={x['name']:x for x in read(r/'SML/results.csv')}
lines=['# Phyllodocida 扩大搜索','', '沿用海洋优先背景库，不使用旧背景面板。普通/D 同坐标，每条简并度≤8，简并限3′半段。S 100–200 bp，M 201–300 bp，L 301–480 bp。引物18/20/22/24 nt；最多1400窗口/档。','', '严格WC背景≤1%，严格GT≤3%；宽松WC≤5%，宽松GT≤10%。GT为单独容许配对模型，不能当作实验非目标扩增率。按D覆盖数、普通覆盖数、双方案完全匹配数依次排序；一物种一票，样本物种多的科权重较高。','', '| 候选 | 跨度 | 原覆盖 | 新覆盖 | 新完全匹配 | WC背景 | GT背景 |','|---|---:|---:|---:|---:|---:|---:|']
groupdetail=[]
for x in rows:
 a=old.get(x['name']);n=int(x['target_species']);bg=int(x['background_species']);lines.append(f"|{x['name']}|{x['span']}|{a['covered'] if a else '—'}/{n}|{x['covered']}/{n} ({100*int(x['covered'])/n:.1f}%)|{x['exact']}/{n}|{x['WC']}/{bg}|{x['GT']}/{bg}|")
 for q in [x['F'],x['R']]:assert deg(q)<=8 and set(q[:len(q)//2])<=set('ACGT')
 assert thermal(x['F'],x['R'])[0]
 gene,band=x['name'].split('-')[1:];band=band.removesuffix('D');meta={z['accession']:z for z in read(r/gene/'background.csv')};tot=collections.Counter(z['group'] for z in meta.values());counts=collections.Counter()
 with (d/gene/band/'independent_hits.tsv').open() as h:
  for z in csv.DictReader(h,delimiter='\t'):
   if z['name']==x['id'] and z['role']=='background':counts[meta[z['accession']]['group'],z['model']]+=int(z['le2last3'])
 for g,n in tot.items():groupdetail.append(dict(name=x['name'],group=g,total=n,WC=counts[g,'WC'],GT=counts[g,'GT']))
controls=[]
for gene in ['16S','28S']:
 for band in 'SML':
  pool=read(r/'SML'/gene/band/'screen.csv');by={x['id']:x for x in pool};n=int(pool[0]['background_species'])
  def passes(x):return int(x['WC'])<=int(n*.01) and int(x['GT'])<=int(n*.03) and int(x['WCloose'])<=int(n*.05) and int(x['GTloose'])<=int(n*.10)
  pairs=[(by[x['id'][:-1]],x) for x in pool if x['id'].endswith('D') and passes(x) and passes(by[x['id'][:-1]])]
  if pairs:
   a,b=max(pairs,key=lambda v:(int(v[1]['covered']),int(v[0]['covered']),sum(int(x['exact']) for x in v)))
   controls.append(dict(gene=gene,band=band,covered=b['covered'],exact=b['exact'],WC=b['WC'],GT=b['GT'],F=b['F'],R=b['R']))
with (d/'selection_controls.csv').open('w') as h:w=csv.DictWriter(h,controls[0]);w.writeheader();w.writerows(controls)
lines+=['','## 同门槛对照','','| 基因/档 | 旧搜索重新筛选的D覆盖 | 新搜索D覆盖 |','|---|---:|---:|']
for c in controls:
 new=next((x for x in rows if x['name']==f"Phyllodocida-{c['gene']}-{c['band']}D"),None)
 lines.append(f"|{c['gene']}/{c['band']}|{c['covered']}|{new['covered'] if new else '—'}|")
lines+=['','此表两侧使用相同的放宽阈值和排序规则；旧搜索与新搜索长度范围不同。']
lines+=['','## 序列（5′→3′）','','| 候选 | F | R |','|---|---|---|']
for x in rows:lines.append(f"|{x['name']}|`{x['F']}`|`{x['R']}`|")
lines+=['','## 各科覆盖','','| 候选 | 科 | 覆盖/参考种数 | 完全匹配 |','|---|---|---:|---:|']
for x in read(d/'families.csv'):
 if x['name'].endswith('D'):lines.append(f"|{x['name']}|{x['family']}|{x['covered']}/{x['species']}|{x['exact']}|")
lines+=['','[逐物种结果](species.csv) · [分组非目标命中](background_groups.csv) · [完整坐标与指标](results.csv) · [FASTA](primers.fasta)','', '分母含不完整参考序列；科名沿用源记录。结果是计算预测，独立Python与C++逐记录一致性已验证。搜索为有限候选搜索，不代表所有可能引物的全局最优。']
with (d/'background_groups.csv').open('w') as h:w=csv.DictWriter(h,groupdetail[0]);w.writeheader();w.writerows(groupdetail)
(d/'README.md').write_text('\n'.join(lines)+'\n');subprocess.run(['pandoc','--standalone','--metadata','title=Expanded Phyllodocida search',str(d/'README.md'),'-o',str(d/'report.html')],check=True)
print('PASS: expanded results, background groups and report written.')
