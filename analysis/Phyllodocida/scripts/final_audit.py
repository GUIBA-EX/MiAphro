from pathlib import Path
import csv,collections,subprocess,json,hashlib
root=Path(__file__).resolve().parents[1]
read=lambda p:list(csv.DictReader(p.open()))
res=read(root/'SML/results.csv');groups=[];missing_bands=[]
for gene in ['16S','28S']:
 meta=read(root/gene/'background.csv');assert len(meta)==len({r['species'] for r in meta})==len({r['accession'] for r in meta})
 groups.extend(dict(gene=gene,group=g,species=n) for g,n in collections.Counter(r['group'] for r in meta).items())
 for band in 'SML':
  d=root/'SML'/gene/band
  assert (d/'background.fasta').read_bytes()==(root/gene/'background.fa').read_bytes()
  rows=[r for r in res if r['name'].removesuffix('D')==f'Phyllodocida-{gene}-{band}'];assert len(rows) in (0,2),(gene,band)
  if not rows:missing_bands.append(gene+'-'+band)
  for r in rows:
   assert r['pass']=='True'
   details=[x for x in read(d/'marine_by_group.csv') if x['id']==r['id']]
   for model in ['WC','GT','WCloose','GTloose']:assert sum(int(x[model]) for x in details)==int(r[model])
with (root/'background_groups.csv').open('w') as h:w=csv.DictWriter(h,groups[0]);w.writeheader();w.writerows(groups)
md=root/'README.md';s=md.read_text();s+='\n## 海水与潮间带背景组成\n\n| 基因 | 门类 | 物种数 |\n|---|---|---:|\n'
for r in groups:s+=f"| {r['gene']} | {r['group']} | {r['species']} |\n"
s+='\n按 WoRMS 属级 isMarine=1 核验，排除 Phyllodocida；这表示海洋优先分类背景，不等同于全部物种均经海洋生境或具体潮间带采集记录核实。16S 海绵参考只有1种，藻类与微生物尚未纳入。原始下载记录仅用于追溯；最终筛选只使用当前 background.fa，未保留旧背景面板结果。\n'
s+='\n未找到普通与D均通过门槛的长度档：'+(', '.join(missing_bands) or '无')+'。\n'
s+='\n## 各档覆盖缺口\n\n'
for r in res:
 if not r['name'].endswith('D'):continue
 fam=[x for x in read(root/'SML/families.csv') if x['name']==r['name']]
 missing=[x['family']+' ('+x['species']+')' for x in fam if int(x['covered'])==0]
 s+='- '+r['name']+' 完全漏检的科级标签：'+(', '.join(missing) or '无')+'。\n'
s+='\n## 非目标阈值比较\n\n[探索性权衡表](SML/tradeoffs.csv) 分别列出仅WC与同时限制WC/GT的0.3%、1%、2%、5%背景阈值。该表未应用宽松错配和分组门槛，不等同于最终入选。每档按D覆盖优先、完全匹配次之选出对照。\n'
md.write_text(s)
subprocess.run(['pandoc','--standalone','--metadata','title=Phyllodocida primer candidates',str(md),'-o',str(root/'report.html')],check=True)
manifest={}
for f in [root/'README.md',root/'SML/results.csv',root/'SML/primers.fasta',root/'background_groups.csv']+[root/g/(role+'.fa') for g in ['16S','28S'] for role in ['target','background']]:manifest[str(f.relative_to(root))]=hashlib.sha256(f.read_bytes()).hexdigest()
(root/'validation.json').write_text(json.dumps(dict(status='PASS',pairs=len(res),missing_bands=missing_bands,independent_raw_scan='Python and C++ per-record agreement checked in finish_sml.py',checksums=manifest),indent=2))
print('PASS: final marine panel, all selected pairs, grouped background totals, checksums.')
