from pathlib import Path
import csv,collections,subprocess,sys,json
from Bio import SeqIO
r=Path(__file__).resolve().parent;sys.path.insert(0,str(r/'scripts'));from chemistry_round2 import deg,thermal,terminal
read=lambda p:list(csv.DictReader(p.open()))
rows=read(r/'Round2/results.csv');old={x['name']:x for x in read(r/'Expanded/results.csv')};lines=['# Nereidiformia · 第二轮扩大搜索','', '沿用现有目标库和海洋背景，普通/D同坐标，简并仅限3′半段且每条≤8。S 100–200 bp，M 201–300 bp，L 301–480 bp。','', '新增26/28 nt长度、最多12个高频单倍型和共识种子；对简并分配采用保留多个分支的搜索，按整条引物的覆盖及完全匹配评分。同窗口保留多个普通与D组合供热力学筛选，每档最多1400个窗口。末端GC条件由最后3位放宽到最后5位，其余热力学要求保持不变。','', '最终在上一轮和本轮候选的并集中选配。严格WC背景≤1%、GT≤3%；宽松WC≤5%、GT≤10%；近缘非目标WC/GT各≤max(2, floor(3%×近缘参考数))。按D覆盖、普通覆盖、合计完全匹配排序。','', '本轮扩大搜索后，16S D最高覆盖62/315（19.7%），28S D最高覆盖45/133（33.8%），提升仍有限。主要缺口是Syllidae：16S各长度最多覆盖19/180，28S最多1/28；尚不作为Nereidiformia亚目通用方案推荐。按覆盖优先选出的候选存在取舍：16S SD完全匹配由31降至15，28S SD由8降至0；16S L普通引物覆盖由58降至48。上一轮结果及其他排序备选保留供比较。','', '## 新旧对比','','|候选|参考跨度|旧覆盖|新覆盖|旧完全匹配|新完全匹配|旧背景WC/GT|新背景WC/GT|','|---|---:|---:|---:|---:|---:|---:|---:|'];grouped=[]
for x in rows:
 a=old[x['name']];n=int(x['target_species']);lines.append(f"|{x['name']}|{x['span']}|{a['covered']}/{n}|{x['covered']}/{n} ({100*int(x['covered'])/n:.1f}%)|{a['exact']}|{x['exact']}|{a['WC']}/{a['GT']}|{x['WC']}/{x['GT']}|")
 for q in [x['F'],x['R']]:assert deg(q)<=8 and terminal(q) and set(q[:len(q)//2])<=set('ACGT')
 assert thermal(x['F'],x['R'])[0]
 gene,band=x['name'].split('-')[1:];band=band.removesuffix('D');bg={z['accession']:z for z in read(r/gene/'background.csv')};ct=collections.Counter()
 # Determine which search supplied this selected pair by sequence and ID, then audit grouped background hits.
 for folder in ['Round2','Expanded']:
  pool=read(r/folder/gene/band/'screen.csv')
  if any(z['id']==x['id'] and z['F']==x['F'] and z['R']==x['R'] for z in pool):break
 with (r/folder/gene/band/'independent_hits.tsv').open() as h:
  for z in csv.DictReader(h,delimiter='\t'):
   if z['name']==x['id'] and z['role']=='background':ct[bg[z['accession']]['group'],z['model']]+=int(z['le2last3'])
 for g,n in collections.Counter(z['group'] for z in bg.values()).items():grouped.append(dict(name=x['name'],group=g,total=n,WC=ct[g,'WC'],GT=ct[g,'GT']))
 for m in ['WC','GT']:assert sum(v[m] for v in grouped if v['name']==x['name'])==int(x[m])
lines+=['','## 序列（5′→3′）','','|候选|F|R|','|---|---|---|']
for x in rows:lines.append(f"|{x['name']}|`{x['F']}`|`{x['R']}`|")
lines+=['','## 每科覆盖','','|候选|科|覆盖/参考物种|完全匹配|','|---|---|---:|---:|']
for x in read(r/'Round2/families.csv'):
 if x['name'].endswith('D'):lines.append(f"|{x['name']}|{x['family']}|{x['covered']}/{x['species']}|{x['exact']}|")
lines+=['','[FASTA](primers.fasta) · [坐标与指标](results.csv) · [逐物种结果](species.csv) · [背景分组](background_groups.csv) · [覆盖/完全匹配排序备选](alternative_objectives.csv)','', '覆盖分母包括全部参考物种，不因窗口缺失而缩小。每端≤2错配且末端3位兼容；完全匹配为双端零错配。不同基因参考库不同，不直接比较百分比作为标记优劣结论。匹配结果经独立Python/C++逐记录复核。']
with (r/'Round2/background_groups.csv').open('w') as h:w=csv.DictWriter(h,grouped[0]);w.writeheader();w.writerows(grouped)
(r/'Round2/README.md').write_text('\n'.join(lines)+'\n');subprocess.run(['pandoc','--standalone','--metadata','title=Nereidiformia round 2',str(r/'Round2/README.md'),'-o',str(r/'Round2/report.html')],check=True)
(r/'Round2/validation.json').write_text(json.dumps(dict(status='PASS',pairs=len(rows),models='independent C++/Python agreement',degeneracy='<=8 per oligo'),indent=2))
print('PASS',len(rows),'pairs')
