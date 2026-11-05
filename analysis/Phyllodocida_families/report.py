from pathlib import Path
import csv,json,sys,subprocess,hashlib
from collections import Counter
from Bio import SeqIO
R=Path(__file__).resolve().parent
sys.path.insert(0,str(R/'Syllidae/scripts'))
from chemistry_round2 import deg,terminal,thermal
read=lambda p:list(csv.DictReader(p.open()))
def write(p,rows):
 if not rows:p.write_text('name\n');return
 with p.open('w') as h:
  w=csv.DictWriter(h,rows[0],lineterminator='\n');w.writeheader();w.writerows(rows)
families=['Syllidae','Nereididae','Hesionidae'];allrows=[];seqs=[];groups=[];genera=[];validation=[];coordinate_refs=[]
base=read(R.parent/'Phyllodocida_suborders/Nereidiformia/Round2/families.csv')
lines=['# 科级16S与28S引物候选','','**后续优化：** [Syllidae 16S与Nereididae 28S扩大搜索报告](Optimization3/report.html)。以下保留首轮科级基线；新报告分别列出均衡方案与条件备选。','', 'Syllidae · Nereididae · Hesionidae','', '作者：Fei Xia','', '以每科为独立目标，沿用海洋背景库并将目标科以外的Phyllodocida参考加入非目标面板。每个命名物种保留一条现有参考；分别重新比对后搜索S（100–200 bp）、M（201–300 bp）、L（301–480 bp）参考跨度。普通与D的四个结合端点完全相同；每条D寡核苷酸简并度≤8，简并仅位于3′半段。','', '## 设计与核验','', '搜索18、20、22、24、26、28 nt引物，起点步长2 nt，以最多12个高频单倍型及共识序列为种子，进行多分支简并搜索。每档最多保留1400个共享坐标窗口，完成目标与海洋背景双链扫描；另将亚目搜索前两轮的同长度已选方案在新的科级目标与背景库上重测，并一同参与选择。最终按D覆盖、普通覆盖、合计完全匹配、较低GT背景命中依次选择。所列候选是本轮搜索中的优选，有限搜索不等同于穷举全部可能序列。','', '目标覆盖要求双端各≤2错配且各末端3位兼容，产物100–500 bp、内部片段≥50 bp；完全匹配要求双端零错配。分母为该科全部参考物种，包含不能跨越该窗口的片段。背景严格模型为每端≤2错配且末端3位兼容；宽松WC为≤3错配、宽松GT为≤2错配，不约束末端。背景扫描允许内部间隔0–1000 bp。GT允许引物方向G:A及T:C兼容，是敏感性模型，不是实测扩增概率。','', '筛选上限保持严格WC≤1%、严格GT≤3%、宽松WC≤5%、宽松GT≤10%；其他Phyllodocida严格WC与GT各≤max(2, floor(3%×近缘参考数))。热力学使用50 nM/变体与200 nM/方向两种浓度情景，最低Tm≥45°C、平均Tm差≤5°C、二级结构Tm<55°C；末端5位至少一个确定G/C，排除可能的末端3连聚与全序列5连聚。每个最终候选均用独立Python算法逐记录复核C++扫描结果。','', '## 参考库','','|科|基因|目标物种|背景物种|','|---|---|---:|---:|']
for x in read(R/'reference_summary.csv'):lines.append(f"|{x['family']}|{x['gene']}|{x['target']}|{x['background']}|")
for fam in families:
 d=R/fam;rows=read(d/'Round2/results.csv');allrows+=rows;lookup={x['name']:x for x in rows}
 lines+=['',f'## {fam}','','|基因/档位|普通覆盖|D覆盖|普通完全匹配|D完全匹配|D严格背景WC/GT|亚目方案在该科覆盖|','|---|---:|---:|---:|---:|---:|---:|']
 for gene in ['16S','28S']:
  bg={x['accession']:x for x in read(d/gene/'background.csv')};target=read(d/gene/'target.csv')
  aligned={x.id.removeprefix('_R_'):x for x in SeqIO.parse(d/gene/'aligned.fa','fasta')}
  control_aligned={x.id.removeprefix('_R_'):x for x in SeqIO.parse(d/'Controls'/gene/'coordinate_reference_alignment.fa','fasta')}
  for ref,pool in {(x['reference'],x['source_pool']) for x in rows if x['name'].split('-')[1]==gene}:
   record=(control_aligned if pool=='Controls' else aligned)[ref];coordinate_refs.append(dict(family=fam,gene=gene,reference=ref,source_pool=pool,orientation='reverse_complement' if record.id.startswith('_R_') else 'original',ungapped_length=len(str(record.seq).replace('-',''))))
  for band in 'SML':
   name=f'{fam}-{gene}-{band}'
   if name not in lookup:lines.append(f'|{gene} {band}|未找到通过候选|—|—|—|—|—|');continue
   a,b=lookup[name],lookup[name+'D'];n=int(a['target_species']);old=next(x for x in base if x['name']==f'Nereidiformia-{gene}-{band}D' and x['family']==fam)
   lines.append(f"|{gene} {band}|{a['covered']}/{n} ({100*int(a['covered'])/n:.1f}%)|{b['covered']}/{n} ({100*int(b['covered'])/n:.1f}%)|{a['exact']}/{n}|{b['exact']}/{n}|{b['WC']}/{b['GT']}|{old['covered']}/{n}|")
   coords=['reference','F_start','F_end','R_start','R_end'];assert all(a[k]==b[k] for k in coords)
   selected={a['id']:a,b['id']:b};counts=Counter();bgcounts=Counter()
   with (d/a['source_pool']/gene/band/'independent_hits.tsv').open() as h:
    for z in csv.DictReader(h,delimiter='\t'):
     if z['name'] in selected and z['role']=='background':bgcounts[z['name'],bg[z['accession']]['group'],z['model']]+=int(z['le2last3'])
   for x in [a,b]:
    for end in 'FR':
     q=x[end];assert len(q)==int(x[end+'_end'])-int(x[end+'_start'])+1
     assert deg(q)<=8 and terminal(q) and set(q[:len(q)//2])<=set('ACGT')
     if not x['name'].endswith('D'):assert deg(q)==1
     seqs.append(f">{x['name']}_{end}\n{q}\n")
    ok,reports=thermal(x['F'],x['R']);assert ok and min(min(z['F_tm_min'],z['R_tm_min']) for z in reports)>=45
    for g,total in Counter(z['group'] for z in bg.values()).items():groups.append(dict(name=x['name'],group=g,total=total,WC=bgcounts[x['id'],g,'WC'],GT=bgcounts[x['id'],g,'GT']))
    for model in ['WC','GT']:assert sum(z[model] for z in groups if z['name']==x['name'])==int(x[model])
    validation.append(dict(name=x['name'],same_coordinates=True,deg_F=deg(x['F']),deg_R=deg(x['R']),thermal_pass=True,independent_verification='PASS'))
 lines+=['','### 引物序列（5′→3′）','','|名称|F|R|简并度F/R|去gap参考坐标F；R|','|---|---|---|---:|---|']
 for x in rows:lines.append(f"|{x['name']}|`{x['F']}`|`{x['R']}`|{deg(x['F'])}/{deg(x['R'])}|{x['reference']} {x['F_start']}–{x['F_end']}；{x['R_start']}–{x['R_end']}|")
 detail=read(d/'Round2/species.csv');ct=Counter()
 for x in detail:
  genus=x['species'].split()[0]
  ct[x['name'],genus,'total']+=1;ct[x['name'],genus,'covered']+=int(x['covered']);ct[x['name'],genus,'exact']+=int(x['exact'])
 for (name,genus,k),n in ct.items():
  if k=='total':genera.append(dict(name=name,genus=genus,total=n,covered=ct[name,genus,'covered'],exact=ct[name,genus,'exact']))
 lines+=['',f'[逐物种结果]({fam}/Round2/species.csv) · [备选排序]({fam}/Round2/alternative_objectives.csv) · [搜索统计]({fam}/Round2/search.csv)']
lines+=['','## 阈值敏感性', '', '[背景阈值诊断](screening_diagnostics.csv)另列严格WC≤2%、GT≤5%的探索结果；该诊断未同时约束宽松模型与近缘非目标，不属于正式通过候选。无限制方案仅用于展示覆盖与特异性的冲突。', '', '## 数据与复现','','[全部序列FASTA](primers.fasta) · [坐标与指标](results.csv) · [逐属覆盖](genera.csv) · [背景类群命中](background_groups.csv) · [检查结果](validation.json)','', '坐标为比对方向参考序列的1-based去gap坐标；参考是否反向互补见[参考方向](reference_coordinates.csv)。新搜索与保留对照可能使用不同指定参考，坐标不可跨参考直接比较。', '', '参考序列、比对、候选扫描及独立复核结果保存在各科目录。`prepare.py`从既有海洋面板拆分；`run.py`调用`controls.py`重测原亚目方案，执行每科独立比对、六档搜索与复核；`report.py`汇总。依赖Python、Biopython、primer3-py、MAFFT、Pandoc及随附C++扫描器。参考库使用既有科级归属，当前结果仅为计算评估，未添加实验验证。']
summary=['## 结果概览','','下表列出各科、各基因按D覆盖最高选出的长度；完整S/M/L结果见后文。最高覆盖不代表最佳实验方案，应同时考虑完全匹配和非目标命中。','','|科|基因|D覆盖最高候选|覆盖/参考物种|完全匹配|严格背景WC/GT|','|---|---|---|---:|---:|---:|']
for fam in families:
 for gene in ['16S','28S']:
  candidates=[x for x in allrows if x['name'].startswith(f'{fam}-{gene}-') and x['name'].endswith('D')]
  if candidates:
   x=max(candidates,key=lambda x:(int(x['covered']),int(x['exact']),-int(x['GT'])))
   summary.append(f"|{fam}|{gene}|{x['name']}|{x['covered']}/{x['target_species']} ({100*int(x['covered'])/int(x['target_species']):.1f}%)|{x['exact']}|{x['WC']}/{x['GT']}|")
  else:summary.append(f'|{fam}|{gene}|无通过候选|—|—|—|')
summary+=['','Syllidae的28S目前未取得实用的科级覆盖；低覆盖候选仅作搜索记录，不作为通用方案推荐。Nereididae的16S短片段与长片段存在覆盖和背景命中的取舍。不同基因参考库不同，不能将百分比直接解释为标记性能差异。','']
idx=lines.index('## 设计与核验');lines[idx:idx]=summary
write(R/'reference_coordinates.csv',coordinate_refs);write(R/'results.csv',allrows);write(R/'background_groups.csv',groups);write(R/'genera.csv',genera);(R/'primers.fasta').write_text(''.join(seqs));(R/'validation.json').write_text(json.dumps(validation,indent=2))
(R/'README.md').write_text('\n'.join(lines)+'\n')
subprocess.run(['pandoc','--standalone','--metadata','title=Family-specific primer candidates',str(R/'README.md'),'-o',str(R/'report.html')],check=True)
print('PASS',len(allrows),'pairs',len(seqs),'oligos')
