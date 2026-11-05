from pathlib import Path
import csv,json,re,subprocess,hashlib,io
from Bio import SeqIO
R=Path(__file__).resolve().parent;O=R/'Optimization3'
read=lambda p:list(csv.DictReader(p.open()))
rows=read(O/'results.csv');data={(x['objective'],x['name']):x for x in rows};old={x['name']:x for fam in ['Syllidae','Nereididae'] for x in read(R/fam/'Round2/results.csv')}
from collections import Counter
genera=Counter()
for x in read(O/'species.csv'):
 key=(x['objective'],x['name'],x['species'].split()[0]);genera[key,'total']+=1;genera[key,'covered']+=int(x['covered']);genera[key,'exact']+=int(x['exact'])
with (O/'genera.csv').open('w') as h:
 w=csv.writer(h);w.writerow(['objective','name','genus','total','covered','exact'])
 for (key,metric),n in genera.items():
  if metric=='total':w.writerow([*key,n,genera[key,'covered'],genera[key,'exact']])
search=[]
for fam,gene in [('Syllidae','16S'),('Nereididae','28S')]:
 for band in 'SML':
  text=(R/fam/f'expanded3_{gene}_{band}.log').read_text();m=re.search(r'coordinate windows (\d+) full combinations windows (\d+)',text);n=re.search(r'selected windows (\d+)',text);assert m and n
  search.append(dict(family=fam,gene=gene,band=band,coordinate_windows=int(m[1]),full_combination_windows=int(m[2]),scanned_windows=int(n[1])))
(O/'search_configuration.json').write_text(json.dumps(dict(primer_lengths=list(range(18,31)),start_step=1,degeneracy_per_oligo_max=8,degeneracy_region='3prime_half',target_background_hashes_verified=True,search=search),indent=2))
from collections import defaultdict
byfile=defaultdict(list);background_rows=[]
for x in rows:
 fam,gene,band=x['name'].split('-');band=band.removesuffix('D');byfile[R/fam/x['source_pool']/gene/band/'independent_hits.tsv'].append(x)
for path,selected in byfile.items():
 fam=selected[0]['name'].split('-')[0];gene=selected[0]['name'].split('-')[1];bg={x['accession']:x for x in read(R/fam/gene/'background.csv')};totals=Counter(x['group'] for x in bg.values());ids={x['id'] for x in selected}
 with path.open() as h:header=h.readline()
 match=subprocess.run(['rg','--no-filename','^('+ '|'.join(re.escape(i) for i in sorted(ids))+r')\t',str(path)],capture_output=True,text=True);assert match.returncode in [0,1]
 counts=Counter()
 for z in csv.DictReader(io.StringIO(header+match.stdout),delimiter='\t'):
  if z['role']=='background':counts[z['name'],bg[z['accession']]['group'],z['model']]+=int(z['le2last3'])
 for x in selected:
  for model in ['WC','GT']:assert sum(counts[x['id'],g,model] for g in totals)==int(x[model])
  assert max(counts[x['id'],'other Phyllodocida',model] for model in ['WC','GT'])<=max(2,int(totals['other Phyllodocida']*.03))
  for group,total in totals.items():background_rows.append(dict(objective=x['objective'],name=x['name'],group=group,total=total,WC=counts[x['id'],group,'WC'],GT=counts[x['id'],group,'GT']))
with (O/'background_groups.csv').open('w') as h:
 w=csv.DictWriter(h,background_rows[0]);w.writeheader();w.writerows(background_rows)
lines=['# Syllidae 16S与Nereididae 28S：扩大搜索','', '作者：Fei Xia','', '沿用同一目标参考与海洋背景，分别搜索S 100–200 bp、M 201–300 bp、L 301–480 bp。普通与D四个结合端点一致，每条简并度≤8，简并仅限3′半段。所有最终列出的均衡、覆盖优先、完全匹配优先与严格匹配限制备选方案均经过独立算法逐记录复核。','', '## 主要结果', '', '原背景阈值下，Syllidae 16S三档均保留现版；Nereididae 28S SD覆盖从28增至29/52，完全匹配保持8种，LD由2增至4/52但覆盖仍低。条件备选中的Syllidae SD达到126/180覆盖、48/180完全匹配，严格WC/GT背景33/42，宽松WC/GT背景1007/601（背景共3675种）；这代表显著的目标匹配收益，同时对允许错配的模型更敏感。', '', '## 选择原则','', '**均衡方案**要求普通与D的覆盖数、完全匹配数四项均不低于原方案，再依次按D覆盖、D完全匹配、普通覆盖、普通完全匹配与较少GT背景命中排序。旧方案一同参与选择，因此没有整体改善时保留旧方案。另列覆盖优先与完全匹配优先方案，展示指标间的取舍，不自动替换均衡方案。','', '## 均衡方案与现版比较','','|名称|覆盖：旧→新|完全匹配：旧→新|严格背景WC/GT：旧→新|来源|','|---|---:|---:|---:|---|']
for x in rows:
 if x['objective']!='balanced':continue
 a=old[x['name']];n=int(x['target_species']);lines.append(f"|{x['name']}|{a['covered']}→{x['covered']}/{n}|{a['exact']}→{x['exact']}/{n}|{a['WC']}/{a['GT']}→{x['WC']}/{x['GT']}|{x['source_pool']}|")
lines+=['','## D方案的单项指标优先与条件备选','','下表覆盖和完全匹配均以全部目标参考物种为分母。strict_profile保留严格WC≤1%、GT≤3%及近缘非目标门槛，要求普通/D覆盖与完全匹配四项不退步，但不再以宽松模型命中作为硬性否决；宽松命中较高时不能将其概括为普遍的低非目标风险。普通配套方案及其完整指标见结果CSV。','','|名称|策略|覆盖|完全匹配|严格背景WC/GT|宽松背景WC/GT|','|---|---|---:|---:|---:|---:|']
for x in rows:
 if x['objective']=='balanced' or not x['name'].endswith('D'):continue
 lines.append(f"|{x['name']}|{x['objective']}|{x['covered']}/{x['target_species']}|{x['exact']}/{x['target_species']}|{x['WC']}/{x['GT']}|{x['WCloose']}/{x['GTloose']}|")
lines+=['','## 均衡方案序列（5′→3′）','','|名称|F|R|指定参考坐标F；R|','|---|---|---|---|']
for x in rows:
 if x['objective']=='balanced':lines.append(f"|{x['name']}|`{x['F']}`|`{x['R']}`|{x['reference']}：{x['F_start']}–{x['F_end']}；{x['R_start']}–{x['R_end']}|")
lines+=['','## 搜索方法','', '引物长度由18/20/22/24/26/28 nt扩展为18–30 nt逐碱基长度，结合起点由2 nt步长加密为1 nt。每个位点使用最多16个高频单倍型及共识种子；简并分支同时保留覆盖与完全匹配优先候选。先对所有满足S/M/L跨度的坐标组合用代表序列快速评分，再综合覆盖、完全匹配、空间分布与上一轮窗口进行完整组合评估，每档最多3000个窗口接受热力学与全背景扫描。该分阶段方法扩大了搜索，但未穷举所有引物序列组合。', '', '同一坐标的寡核苷酸候选保留覆盖和完全匹配两种排序；进入背景扫描的首选配对仍优先选取热力学合格的高覆盖组合。单项指标优先结果是已扫描候选池中的最优，不代表序列空间的全局最优。','', '目标覆盖为双端各≤2错配、末端3位兼容、产物100–500 bp且内部≥50 bp；完全匹配为双端零错配。原阈值方案背景严格WC≤1%、GT≤3%，宽松WC≤5%、GT≤10%；strict_profile取消后两项宽松模型硬阈值，其余阈值不变。其他Phyllodocida严格WC和GT各≤max(2, floor(3%×近缘参考数))。严格模型要求末端3位兼容，宽松WC允许每端≤3错配且不约束末端，宽松GT允许≤2错配且不约束末端。背景允许内部间隔0–1000 bp。GT为敏感性模型，其命中不是实测非目标扩增概率。','', '热力学与末端规则保持不变：50 nM/变体和200 nM/方向两种浓度情景，最低Tm≥45°C、平均Tm差≤5°C、二级结构Tm<55°C；末端5位至少一个确定G/C，排除可能的末端3连聚与全序列5连聚。坐标是指定参考按比对方向去gap后的1-based坐标；跨参考不能直接比较。','', '## 局部候选细化', '', '每档另外挑选最多24个有希望的坐标，使用最多32个高频种子及共识，并对普通种子尝试两轮向位点主要碱基的单碱基改进。简并候选分别保留覆盖与完全匹配前列分支；对每一类候选、每一目标保留最多5个热力学合格配对进行背景扫描。最终允许同一坐标内通过的任意普通方案与任意D方案组合，减少单次预选配对造成的遗漏。', '', '## 文件','', '[均衡序列](primers.fasta) · [覆盖优先序列](coverage_primers.fasta) · [完全匹配优先序列](exact_primers.fasta) · [严格匹配限制备选序列](strict_profile_primers.fasta) · [全部指标](results.csv) · [严格背景分组](background_groups.csv) · [逐物种匹配](species.csv) · [逐属覆盖](genera.csv) · [搜索汇总](search_summary.csv) · [搜索配置与窗口数](search_configuration.json) · [独立复核](validation.json)','', '来源标记：Refine3为本轮局部细化，Expanded3为本轮扩大搜索，Round2为上一轮科级搜索，Controls为原亚目候选在当前科级面板上的重测。运行`run_expanded3.py`完成搜索，随后运行`refine_expanded3.py`；再依次运行`finish_expanded3.py`和`report_expanded3.py`完成独立复核及报告。未添加实验验证。']
(O/'README.md').write_text('\n'.join(lines)+'\n');subprocess.run(['pandoc','--standalone','--metadata','title=Expanded family primer search',str(O/'README.md'),'-o',str(O/'report.html')],check=True)
html=(O/'report.html').read_text();html=html.replace('</head>','<style>body{max-width:1450px;margin:32px auto;padding:0 24px;font:15px/1.65 system-ui,sans-serif;color:#183247}h1,h2{color:#124d62}table{display:block;overflow-x:auto;width:100%;border-collapse:collapse;margin:20px 0}th,td{padding:8px 12px;border-bottom:1px solid #dce5e9;text-align:left}th{background:#eef5f7}code{white-space:nowrap;font-size:13px}a{color:#086a87}</style></head>');(O/'report.html').write_text(html)
for objective in ['balanced','coverage','exact','strict_profile']:
 path=O/('primers.fasta' if objective=='balanced' else objective+'_primers.fasta');seq={x.id:str(x.seq) for x in SeqIO.parse(path,'fasta')};rs=[x for x in rows if x['objective']==objective];assert len(seq)==len(rs)*2
 for x in rs:
  for e in 'FR':assert seq[x['name']+'_'+e]==x[e]
for link in re.findall(r'\]\(([^)]+)\)',(O/'README.md').read_text()):assert (O/link).exists()
(O/'SHA256SUMS').write_text(''.join(hashlib.sha256(p.read_bytes()).hexdigest()+'  '+p.name+'\n' for p in sorted(O.iterdir()) if p.is_file() and p.name!='SHA256SUMS'))
print('PASS: FASTA, links, report and checksums')
