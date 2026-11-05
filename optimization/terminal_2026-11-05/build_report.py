from pathlib import Path
import json,csv
O=Path(__file__).resolve().parent;R=O.parents[1];cs=json.loads((O/'comparison.json').read_text());old={x['name']:x for x in csv.DictReader((O/'terminal_audit_baseline.csv').open())};new={x['name']:x for x in csv.DictReader((O/'current_terminal/summary.csv').open())}
for en in [False,True]:
 f='README.en.md' if en else 'README.md'
 lines=['# MiPhyllo · Joint terminal degeneracy optimization' if en else '# MiPhyllo · 末端简并联合优化','',f'[Current catalogue](../../{f}) · [Comparison data](comparison.json) · [Updated sequences](selected_results.csv) · [Previous PDF](pre_optimization_preprint.pdf)','',
 'Nine D pairs were evaluated. Four were replaced; five were retained. Ordinary sequences, lengths and all four endpoints are unchanged. Each oligo has degeneracy ≤8, restricted to its 3′ half. Replacement requires no loss in target coverage or exact matching and improvement in at least one, plus the original thermodynamic, background and close-relative gates.' if en else '共检查9对D：4对替换、5对保留。普通引物、长度和四个结合端点均不变；每条简并度≤8，简并仅位于3′半段。替换要求覆盖与完全匹配均不下降、至少一项提高，并通过原理化、背景及近缘筛选条件。','',
 'Joint redistribution resolves previous ambiguities and reallocates them simultaneously. Seeds include canonical expansions and frequent target-window suffixes. Up to three 3′-half positions are jointly assigned masks; terminal-triplet masks are also enumerated within budget. Coverage- and exact-oriented shortlists are combined and re-scanned on raw sequences. This is an expanded staged search, not a global optimum proof.' if en else '联合重分配允许先解除原简并，再同时分配新位点。种子包括原引物的具体展开序列及高频目标窗口后半段；最多3个3′半段位点联合分配简并，并枚举符合预算的末三位组合。覆盖优先与完全匹配优先短名单交叉配对，再扫原始序列。这是扩大的分阶段搜索，不构成全局最优证明。','',
 '|D pair|F/R candidates searched|Coverage|Exact|Strict WC/GT|Loose WC/GT|Terminal exclusions|Decision|','|---|---:|---:|---:|---|---|---:|---|']
 for x in cs:
  b,a=x['before'],x['after'];n=json.loads((O/x['name']/'search_size.json').read_text());lines.append(f"|{x['name'].removeprefix('MiPhyllo-')}|{n['F']:,}/{n['R']:,}|{b['covered']}→{a['covered']}/{x['target_N']}|{b['exact']}→{a['exact']}/{x['target_N']}|{b['WC']}/{b['GT']}→{a['WC']}/{a['GT']}|{b['WCloose']}/{b['GTloose']}→{a['WCloose']}/{a['GTloose']}|{old[x['name']]['terminal_excluded']}→{new[x['name']]['terminal_excluded']}|{x['status']}|")
 lines+=['','Syllidae 16S SD: 74,215 forward and 82,665 reverse allocations were generated. No shortlisted pair improved coverage/exact matching without a loss in the other metric. Its sequence is retained; terminal mismatch exclusions remain 16 species.' if en else 'Syllidae 16S SD生成74,215个F和82,665个R分配候选；配对短名单中未找到覆盖和完全匹配均不下降且至少一项提高的方案。因此保留原序列，末端约束排除仍为16种。','',
 '“Terminal exclusions” counts reference species gained by removing only terminal compatibility from the same scan; it is not a measure of experimental amplification. Every species counted as covered retains compatible terminal three bases at both ends.' if en else '“末端排除”表示同一扫描仅去掉末端兼容约束后增加的参考物种数，不表示实验扩增结果。所有计入覆盖的物种均有双端末三位兼容的匹配组合。','',
 '[Current terminal audit](current_terminal/summary.csv) · [Species-level terminal audit](current_terminal/species.csv) · [Updated species outcomes](selected_species.csv)','',
 '## Reproduction','',
 'Install biopython and primer3-py. `joint_search.py` regenerates expanded candidates using the frozen `baseline_catalog.csv`, the preserved initial `candidates.json` lists and repository reference panels. Remove only the generated `joint_pairs.json` files in a temporary copy to rerun search. `background.py` scans target/background panels with the compiled repository scanner, applies full gates and independently checks selected pairs in Python. `audit_terminal_current.py` audits the published current catalogue. Candidate-generation and screening details are retained per pair.','',
 'No experimental results are added. Background counts are not amplification probabilities; modest increases in some accepted candidates are shown explicitly.' if en else '未新增实验结果。背景命中数不是扩增概率；部分替换方案的小幅背景计数增加已逐项列出。','']
 (O/f).write_text('\n'.join(lines))
# Figure drawn only from numerical comparison records.
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
fig,axs=plt.subplots(1,2,figsize=(12,6),sharey=True)
for ax,key,title in zip(axs,['covered','exact'],['Target coverage','Exact matching']):
 for i,x in enumerate(cs):
  b=100*x['before'][key]/int(x['target_N']);a=100*x['after'][key]/int(x['target_N']);ax.plot([b,a],[i,i],color='#8aadb1',lw=2);ax.scatter(b,i,facecolor='white',edgecolor='#758991',s=55,zorder=3);ax.scatter(a,i,color='#007c83' if x['status']=='replace' else '#9ba9ac',s=24,zorder=4)
 ax.set_title(title);ax.set_xlim(0,105);ax.set_xlabel('Reference species (%)');ax.set_yticks(range(len(cs)),[x['name'].removeprefix('MiPhyllo-') for x in cs]);ax.grid(axis='x',alpha=.15);ax.spines[['right','top']].set_visible(False)
axs[0].invert_yaxis();fig.suptitle('MiPhyllo · terminal optimization',fontsize=16);fig.text(.6,.02,'Open: before   •   Filled: after   •   Teal: accepted replacement',ha='center',fontsize=10);fig.tight_layout(rect=[0,.05,1,.94]);fig.savefig(R/'figures/terminal_optimization.svg');fig.savefig(R/'figures/terminal_optimization.png',dpi=160)
p=R/'figures/terminal_optimization.svg';p.write_text('\n'.join(z.rstrip() for z in p.read_text().splitlines())+'\n')
