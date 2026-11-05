# MiPhyllo · 末端简并联合优化

[Current catalogue](../../README.md) · [Comparison data](comparison.json) · [Updated sequences](selected_results.csv) · [Previous PDF](pre_optimization_preprint.pdf)

共检查9对D：4对替换、5对保留。普通引物、长度和四个结合端点均不变；每条简并度≤8，简并仅位于3′半段。替换要求覆盖与完全匹配均不下降、至少一项提高，并通过原理化、背景及近缘筛选条件。

联合重分配允许先解除原简并，再同时分配新位点。种子包括原引物的具体展开序列及高频目标窗口后半段；最多3个3′半段位点联合分配简并，并枚举符合预算的末三位组合。覆盖优先与完全匹配优先短名单交叉配对，再扫原始序列。这是扩大的分阶段搜索，不构成全局最优证明。

|D pair|F/R candidates searched|Coverage|Exact|Strict WC/GT|Loose WC/GT|Terminal exclusions|Decision|
|---|---:|---:|---:|---|---|---:|---|
|Aphroditiformia-28S-MD|72,703/100,464|58→58/68|53→53/68|1/2→1/2|17/2→17/2|1→1|retain|
|Glyceriformia-16S-SD|63,957/43,835|14→15/20|1→1/20|14/33→13/33|113/104→107/105|1→0|replace|
|Nereididae-16S-LD|119,190/39,861|58→58/61|11→11/61|30/58→30/58|166/285→166/285|1→1|retain|
|Nereididae-28S-SD|96,492/80,430|29→30/52|8→19/52|0/0→0/0|0/0→0/0|7→5|replace|
|Phyllodociformia-16S-LD|107,108/99,790|40→41/59|3→3/59|0/1→0/0|76/29→80/31|2→2|replace|
|Phyllodociformia-16S-SD|4,160/56,040|49→50/59|11→11/59|0/8→0/12|81/11→82/15|1→0|replace|
|Phyllodociformia-28S-MD|11,689/86,668|52→52/63|5→5/63|21/26→21/26|95/42→95/42|5→5|retain|
|Phyllodociformia-28S-SD|54,053/15,417|51→51/63|4→4/63|7/23→7/23|53/53→53/53|5→5|retain|
|Syllidae-16S-SD|74,215/82,665|110→110/180|7→7/180|1/2→1/2|15/2→15/2|16→16|retain|

Syllidae 16S SD生成74,215个F和82,665个R分配候选；配对短名单中未找到覆盖和完全匹配均不下降且至少一项提高的方案。因此保留原序列，末端约束排除仍为16种。

“末端排除”表示同一扫描仅去掉末端兼容约束后增加的参考物种数，不表示实验扩增结果。所有计入覆盖的物种均有双端末三位兼容的匹配组合。

[Current terminal audit](current_terminal/summary.csv) · [Species-level terminal audit](current_terminal/species.csv) · [Updated species outcomes](selected_species.csv)

## Reproduction

Install biopython and primer3-py. `joint_search.py` regenerates expanded candidates using the frozen `baseline_catalog.csv`, the preserved initial `candidates.json` lists and repository reference panels. Remove only the generated `joint_pairs.json` files in a temporary copy to rerun search. `background.py` scans target/background panels with the compiled repository scanner, applies full gates and independently checks selected pairs in Python. `audit_terminal_current.py` audits the published current catalogue. Candidate-generation and screening details are retained per pair.

未新增实验结果。背景命中数不是扩增概率；部分替换方案的小幅背景计数增加已逐项列出。
