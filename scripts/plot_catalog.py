from pathlib import Path
import csv
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
R=Path(__file__).resolve().parents[1];rows=list(csv.DictReader((R/'catalog.csv').open()));taxa=list(dict.fromkeys(x['target_taxon'] for x in rows));fig,axes=plt.subplots(1,2,figsize=(13,6),sharey=True)
colors={'S':'#078b8a','M':'#df941a','L':'#7352ab'};markers={'S':'o','M':'s','L':'^'}
for ax,gene in zip(axes,['16S','28S']):
 for band,offset in [('S',-.18),('M',0),('L',.18)]:
  rs=[x for x in rows if x['gene']==gene and x['band']==band and x['variant']=='D'];ax.scatter([100*int(x['covered'])/int(x['target_species']) for x in rs],[taxa.index(x['target_taxon'])+offset for x in rs],color=colors[band],marker=markers[band],label=band,s=50,zorder=3)
 ax.set_xlim(-2,102);ax.set_xticks([0,25,50,75,100]);ax.grid(axis='x',alpha=.2);ax.set_title(gene);ax.set_xlabel('D coverage within each reference panel (%)');ax.set_yticks(range(len(taxa)),taxa);ax.spines[['top','right']].set_visible(False)
axes[0].invert_yaxis();fig.legend(*axes[0].get_legend_handles_labels(),loc='upper center',bbox_to_anchor=(.55,.94),ncol=3,frameon=False);fig.suptitle('MiPhyllo · current candidates',fontsize=16);fig.tight_layout(rect=[0,0,1,.9]);(R/'figures').mkdir(exist_ok=True);fig.savefig(R/'figures/catalog_coverage.svg');fig.savefig(R/'figures/catalog_coverage.png',dpi=150);plt.close(fig)

p=R/'figures/catalog_coverage.svg';p.write_text('\n'.join(line.rstrip() for line in p.read_text().splitlines())+'\n')
