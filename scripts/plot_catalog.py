from pathlib import Path
import csv
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D
R=Path(__file__).resolve().parents[1];rows=list(csv.DictReader((R/'catalog.csv').open()));taxa=list(dict.fromkeys(x['target_taxon'] for x in rows))
fig,axes=plt.subplots(1,2,figsize=(12,5.7),sharey=True)
colors={'S':'#078b8a','M':'#c38213','L':'#7352ab'};markers={'S':'o','M':'s','L':'^'}
for ax,gene in zip(axes,['16S','28S']):
 for band,offset in [('S',-.23),('M',0),('L',.23)]:
  rs=[x for x in rows if x['gene']==gene and x['band']==band and x['variant']=='D']
  for x in rs:
   y=taxa.index(x['target_taxon'])+offset;c=100*int(x['covered'])/int(x['target_species']);e=100*int(x['exact'])/int(x['target_species'])
   ax.plot([e,c],[y,y],color='#bbc6cc',lw=1.2,zorder=1)
   ax.scatter(c,y,color=colors[band],marker=markers[band],s=46,zorder=3)
   ax.scatter(e,y,facecolor='white',edgecolor=colors[band],marker=markers[band],s=46,zorder=3)
 ax.axvline(50,color='#c5ced1',ls=':',lw=1);ax.set_xlim(-3,103);ax.set_xticks([0,25,50,75,100]);ax.grid(axis='x',alpha=.15);ax.set_title(gene,fontweight='bold');ax.set_xlabel('Reference species (%)');ax.set_yticks(range(len(taxa)),taxa);ax.spines[['top','right']].set_visible(False)
axes[0].invert_yaxis()
handles=[Line2D([],[],color=colors[b],marker=markers[b],linestyle='none',label=b) for b in colors]+[Line2D([],[],color='#4b5861',marker='o',linestyle='none',label='D coverage'),Line2D([],[],color='#4b5861',marker='o',markerfacecolor='white',linestyle='none',label='Exact match')]
fig.legend(handles=handles,loc='upper center',bbox_to_anchor=(.57,.92),ncol=5,frameon=False)
fig.suptitle('MiPhyllo · curated ordinary/D sets',fontsize=16)
fig.text(.53,.015,'24 retained sets · missing lengths archived or unavailable · each panel has its own denominator',ha='center',fontsize=9,color='#4b5861')
fig.tight_layout(rect=[0,.05,1,.86]);fig.savefig(R/'figures/catalog_coverage.svg');fig.savefig(R/'figures/catalog_coverage.png',dpi=180);plt.close(fig)
p=R/'figures/catalog_coverage.svg';p.write_text('\n'.join(x.rstrip() for x in p.read_text().splitlines())+'\n')
