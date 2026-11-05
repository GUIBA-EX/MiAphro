"""Build coverage heatmap and reference genus/species census from cached WoRMS mappings."""
from pathlib import Path
import csv
from html import escape
root=Path(__file__).resolve().parents[1]
def read(p):return list(csv.DictReader((root/p).open()))
refs=read('taxonomy/references.csv');mapping={(r['gene'],r['accession']):r for r in refs};families=sorted({r['worms_family'] for r in refs}|{'Iphionidae'})
def start(w,h):return [f'<svg xmlns="http://www.w3.org/2000/svg" width="{w}" height="{h}" viewBox="0 0 {w} {h}"><rect width="100%" height="100%" fill="white"/><style>text{{font-family:Arial,sans-serif;fill:#172b27;font-size:16px}}</style>']
def txt(p,x,y,s):p.append(f'<text x="{x}" y="{y}">{escape(s)}</text>')
p=start(1240,540);stats=[]
for gi,gene in enumerate(['16S','28S']):
 x=190+gi*520;txt(p,x,30,gene+' · covered / reference species')
 for j,n in enumerate(['S','SD','M','MD','L','LD']):txt(p,x+j*83+27,65,n)
 for i,fam in enumerate(families):
  y=85+i*58
  if gi==0:txt(p,12,y+30,fam)
  for j,n in enumerate(['S','SD','M','MD','L','LD']):
   rr=[r for r in read('SML/species.csv') if r['name']==f'MiAphro-{gene}-{n}' and mapping[gene,r['accession']]['worms_family']==fam];total=len(rr);c=sum(int(r['covered']) for r in rr);v=c/total if total else 0
   rgb=tuple(round(a+(b-a)*v) for a,b in zip((239,246,244),(43,136,120)));fill='#%02x%02x%02x'%rgb if total else '#dddddd'
   p.append(f'<rect x="{x+j*83}" y="{y}" width="79" height="54" fill="{fill}"/>');txt(p,x+j*83+12,y+23,f'{100*v:.0f}%' if total else 'NA');txt(p,x+j*83+12,y+44,f'{c}/{total}')
   stats.append(dict(gene=gene,primer=n,family=fam,covered=c,species=total))
txt(p,12,470,'Color scale: 0% (pale) to 100% (dark teal); grey = no reference species.');txt(p,12,500,'Families follow WoRMS genus assignments (5 November 2026). Terminal triplet compatible; ≤2 mismatches/end.')
p.append('</svg>');(root/'figures/family_coverage.svg').write_text('\n'.join(p))
p=start(1050,500);txt(p,15,30,'Reference sampling by WoRMS family');labels=['16S genera','16S species','28S genera','28S species','Union genera','Union species']
for j,label in enumerate(labels):txt(p,200+j*138,70,label)
census=[]
for i,fam in enumerate(families):
 y=95+i*53;txt(p,15,y+27,fam);values=[]
 for gene in ['16S','28S','Union']:
  rr=[r for r in refs if r['worms_family']==fam and (gene=='Union' or r['gene']==gene)];g=len({r['worms_genus'] for r in rr});s=len({r['species'] for r in rr});values += [g,s];census.append(dict(family=fam,panel=gene,genera=g,species=s,genera_names='; '.join(sorted({r['worms_genus'] for r in rr}))))
 for j,v in enumerate(values):
  p.append(f'<rect x="{195+j*138}" y="{y}" width="131" height="47" fill="'+('#e3efeb' if j%2==0 else '#edf2f8')+'"/>');txt(p,245+j*138,y+28,str(v))
txt(p,15,450,'Union counts unique taxa across markers. These are sampled counts, not percentages of all WoRMS taxa.');txt(p,15,477,'Genus and family assignments checked against WoRMS; species names retain source-record binomials.')
p.append('</svg>');(root/'figures/reference_taxonomy.svg').write_text('\n'.join(p))
for name,rows in [('family_coverage.csv',stats),('reference_coverage.csv',census)]:
 with (root/'taxonomy'/name).open('w') as h:
  w=csv.DictWriter(h,rows[0]);w.writeheader();w.writerows(rows)
print('Generated WoRMS family coverage and reference census.')
