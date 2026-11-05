"""Plot primer coverage and genus sampling using the cached WoRMS census."""
from pathlib import Path
import csv
from html import escape
root=Path(__file__).resolve().parents[1]
def read(p):return list(csv.DictReader((root/p).open()))
refs=read('taxonomy/references.csv');genera=read('taxonomy/worms_accepted_genera.csv');families=sorted({r['family'] for r in genera});mapping={(r['gene'],r['accession']):r for r in refs}
def start(w,h):return [f'<svg xmlns="http://www.w3.org/2000/svg" width="{w}" height="{h}" viewBox="0 0 {w} {h}"><rect width="100%" height="100%" fill="white"/><style>text{{font-family:Arial,sans-serif;fill:#172b27;font-size:17px}}.small{{font-size:14px}}</style>']
def txt(p,x,y,s,small=False):p.append(f'<text x="{x}" y="{y}"'+(' class="small"' if small else '')+f'>{escape(str(s))}</text>')
def rect(p,x,y,w,h,c):p.append(f'<rect x="{x}" y="{y}" width="{w}" height="{h}" rx="3" fill="{c}"/>')
def col(v):return '#%02x%02x%02x'%tuple(round(a+(b-a)*v) for a,b in zip((241,247,245),(81,166,145)))
p=start(1230,535);stats=[];species=read('SML/species.csv')
for gi,gene in enumerate(['16S','28S']):
 x=190+gi*520;txt(p,x,30,gene+' · predicted species coverage')
 for j,n in enumerate(['S','SD','M','MD','L','LD']):txt(p,x+j*83+25,65,n)
 for i,fam in enumerate(families):
  y=85+i*55
  if gi==0:txt(p,12,y+31,fam)
  for j,n in enumerate(['S','SD','M','MD','L','LD']):
   rr=[r for r in species if r['name']==f'MiAphro-{gene}-{n}' and mapping[gene,r['accession']]['worms_family']==fam];total=len(rr);c=sum(int(r['covered']) for r in rr);v=c/total if total else 0
   rect(p,x+j*83,y,79,51,col(v) if total else '#e4e4e4');txt(p,x+j*83+12,y+22,f'{100*v:.0f}%' if total else 'NA');txt(p,x+j*83+12,y+42,f'{c}/{total}',True)
   stats.append(dict(gene=gene,primer=n,family=fam,covered=c,species=total))
for i in range(101):rect(p,190+i*3,440,3.1,14,col(i/100))
txt(p,185,477,'0%',True);txt(p,460,477,'100%',True);txt(p,530,453,'Grey: no reference species. Cells: coverage and covered/total.',True)
txt(p,12,512,'WoRMS family assignments · ≤2 mismatches/end; compatible terminal triplets · All panel references included',True)
p.append('</svg>');(root/'figures/family_coverage.svg').write_text('\n'.join(p))
p=start(1230,570);census=[];colors=['#368e83','#497cba','#8056a0']
for j,label in enumerate(['16S','28S','Union']):rect(p,210+j*210,16,22,13,colors[j]);txt(p,242+j*210,29,label)
for i,fam in enumerate(families):
 y=65+i*70;den={r['aphia_id'] for r in genera if r['family']==fam};txt(p,12,y+18,fam);txt(p,12,y+39,f'{len(den)} accepted genera',True)
 for j,gene in enumerate(['16S','28S','Union']):
  rr=[r for r in refs if r['worms_family']==fam and (gene=='Union' or r['gene']==gene)];ids={r['aphia_id'] for r in rr};assert ids<=den,(fam,ids-den)
  count=len(ids);sp=len({r['species'] for r in rr});pct=100*count/len(den);yy=y+j*19;rect(p,215,yy,690,13,'#edf0f2');rect(p,215,yy,690*pct/100,13,colors[j]);txt(p,925,yy+12,f'{count}/{len(den)} genera ({pct:.1f}%)',True)
  census.append(dict(family=fam,panel=gene,genera=count,worms_genera=len(den),genus_coverage_percent=round(pct,3),species=sp,genera_names='; '.join(sorted({r['worms_genus'] for r in rr}))))
for v in [0,25,50,75,100]:txt(p,210+690*v/100,510,f'{v}%',True)
txt(p,12,544,'Accepted genus rank, all habitat/extinction flags; synonyms and subgenera excluded. WoRMS snapshot: 5 November 2026.',True)
p.append('</svg>');(root/'figures/reference_taxonomy.svg').write_text('\n'.join(p))
for name,rows in [('family_coverage.csv',stats),('reference_coverage.csv',census)]:
 with (root/'taxonomy'/name).open('w') as h:w=csv.DictWriter(h,rows[0],lineterminator='\n');w.writeheader();w.writerows(rows)
print('PASS: all sampled genus IDs occur in the accepted genus census.')
