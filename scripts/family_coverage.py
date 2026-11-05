"""Render family-specific coverage bars from the released species counts."""
from pathlib import Path
import csv
from html import escape
root=Path(__file__).resolve().parents[1]
rows=list(csv.DictReader((root/'SML/families.csv').open()))
families=sorted({r['family'] for r in rows})
parts=['<svg xmlns="http://www.w3.org/2000/svg" width="1200" height="1150" viewBox="0 0 1200 1150">','<rect width="1200" height="1150" fill="white"/>','<style>text{font-family:Arial,sans-serif;fill:#172b27;font-size:14px}.small{font-size:12px}</style>']
def text(x,y,s,cls=''):
 parts.append(f'<text x="{x}" y="{y}" class="{cls}">{escape(s)}</text>')
colors=['#a4ccc5','#247c70','#a9c5e3','#396da1','#d2bee3','#80549f']
for g,gene in enumerate(['16S','28S']):
 x=25+g*600;text(x,28,gene+' · family-level species coverage')
 for k,family in enumerate(families):
  y=64+k*147;rr={r['name'].split('-')[-1]:r for r in rows if r['family']==family and f'-{gene}-' in r['name']}
  text(x,y,family+' (n = '+rr['S']['species']+')')
  for j,name in enumerate(['S','SD','M','MD','L','LD']):
   r=rr[name];n,c,e=map(int,[r['species'],r['covered'],r['exact']]);assert 0<=e<=c<=n
   yy=y+12+j*19;text(x,yy+11,name,'small')
   parts.append(f'<rect x="{x+38}" y="{yy}" width="350" height="13" fill="#edf1f0"/><rect x="{x+38}" y="{yy}" width="{350*c/n:.3f}" height="13" fill="{colors[j]}"/>')
   text(x+398,yy+11,f'{c}/{n} ({100*c/n:.1f}%)','small')
 text(x,1110,'Bar scale: 0–100% within each family','small')
text(25,1140,'Coverage: ≤2 mismatches per end; terminal 3 bases compatible. One reference per species; all panel species included.','small')
parts.append('</svg>');(root/'figures/family_coverage.svg').write_text('\n'.join(parts))
print('PASS: 84 family/primer counts checked; family coverage figure generated.')
