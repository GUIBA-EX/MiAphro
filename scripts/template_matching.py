"""Draw aligned primer/template letters for every window-evaluable species."""
from xml.etree import ElementTree as ET
from Bio.Seq import Seq
from Bio.Data.IUPACData import ambiguous_dna_values
from pathlib import Path
ROOT=Path(__file__).resolve().parents[1]
def element(g,tag,**attrs):
    return ET.SubElement(g,tag,{k.replace('_','-'):str(v) for k,v in attrs.items()})
def text(g,x,y,label,size,**kw):
    element(g,'text',x=x,y=y,font_size=size,**kw).text=label

OUT = ROOT/'figures'/'template_matching'
COLORS = {'match':'#d7edf4','mismatch':'#ffc895','unknown':'#e3e7e9'}


def render(gene,pair,records,alignment,cols,preview=False):
    records=sorted(records,key=lambda r:(r['family'],r['species'],r['accession']))
    if preview:
        families=set()
        selected=[]
        for r in records:
            if r['family'] not in families:
                selected.append(r)
                families.add(r['family'])
            if len(selected)==8:break
        records=selected
    width=1330
    n=len(records)
    height=129+21*n+16
    svg=ET.Element('svg',xmlns='http://www.w3.org/2000/svg',width=str(width),height=str(height),viewBox=f'0 0 {width} {height}')
    element(svg,'rect',width='100%',height='100%',fill='#f7fbff')
    g=element(svg,'g',font_family='Arial,sans-serif',fill='#173647')
    text(g,15,25,f"{pair['name']}  |  Primer–template base matching",21,font_weight='bold')
    text(g,15,50,f"{n} species shown of {pair['target_species']} panel species  ·  Blue: compatible  ·  Orange: mismatch  ·  Grey: gap/unknown  ·  Box: final 3 nt  ·  Primer-oriented 5′ → 3′",15)
    text(g,15,89,'Species / accession',15);text(g,358,89,'Raw hit',14)
    for side,x in [('F',438),('R',891)]:
        a,b=(int(pair[side+'_'+k]) for k in ('start','end'))
        text(g,x,87,f'{side}  {a}–{b}',16,font_weight='bold')
        q=pair[side]
        def bases(word,y,primer=False):
            paths={}
            for j,c in enumerate(word):
                status='unknown' if c not in 'ACGT' else 'match' if c in ambiguous_dna_values[q[j]] else 'mismatch'
                color='#c2e2eb' if primer else COLORS[status]
                paths.setdefault(color,[]).append(f'M{x+18*j},{y}h17v19h-17z')
            for color,path in paths.items():element(g,'path',d=''.join(path),fill=color)
            label=element(g,'text',x=' '.join(str(x+18*j+8.5) for j in range(len(word))),y=y+15,font_family='monospace',font_size=14,text_anchor='middle')
            label.text=word
        bases(q,97,True)
        for i,r in enumerate(records):
            y=125+21*i
            template=''.join(alignment[r['accession']][c] for c in cols[a-1:b])
            if side=='R':template=str(Seq(template).reverse_complement())
            bases(template,y)
            if side=='F':
                text(g,15,y+14,r['species']+' · '+r['accession'],12)
                text(g,358,y+14,'Exact' if r['exact']=='1' else 'Hit' if r['covered']=='1' else 'Miss',13)
        for top,h in [(97,19),(125,21*(n-1)+19)]:
            element(g,'rect',x=x+18*(len(q)-3),y=top,width=53,height=h,fill='none',stroke='#2a505f',stroke_width=1.4)
    suffix='_preview' if preview else ''
    ET.ElementTree(svg).write(OUT/f"{pair['name']}{suffix}.svg",encoding='utf-8',xml_declaration=True)


