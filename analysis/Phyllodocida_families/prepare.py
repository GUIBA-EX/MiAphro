from pathlib import Path
import csv,shutil,json
from Bio import SeqIO
src=Path(__file__).resolve().parent.parent/'Phyllodocida_suborders/Nereidiformia';root=Path(__file__).resolve().parent
read=lambda p:list(csv.DictReader(p.open()))
def write(p,rows):
 with p.open('w') as h:
  w=csv.DictWriter(h,rows[0]);w.writeheader();w.writerows(rows)
summary=[]
for fam in ['Syllidae','Nereididae','Hesionidae']:
 d=root/fam;(d/'scripts').mkdir(parents=True,exist_ok=True)
 for name in ['chemistry_round2.py','independent.py','raw_scan.cpp','scan','design_round2.py','finish_round2.py']:
  shutil.copy2(src/'scripts'/name,d/'scripts'/name)
 p=d/'scripts/finish_round2.py';s=p.read_text().replace("(ROOT/'Expanded'/gene/band,OUT/gene/band)","(OUT/gene/band,ROOT/'Controls'/gene/band)").replace("f'Nereidiformia-{gene}-{band}{suffix}'",f"f'{fam}-{{gene}}-{{band}}{{suffix}}'");s=s.replace("passing=[v for v in pool.values() if all(r['pass']=='True' for r in v[:2])]", "passing=[v for v in pool.values() if all(r['pass']=='True' and int(r['covered'])>0 for r in v[:2])]");s=s.replace("dict(name=name,span=span,**r)","dict(name=name,span=span,source_pool=path.parents[1].name,**r)");p.write_text(s)
 for gene in ['16S','28S']:
  dest=d/gene;dest.mkdir(exist_ok=True);control=d/'Controls'/gene;control.mkdir(parents=True,exist_ok=True);shutil.copyfile(src/gene/'aligned.fa',control/'coordinate_reference_alignment.fa');rows=read(src/gene/'target.csv');target=[x for x in rows if x['family']==fam];other=[dict(**x,group='other Phyllodocida') for x in rows if x['family']!=fam];bg=read(src/gene/'background.csv')+other
  assert len({x['species'] for x in bg})==len(bg)
  write(dest/'target.csv',target);write(dest/'background.csv',bg)
  seq={x.id:x for x in SeqIO.parse(src/gene/'target.fa','fasta')};seq.update({x.id:x for x in SeqIO.parse(src/gene/'background.fa','fasta')})
  for role,rs in [('target',target),('background',bg)]:SeqIO.write([seq[x['accession']] for x in rs],dest/(role+'.fa'),'fasta')
  summary.append(dict(family=fam,gene=gene,target=len(target),background=len(bg)))
write(root/'reference_summary.csv',summary);print(json.dumps(summary,indent=2))
