"""Rank shared-endpoint pairs on the frozen panel and independently check selected hits."""
import csv
from collections import Counter
from pathlib import Path
from Bio import SeqIO
from independent import flags
ROOT=Path(__file__).resolve().parents[1]
OUT=ROOT/'Expanded'
def read(path):
 with path.open() as f:return list(csv.DictReader(f))
def write(name,rows):
 if not rows:(OUT/name).write_text('name\n');return
 with (OUT/name).open('w') as h:
  w=csv.DictWriter(h,rows[0],lineterminator='\n');w.writeheader();w.writerows(rows)
alternatives=[];results=[];families=[];details=[];summary=[];sequence=[];text=[]
for gene in ('16S','28S'):
 src=ROOT/gene
 panel={x['accession']:x for x in read(src/'target.csv')}
 libs={role:{r.id:str(r.seq).upper() for r in SeqIO.parse(src/(role+'.fa'),'fasta')} for role in ['target','background']}
 for band,(low,high) in zip('SML',[(100,200),(201,300),(301,480)]):
  pool={}
  for path in (OUT/gene/band,):
   rs=read(path/'screen.csv');lookup={r['id']:r for r in rs}
   for a in rs:
    if a['id'].endswith('D'):continue
    b=lookup[a['id']+'D'];span=int(a['R_end'])-int(a['F_start'])+1
    if not low<=span<=high:continue
    key=tuple(a[k] for k in ('F_start','F_end','R_start','R_end'))
    assert all(a[k]==b[k] for k in ('F_start','F_end','R_start','R_end'))
    pool[key+(a['F'],a['R'],b['F'],b['R'])]=a,b,path,span
  passing=[v for v in pool.values() if all(r['pass']=='True' for r in v[:2])]
  summary.append(dict(gene=gene,band=band,windows=len({k[:4] for k in pool}),passing=len({tuple(v[0][k] for k in ("F_start","F_end","R_start","R_end")) for v in passing})))
  if not passing:
   text.append(f'| {gene} | {band} | — | 未通过 | — | — | — |');continue
  def rank(v):
   a,b=v[:2]
   return (int(b['covered']),int(a['covered']),sum(int(r['exact']) for r in [a,b]),-max(int(r['GT']) for r in [a,b]))
  for objective,measure in [('coverage',lambda v:sum(int(r['covered']) for r in v[:2])),('exact',lambda v:sum(int(r['exact']) for r in v[:2]))]:
   x,y,_,length=max(passing,key=lambda v:(measure(v),rank(v)))
   for row in (x,y):alternatives.append(dict(gene=gene,band=band,objective=objective,span=length,**row))
  a,b,path,span=max(passing,key=rank)
  with (path/'independent_hits.tsv').open() as h:
   baseline={(x['name'],x['role'],x['model'],x['accession']):x for x in csv.DictReader(h,delimiter='\t') if x['name'] in [a['id'],b['id']]}
  for r,suffix in [(a,''),(b,'D')]:
   name=f'Nereidiformia-{gene}-{band}{suffix}';counts=Counter();fam=Counter()
   for role,lib in libs.items():
    for acc,s in lib.items():
     for gt in ([False] if role=='target' else [False,True]):
      model='GT' if gt else 'WC';e,c,l=flags(s,r['F'],r['R'],role=='target',gt)
      x=baseline.get((r['id'],role,model,acc),{})
      assert [e,c,l]==[x.get(k)=='1' for k in ['exact','le2last3','le3' if role=='background' and not gt else 'le2']]
      if role=='target':
       f=panel[acc]['family'];fam[f,'species']+=1;fam[f,'covered']+=c;fam[f,'exact']+=e
       details.append(dict(name=name,**panel[acc],covered=int(c),exact=int(e)))
       counts['covered']+=c;counts['exact']+=e
      else:counts[model]+=c;counts[model+'loose']+=l
   assert all(int(r[k])==v for k,v in counts.items())
   families.extend(dict(name=name,family=f,**{k:fam[f,k] for k in ('species','covered','exact')}) for f,k in fam if k=='species')
   results.append(dict(name=name,span=span,**r));sequence.extend(f'>{name}_{s}\n{r[s]}\n' for s in 'FR')
  text.append(f"| {gene} | {band} | {span} | {a['covered']} / {b['covered']} | {a['exact']} / {b['exact']} | {a['WC']},{a['GT']} / {b['WC']},{b['GT']} | {a['WCloose']},{a['GTloose']} / {b['WCloose']},{b['GTloose']} |")
write('alternative_objectives.csv',alternatives);write('results.csv',results);write('families.csv',families);write('species.csv',details);write('search.csv',summary)
(OUT/'primers.fasta').write_text(''.join(sequence))
