from pathlib import Path
import csv
from collections import Counter
from Bio import SeqIO
from independent import flags
p=Path(__file__).resolve().parents[1];out=[];fr=[];fa=[]
for g,id in [('16S','p310'),('28S','p487')]:
 d=p/g;rr=list(csv.DictReader(open(d/'screen.csv')));panel={x['accession']:x for x in csv.DictReader(open(d/'target.csv'))}
 for suffix in ['','D']:
  r=next(x for x in rr if x['id']==id+suffix);name='MiAphro-'+g+suffix;hits=Counter();fam=Counter()
  for role in ['target','background']:
   for record in SeqIO.parse(d/(role+'.fa'),'fasta'):
    for gt in ([False] if role=='target' else [False,True]):
     e,c,b=flags(str(record.seq).upper(),r['F'],r['R'],role=='target',gt)
     if role=='target':
      hits['exact']+=e;hits['covered']+=c;k=panel[record.id]['family'];fam[k,'n']+=1;fam[k,'hit']+=c;fam[k,'exact']+=e
     else:
      k='GT' if gt else 'WC';hits[k]+=c;hits[k+'loose']+=b
  assert all(hits[k]==int(r[k]) for k in hits)
  out.append(dict(name=name,**r));fr.extend([name,k,fam[k,'n'],fam[k,'hit'],fam[k,'exact']] for k,typ in fam if typ=='n');fa.extend(f'>{name}_{s}\n{r[s]}\n' for s in 'FR')
outdir=p/'baseline';outdir.mkdir(exist_ok=True)
with open(outdir/'results.csv','w') as f:
 w=csv.DictWriter(f,out[0].keys());w.writeheader();w.writerows(out)
with open(outdir/'families.csv','w') as f:
 w=csv.writer(f);w.writerow(['name','family','species','covered','exact']);w.writerows(fr)
(outdir/'primers.fasta').write_text(''.join(fa))
