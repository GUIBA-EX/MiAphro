from pathlib import Path
from Bio import SeqIO
import sys
import csv,shutil,subprocess,collections,concurrent.futures,hashlib
root=Path(__file__).resolve().parents[1]
def run(gb):
 gene,band=gb;d=root/'SML'/gene/band;out=d/'marine';out.mkdir(exist_ok=True)
 digest=hashlib.sha256((root/gene/'marine_background.fa').read_bytes()).hexdigest()
 stamp=d/'marine_background.sha256'
 if stamp.exists() and stamp.read_text()==digest:return
 new={r.id:str(r.seq) for r in SeqIO.parse(root/gene/'marine_background.fa','fasta')}
 old={r.id:str(r.seq) for r in SeqIO.parse(out/'background.fasta','fasta')} if (out/'background.fasta').exists() else {}
 changed={k for k,v in new.items() if old.get(k)!=v}
 if changed:
  delta=out/'delta';delta.mkdir(exist_ok=True)
  shutil.copyfile(d/'independent_primers.tsv',delta/'independent_primers.tsv');(delta/'target.fasta').write_text('');(delta/'background.fasta').write_text(''.join('>'+k+'\n'+new[k]+'\n' for k in sorted(changed)))
  subprocess.run([str(root/'scripts/scan'),str(delta)],stderr=subprocess.DEVNULL,check=True)
  keep=[]
  if (out/'independent_hits.tsv').exists():
   with (out/'independent_hits.tsv').open() as h:
    next(h);keep=[line for line in h if line.split('\t')[3] in new and line.split('\t')[3] not in changed]
  with (delta/'independent_hits.tsv').open() as h:head=next(h);keep+=list(h)
  (out/'independent_hits.tsv').write_text(head+''.join(keep));shutil.rmtree(delta)
 shutil.copyfile(root/gene/'marine_background.fa',out/'background.fasta')
 shutil.copyfile(d/'independent_hits.tsv',d/'temporary_hits.tsv')
 rows=list(csv.DictReader((d/'screen.csv').open()));stats=collections.defaultdict(collections.Counter);bg=list(csv.DictReader((root/gene/'marine_background.csv').open()));groups={r['accession']:r['group'] for r in bg};den=collections.Counter(groups.values());by=collections.defaultdict(collections.Counter)
 with (out/'independent_hits.tsv').open() as h:
  for r in csv.DictReader(h,delimiter='\t'):
   key=r['name'];model=r['model'];c=int(r['le2last3']);l=int(r['le2'] if model=='GT' else r['le3']);stats[key][model]+=c;stats[key][model+'loose']+=l;by[key,groups[r['accession']]][model]+=c;by[key,groups[r['accession']]][model+'loose']+=l
 for r in rows:
  for k in ['WC','GT','WCloose','GTloose']:r[k]=stats[r['id']][k]
  r['background_species']=len(bg);r['pass']=max(r['WC'],r['GT'])<=10 and max(r['WCloose'],r['GTloose'])<=50 and all(max(by[r['id'],g]['WC'],by[r['id'],g]['GT'])<=max(2,int(n*.01)) for g,n in den.items())
 with (d/'screen.csv').open('w') as h:w=csv.DictWriter(h,rows[0]);w.writeheader();w.writerows(rows)
 with (d/'independent_hits.tsv').open('w') as h:
  with (d/'temporary_hits.tsv').open() as a:
   h.write(next(a));h.writelines(line for line in a if '\ttarget\t' in line)
  with (out/'independent_hits.tsv').open() as a:next(a);h.writelines(a)
 (d/'temporary_hits.tsv').unlink()
 detail=[dict(id=r['id'],group=g,total=n,**{k:by[r['id'],g][k] for k in ['WC','GT','WCloose','GTloose']}) for r in rows for g,n in den.items()]
 with (d/'marine_by_group.csv').open('w') as h:w=csv.DictWriter(h,detail[0]);w.writeheader();w.writerows(detail)
 shutil.copyfile(root/gene/'marine_background.fa',d/'background.fasta')
 stamp.write_text(digest)
 print('Marine scan',gene,band,'passing oligos pairs',sum(r['pass'] for r in rows),flush=True)
with concurrent.futures.ThreadPoolExecutor(max_workers=3) as pool:list(pool.map(run,[(g,b) for g in (sys.argv[1:] or ['16S','28S']) for b in 'SML']))
# Use only the marine design background in the final package.
for gene in ['16S','28S']:
 for ext in ['fa','csv']:
  d=root/gene;shutil.copyfile(d/('marine_background.'+ext),d/('background.'+ext))
