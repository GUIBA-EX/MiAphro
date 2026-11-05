from pathlib import Path
from collections import Counter,defaultdict
import csv,subprocess,shutil,concurrent.futures
from Bio import SeqIO
R=Path(__file__).resolve().parent;P=R.parent/'Phyllodocida_suborders/Nereidiformia'
read=lambda p:list(csv.DictReader(p.open()))
def run(job):
 fam,gene,band=job;root=R/fam;src=root/gene;out=root/'Controls'/gene/band;out.mkdir(parents=True,exist_ok=True)
 pairs={}
 for i,version in enumerate(['Expanded','Round2']):
  rows={x['name']:x for x in read(P/version/'results.csv')}
  for suffix in ['', 'D']:
   x=rows[f'Nereidiformia-{gene}-{band}{suffix}'];pairs[f'c{i}{suffix}']=x
 for role in ['target','background']:shutil.copyfile(src/(role+'.fa'),out/(role+'.fasta'))
 (out/'independent_primers.tsv').write_text(''.join(f"{k}\t{x['F']}\t{x['R']}\n" for k,x in pairs.items()))
 subprocess.run([str(root/'scripts/scan'),str(out)],stderr=subprocess.DEVNULL,check=True)
 bg={x['accession']:x for x in read(src/'background.csv')};near=sum(x['group']=='other Phyllodocida' for x in bg.values());n=len(bg);N=len(read(src/'target.csv'));stats=defaultdict(Counter)
 with (out/'independent_hits.tsv').open() as h:
  for x in csv.DictReader(h,delimiter='\t'):
   s=stats[x['name']]
   if x['role']=='target':s['covered']+=int(x['le2last3']);s['exact']+=int(x['exact'])
   else:
    m=x['model'];s[m]+=int(x['le2last3']);s[m+'loose']+=int(x['le2'] if m=='GT' else x['le3'])
    if bg[x['accession']]['group']=='other Phyllodocida':s['near'+m]+=int(x['le2last3'])
 fields=['id','reference','F_start','F_end','R_start','R_end','F','R','target_species','covered','exact','background_species','WC','GT','WCloose','GTloose','pass']
 with (out/'screen.csv').open('w') as h:
  w=csv.DictWriter(h,fields,lineterminator='\n');w.writeheader()
  for key,x in pairs.items():
   s=stats[key];passed=s['WC']<=int(n*.01) and s['GT']<=int(n*.03) and s['WCloose']<=int(n*.05) and s['GTloose']<=int(n*.10) and max(s['nearWC'],s['nearGT'])<=max(2,int(near*.03))
   z={k:x[k] for k in ['reference','F_start','F_end','R_start','R_end','F','R']};z.update(id=key,target_species=N,background_species=n,**{k:s[k] for k in ['covered','exact','WC','GT','WCloose','GTloose']});z['pass']=passed;w.writerow(z)
 print('Control',fam,gene,band,flush=True)
with concurrent.futures.ThreadPoolExecutor(max_workers=2) as pool:list(pool.map(run,[(f,g,b) for f in ['Syllidae','Nereididae','Hesionidae'] for g in ['16S','28S'] for b in 'SML']))
