"""Finite one-base search for low-margin pairs; no catalog replacements.

Both 3-prime triplets are fixed. Enumerates a single substitution in one end,
retaining the existing other end. Search is heuristic, not exhaustive.
"""
import csv,itertools,json,subprocess,tempfile,shutil,sys
from pathlib import Path
from collections import Counter
from Bio import SeqIO
ROOT=Path(__file__).resolve().parents[1]
sys.path.insert(0,str(ROOT/'optimization/terminal_2026-11-05'))
sys.path.insert(0,str(ROOT/'scripts'))
from chemistry_rules import thermal, oligo_thermal, terminal
from independent import flags

def read(p):
 with p.open() as h:return list(csv.DictReader(h))

def main():
 catalog=read(ROOT/'catalog.csv'); structures=read(ROOT/'audit/review_20261107/all_pair_structures.csv')
 selected={r['name'] for r in structures if float(r['margin'])<5}
 out=ROOT/'audit/review_20261107';summary=[];allc=[]
 with tempfile.TemporaryDirectory(prefix='miphyllo-lowmargin-') as tmp:
  tmp=Path(tmp);binary=tmp/'scanner';subprocess.run(['c++','-O3','-std=c++17',str(ROOT/'scripts/raw_scan.cpp'),'-o',str(binary)],check=True)
  for row in catalog:
   if row['name'] not in selected:continue
   old=min(float(r['margin']) for r in structures if r['name']==row['name'])
   options=set()
   for end in ('F','R'):
    q=row[end]
    for pos in range(len(q)-3):
     for base in 'ACGT':
      if base==q[pos]:continue
      seq=q[:pos]+base+q[pos+1:]
      if not terminal(seq):continue
      f,r=(seq,row['R']) if end=='F' else (row['F'],seq)
      ok,sc=thermal(f,r)
      if not ok:continue
      margin=min(min(x['F_tm_min'],x['R_tm_min'])-x['structure_max'] for x in sc)
      if min(min(x['F_tm_min'],x['R_tm_min']) for x in sc)<45 or max(max(x['F_tm_max'],x['R_tm_max']) for x in sc)>75 or margin<5:continue
      options.add((f,r,margin))
   options=sorted(options);d=tmp/row['name'];d.mkdir();panel=ROOT/'panels'/row['target_taxon']/'references'/row['gene']
   for role in ('target','background'):shutil.copyfile(panel/(role+'.fa'),d/(role+'.fasta'))
   (d/'independent_primers.tsv').write_text(''.join(f'p{i}\t{f}\t{r}\n' for i,(f,r,m) in enumerate(options)))
   if not options:
    summary.append(dict(name=row['name'],baseline_margin=old,thermally_eligible=0,passes_all=0));continue
   subprocess.run([str(binary),str(d)],stdout=subprocess.DEVNULL,stderr=subprocess.DEVNULL,check=True)
   count={f'p{i}':Counter() for i in range(len(options))};bgmeta={x['accession']:x for x in read(panel/'background.csv')};near=sum(x.get('group')=='other Phyllodocida' for x in bgmeta.values())
   with (d/'independent_hits.tsv').open() as h:
    for x in csv.DictReader(h,delimiter='\t'):
     c=count[x['name']]
     if x['role']=='target':c['covered']+=int(x['le2last3']);c['exact']+=int(x['exact'])
     else:
      model=x['model'];c[model]+=int(x['le2last3']);c[model+'loose']+=int(x['le3' if model=='WC' else 'le2'])
      if bgmeta[x['accession']].get('group')=='other Phyllodocida':c['near'+model]+=int(x['le2last3'])
   passed=[]
   for i,(f,r,margin) in enumerate(options):
    c=count[f'p{i}'];n=len(bgmeta)
    if row['target_taxon']=='Aphroditiformia':bg=max(c['WC'],c['GT'])<=10 and max(c['WCloose'],c['GTloose'])<=50
    else:bg=c['WC']<=int(n*.01) and c['GT']<=int(n*.03) and c['WCloose']<=int(n*.05) and c['GTloose']<=int(n*.10) and max(c['nearWC'],c['nearGT'])<=max(2,int(near*.03))
    good=bg and c['covered']>=int(row['covered']) and c['exact']>=int(row['exact'])
    item=dict(name=row['name'],F=f,R=r,margin=margin,covered=c['covered'],exact=c['exact'],WC=c['WC'],GT=c['GT'],WCloose=c['WCloose'],GTloose=c['GTloose'],nearWC=c['nearWC'],nearGT=c['nearGT'],background_pass=int(bg),target_no_loss=int(c['covered']>=int(row['covered']) and c['exact']>=int(row['exact'])),passes_all=int(good))
    allc.append(item)
    if good:passed.append(item)
   summary.append(dict(name=row['name'],baseline_margin=old,thermally_eligible=len(options),passes_all=len(passed)))
   print(row['name'],'thermal candidates',len(options),'pass',len(passed),flush=True)
 for name,data in [('low_margin_search_summary.csv',summary),('low_margin_search_candidates.csv',allc)]:
  if data:
   with (out/name).open('w') as h:w=csv.DictWriter(h,fieldnames=list(data[0]),lineterminator='\n');w.writeheader();w.writerows(data)
 print('complete',len(summary))
if __name__=='__main__':main()
