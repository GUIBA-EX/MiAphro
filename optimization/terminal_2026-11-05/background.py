from pathlib import Path
import csv,json,shutil,subprocess,sys,math
from collections import Counter
from Bio import SeqIO
R=Path(__file__).resolve().parents[2];O=Path(__file__).resolve().parent
sys.path.insert(0,str(R/'scripts'));from independent import flags
read=lambda p:list(csv.DictReader(p.open()));rows=read(O/'baseline_catalog.csv');summary=[]
for d in sorted(O.glob('MiPhyllo-*')):
 if not (d/'joint_pairs.json').exists():continue
 x=next(z for z in rows if z['name']==d.name);allc=json.loads((d/'joint_pairs.json').read_text());cs=[z for z in allc if z['thermal_pass'] or z['baseline']];p=R/'panels'/x['target_taxon']/'references'/x['gene'];scan=d/'background_scan';scan.mkdir(exist_ok=True)
 for role in ['target','background']:shutil.copyfile(p/(role+'.fa'),scan/(role+'.fasta'))
 (scan/'independent_primers.tsv').write_text(''.join(f"p{i}\t{z['F']}\t{z['R']}\n" for i,z in enumerate(cs)))
 subprocess.run([str(R/'scripts/scan'),str(scan)],check=True,stderr=subprocess.DEVNULL)
 meta={z['accession']:z for z in read(p/'background.csv')};near=sum(z.get('group')=='other Phyllodocida' for z in meta.values());stats={f'p{i}':Counter() for i in range(len(cs))}
 for h in csv.DictReader((scan/'independent_hits.tsv').open(),delimiter='\t'):
  s=stats[h['name']]
  if h['role']=='target':s['covered']+=int(h['le2last3']);s['exact']+=int(h['exact'])
  else:
   m=h['model'];s[m]+=int(h['le2last3']);s[m+'loose']+=int(h['le3' if m=='WC' else 'le2'])
   if meta[h['accession']].get('group')=='other Phyllodocida':s['near'+m]+=int(h['le2last3'])
 evaluated=[]
 for i,z in enumerate(cs):
  s=stats[f'p{i}'];assert (s['covered'],s['exact'])==(z['covered'],z['exact'])
  if x['target_taxon']=='Aphroditiformia':ok=max(s['WC'],s['GT'])<=10 and max(s['WCloose'],s['GTloose'])<=50
  else:ok=s['WC']<=int(len(meta)*.01) and s['GT']<=int(len(meta)*.03) and s['WCloose']<=int(len(meta)*.05) and s['GTloose']<=int(len(meta)*.1) and max(s['nearWC'],s['nearGT'])<=max(2,int(near*.03))
  evaluated.append(dict(z,**{k:s[k] for k in ['WC','GT','WCloose','GTloose','nearWC','nearGT']},background_pass=ok))
 baseline=next(z for z in evaluated if z['baseline']);assert all(int(x[k])==baseline[k] for k in ['covered','exact','WC','GT','WCloose','GTloose'])
 passed=[z for z in evaluated if z['background_pass'] and z['thermal_pass'] and not z['baseline']]
 chosen=max(passed,key=lambda z:(z['covered'],z['exact'],-z['GT'],-z['WCloose'],-z['GTloose'])) if passed else baseline
 # Independently re-scan chosen pair's complete target and background panels in Python.
 counts=Counter()
 for role in ['target','background']:
  for rec in SeqIO.parse(p/(role+'.fa'),'fasta'):
   for gt in ([False] if role=='target' else [False,True]):
    e,c,l=flags(str(rec.seq).upper(),chosen['F'],chosen['R'],role=='target',gt)
    if role=='target':counts['covered']+=c;counts['exact']+=e
    else:m='GT' if gt else 'WC';counts[m]+=c;counts[m+'loose']+=l
 assert all(counts[k]==chosen[k] for k in ['covered','exact','WC','GT','WCloose','GTloose'])
 (d/'background_results.json').write_text(json.dumps(evaluated,indent=2));(d/'selected.json').write_text(json.dumps(chosen,indent=2))
 summary.append(dict(name=x['name'],status='retain' if chosen['baseline'] else 'replace',before=baseline,after=chosen,target_N=x['target_species'],background_N=x['background_species'],target_candidates=len(allc)-1,thermal_candidates=len(cs)-1,background_passed=len(passed),independent_status='PASS'))
 print(x['name'],summary[-1]['status'],(baseline['covered'],baseline['exact']),'->',(chosen['covered'],chosen['exact']), 'bg',*[chosen[k] for k in ['WC','GT','WCloose','GTloose']],flush=True)
(O/'comparison.json').write_text(json.dumps(summary,indent=2))
