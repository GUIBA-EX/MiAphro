"""Independently verify low-margin candidates and accession-level gains/losses."""
import csv,json,sys
from pathlib import Path
from Bio import SeqIO
ROOT=Path(__file__).resolve().parents[1];sys.path.insert(0,str(ROOT/'scripts'))
from independent import flags

def read(p):
 with p.open() as f:return list(csv.DictReader(f))

def main():
 out=ROOT/'audit/review_20261107';catalog={r['name']:r for r in read(ROOT/'catalog.csv')};res=[]
 for c in read(out/'low_margin_search_candidates.csv'):
  if c['passes_all']!='1':continue
  orig=catalog[c['name']];panel=ROOT/'panels'/orig['target_taxon']/'references'/orig['gene'];sets={}
  for label,pair in [('original',orig),('candidate',c)]:
   for role in ('target','background'):
    for rec in SeqIO.parse(panel/(role+'.fa'),'fasta'):
     for gt in ([False] if role=='target' else [False,True]):
      ex,strict,loose=flags(str(rec.seq).upper(),pair['F'],pair['R'],role=='target',gt)
      for metric,state in [('exact',ex),('strict',strict),('loose',loose)]:
       if state:sets.setdefault((label,role,'GT' if gt else 'WC',metric),set()).add(rec.id)
  def get(a,b,c,d):return sets.get((a,b,c,d),set())
  counts={'covered':len(get('candidate','target','WC','strict')),'exact':len(get('candidate','target','WC','exact')),'WC':len(get('candidate','background','WC','strict')),'GT':len(get('candidate','background','GT','strict')),'WCloose':len(get('candidate','background','WC','loose')),'GTloose':len(get('candidate','background','GT','loose'))}
  if any(counts[k]!=int(c[k]) for k in counts):raise ValueError((c['name'],counts))
  near={r['accession'] for r in read(panel/'background.csv') if r.get('group')=='other Phyllodocida'}
  for model in ('WC','GT'):
   got=len(get('candidate','background',model,'strict') & near)
   if got!=int(c['near'+model]):raise ValueError((c['name'],'near'+model,got,c['near'+model]))
   counts['near'+model]=got
  compare={}
  for metric in ('strict','exact'):
   a=get('original','target','WC',metric);b=get('candidate','target','WC',metric)
   compare[metric]={'lost':sorted(a-b),'gained':sorted(b-a),'retained':len(a&b)}
  res.append({'name':c['name'],'F':c['F'],'R':c['R'],'candidate_counts':counts,'target_accession_comparison':compare,'status':'PASS','same_target_set':all(not x['lost'] and not x['gained'] for x in compare.values())})
 (out/'low_margin_independent_verification.json').write_text(json.dumps(res,indent=2)+'\n');print(json.dumps(res,indent=2))
if __name__=='__main__':main()
