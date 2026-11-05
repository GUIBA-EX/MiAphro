from pathlib import Path
import csv,math
root=Path(__file__).resolve().parents[1];out=[]
for gene in ['16S','28S']:
 for band in 'SML':
  rows=list(csv.DictReader((root/'SML'/gene/band/'screen.csv').open()));by={r['id']:r for r in rows};n=int(rows[0]['background_species'])
  for model in ['WC','WC_and_GT']:
   for fraction in [.003,.01,.02,.05]:
    cap=math.floor(n*fraction);pool=[]
    for r in rows:
     if not r['id'].endswith('D'):continue
     a=by[r['id'][:-1]]
     if all(max(int(x[k]) for k in (['WC'] if model=='WC' else ['WC','GT']))<=cap for x in [a,r]):pool.append((a,r))
    if not pool:continue
    a,r=max(pool,key=lambda v:(int(v[1]['covered']),int(v[1]['exact']),int(v[0]['covered'])))
    for x,label in [(a,'ordinary'),(r,'D')]:out.append(dict(gene=gene,band=band,model=model,strict_cap=cap,cap_fraction=fraction,variant=label,**x))
with (root/'SML/tradeoffs.csv').open('w') as h:w=csv.DictWriter(h,out[0]);w.writeheader();w.writerows(out)
print('Tradeoff scenarios written; loose/per-group gates are not applied to these exploratory alternatives.')
