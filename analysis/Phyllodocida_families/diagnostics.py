from pathlib import Path
import csv
R=Path(__file__).resolve().parent
rows=[]
for p in sorted(R.glob('*/Round2/*/*/screen.csv')):
 raw=list(csv.DictReader(p.open()))
 control=p.parents[3]/'Controls'/p.parts[-3]/p.parts[-2]/'screen.csv'
 if control.exists():raw+=list(csv.DictReader(control.open()))
 lookup={x['id']:x for x in raw};family=p.parts[-5];gene=p.parts[-3];band=p.parts[-2]
 for scenario,predicate in [('current',lambda x:x['pass']=='True'),('strict_WC_2pct_GT_5pct',lambda x:int(x['WC'])<=int(int(x['background_species'])*.02) and int(x['GT'])<=int(int(x['background_species'])*.05)),('unrestricted',lambda x:True)]:
  candidates=[x for x in raw if x['id'].endswith('D') and predicate(x) and predicate(lookup[x['id'][:-1]])]
  if not candidates:
   rows.append(dict(family=family,gene=gene,band=band,scenario=scenario,id='',covered='',exact='',target_species='',WC='',GT='',background_species=''));continue
  x=max(candidates,key=lambda x:(int(x['covered']),int(lookup[x['id'][:-1]]['covered']),int(x['exact'])+int(lookup[x['id'][:-1]]['exact']),-max(int(x['GT']),int(lookup[x['id'][:-1]]['GT']))))
  rows.append(dict(family=family,gene=gene,band=band,scenario=scenario,**{k:x[k] for k in ['id','covered','exact','target_species','WC','GT','background_species']}))
with (R/'screening_diagnostics.csv').open('w') as h:
 w=csv.DictWriter(h,rows[0]);w.writeheader();w.writerows(rows)
print('Diagnostics',len(rows))
