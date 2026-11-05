from pathlib import Path
import csv
R=Path(__file__).resolve().parent
read=lambda p:list(csv.DictReader(p.open()))
for fam,gene in [('Syllidae','16S'),('Nereididae','28S')]:
 r=R/fam;old={x['name']:x for x in read(r/'Round2/results.csv')}
 for band in 'SML':
  if not (r/'Expanded3'/gene/band/'complete.flag').exists():continue
  pool=[]
  for src in ['Round2','Controls','Expanded3','Refine3']:
   p=r/src/gene/band/'screen.csv'
   if not p.exists():continue
   rs=[x for x in read(p) if x['pass']=='True' and int(x['covered'])>0]
   for a in rs:
    if a['id'].endswith('D'):continue
    for b in rs:
     if b['id'].endswith('D') and all(a[k]==b[k] for k in ['reference','F_start','F_end','R_start','R_end']):pool.append((a,b,src))
  a0,b0=old[f'{fam}-{gene}-{band}'],old[f'{fam}-{gene}-{band}D'];rank=lambda v:(int(v[1]['covered']),int(v[1]['exact']),int(v[0]['covered']),int(v[0]['exact']),-max(int(v[0]['GT']),int(v[1]['GT'])))
  bal=[v for v in pool if all(int(v[i][k])>=int([a0,b0][i][k]) for i in [0,1] for k in ['covered','exact'])]
  for mode,v in [('balanced',max(bal,key=rank)),('coverage',max(pool,key=rank)),('exact',max(pool,key=lambda v:(int(v[1]['exact']),*rank(v))))]:
   print(fam,gene,band,mode,v[2],[(x['covered'],x['exact'],x['WC'],x['GT']) for x in v[:2]])
