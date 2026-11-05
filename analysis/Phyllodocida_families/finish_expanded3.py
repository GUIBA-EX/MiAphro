from pathlib import Path
import csv,sys,json,subprocess,io,re,time
from collections import Counter
from Bio import SeqIO
R=Path(__file__).resolve().parent;O=R/'Optimization3';O.mkdir(exist_ok=True)
sys.path.insert(0,str(R/'Syllidae/scripts'))
from independent import flags
from chemistry_round2 import deg,terminal,thermal
read=lambda p:list(csv.DictReader(p.open()))
def write(p,rows):
 with p.open('w') as h:
  w=csv.DictWriter(h,rows[0],lineterminator='\n');w.writeheader();w.writerows(rows)
results=[];details=[];validations=[];summary=[];memo={};near_memo={};seqs={k:[] for k in ['balanced','coverage','exact','strict_profile']}
for fam,gene in [('Nereididae','28S'),('Syllidae','16S')]:
 root=R/fam;metadata={x['accession']:x for x in read(root/gene/'target.csv')};bg={x['accession']:x for x in read(root/gene/'background.csv')};libs={role:{x.id:str(x.seq).upper() for x in SeqIO.parse(root/gene/(role+'.fa'),'fasta')} for role in ['target','background']};old={x['name']:x for x in read(root/'Round2/results.csv')}
 for band in 'SML':
  while not (root/'Refine3'/gene/band/'complete.flag').exists():time.sleep(5)
  pool=[];allpool=[]
  for source in ['Round2','Controls','Expanded3','Refine3']:
   path=root/source/gene/band;rs=read(path/'screen.csv');look={x['id']:x for x in rs}
   for a in rs:
    if a['id'].endswith('D'):continue
    partners=([b for b in rs if b['id'].endswith('D') and all(a[k]==b[k] for k in ['reference','F_start','F_end','R_start','R_end'])] if source=='Refine3' else [look[a['id']+'D']])
    for b in partners:
     assert all(a[k]==b[k] for k in ['reference','F_start','F_end','R_start','R_end'])
     if all(int(x['covered'])>0 for x in [a,b]):
      allpool.append((a,b,path))
      if all(x['pass']=='True' for x in [a,b]):pool.append((a,b,path))
  a0,b0=old[f'{fam}-{gene}-{band}'],old[f'{fam}-{gene}-{band}D']
  def rank(v):
   a,b,p=v;return (int(b['covered']),int(b['exact']),int(a['covered']),int(a['exact']),-max(int(a['GT']),int(b['GT'])))
  balanced=[v for v in pool if all(int(v[i][k])>=int([a0,b0][i][k]) for i in [0,1] for k in ['covered','exact'])]
  assert balanced,'Baseline unexpectedly missing'
  chosen={'balanced':max(balanced,key=rank),'coverage':max(pool,key=rank),'exact':max(pool,key=lambda v:(int(v[1]['exact']),*rank(v)))}
  # Additional conditional profile: retain strict global and close-relative gates,
  # but report loose-model hits diagnostically rather than use them as a veto.
  strict_candidates=[v for v in allpool if all(int(v[i][k])>=int([a0,b0][i][k]) for i in [0,1] for k in ['covered','exact']) and all(int(x['WC'])<=int(len(bg)*.01) and int(x['GT'])<=int(len(bg)*.03) for x in v[:2])]
  near_limit=max(2,int(sum(z['group']=='other Phyllodocida' for z in bg.values())*.03))
  for v in sorted(strict_candidates,key=rank,reverse=True):
   near_ok=True
   for x in v[:2]:
    path=v[2]
    near_key=(str(path),x['id'])
    if near_key not in near_memo:
     with (path/'independent_hits.tsv').open() as h:header=h.readline()
     matched=subprocess.run(['rg','--no-filename','^'+re.escape(x['id'])+r'\t',str(path/'independent_hits.tsv')],capture_output=True,text=True)
     assert matched.returncode in [0,1]
     counts=Counter()
     for z in csv.DictReader(io.StringIO(header+matched.stdout),delimiter='\t'):
      if z['role']=='background' and bg[z['accession']]['group']=='other Phyllodocida':counts[z['model']]+=int(z['le2last3'])
     near_memo[near_key]=max(counts.values(),default=0)
    if near_memo[near_key]>near_limit:near_ok=False;break
   if near_ok:chosen['strict_profile']=v;break
  assert 'strict_profile' in chosen
  summary.append(dict(family=fam,gene=gene,band=band,passing_pairs=len(pool),noninferior_pairs=len(balanced)))
  for objective,(a,b,path) in chosen.items():
   for x in [a,b]:
    suffix='D' if x['id'].endswith('D') else '';name=f'{fam}-{gene}-{band}{suffix}';key=(fam,gene,x['F'],x['R'])
    for end in 'FR':
     q=x[end];assert len(q)==int(x[end+'_end'])-int(x[end+'_start'])+1
     assert deg(q)<=8 and terminal(q) and set(q[:len(q)//2])<=set('ACGT')
     if not suffix:assert deg(q)==1
    ok,tr=thermal(x['F'],x['R']);assert ok and min(min(t['F_tm_min'],t['R_tm_min']) for t in tr)>=45
    if key not in memo:
     with (path/'independent_hits.tsv').open() as h:header=h.readline()
     matched=subprocess.run(['rg','--no-filename','^'+re.escape(x['id'])+r'\t',str(path/'independent_hits.tsv')],capture_output=True,text=True)
     assert matched.returncode in [0,1]
     baseline={(z['role'],z['model'],z['accession']):z for z in csv.DictReader(io.StringIO(header+matched.stdout),delimiter='\t')}
     ct=Counter();det=[];groups=Counter()
     for role,seq in libs.items():
      for acc,s in seq.items():
       for gt in ([False] if role=='target' else [False,True]):
        model='GT' if gt else 'WC';e,c,l=flags(s,x['F'],x['R'],role=='target',gt);z=baseline.get((role,model,acc),{})
        assert [e,c,l]==[z.get(k)=='1' for k in ['exact','le2last3','le3' if role=='background' and not gt else 'le2']],(name,acc)
        if role=='target':ct['covered']+=c;ct['exact']+=e;det.append(dict(**metadata[acc],covered=int(c),exact=int(e)))
        else:ct[model]+=c;ct[model+'loose']+=l;groups[bg[acc]['group'],model]+=c
     memo[key]=(ct,det,groups)
    ct,det,groups=memo[key];assert all(int(x[k])==ct[k] for k in ['covered','exact','WC','GT','WCloose','GTloose'])
    results.append(dict(objective=objective,name=name,source_pool=path.parents[1].name,**x));details.extend(dict(objective=objective,name=name,**z) for z in det)
    validations.append(dict(objective=objective,name=name,independent='PASS',shared_coordinates=True,deg_F=deg(x['F']),deg_R=deg(x['R']),thermal='PASS'))
    seqs[objective].extend(f'>{name}_{end}\n{x[end]}\n' for end in 'FR')
  print('Verified',fam,gene,band,flush=True)
write(O/'results.csv',results);write(O/'species.csv',details);write(O/'search_summary.csv',summary);(O/'validation.json').write_text(json.dumps(validations,indent=2))
for objective,seq in seqs.items():(O/('primers.fasta' if objective=='balanced' else objective+'_primers.fasta')).write_text(''.join(seq))
print('PASS',len(memo),'unique pairs',len(results),'objective records',flush=True)
