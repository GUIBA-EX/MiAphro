from pathlib import Path
from collections import Counter,defaultdict
import csv,subprocess,shutil,sys,time
from Bio import SeqIO
from Bio.Seq import Seq
from chemistry import I,thermal,terminal,deg
root=Path(__file__).resolve().parents[1]
gene=sys.argv[1];source=root/gene
band=sys.argv[2] if len(sys.argv)>2 else None
limits={'S':(100,200),'M':(201,300),'L':(301,480)}
low,high=limits[band] if band else (135,330)
d=root/'Expanded'/gene/band if band else source
d.mkdir(parents=True,exist_ok=True)
if (d/'complete.flag').exists():sys.exit(0)
try:(d/'running.lock').mkdir()
except FileExistsError:
 while not (d/'complete.flag').exists():time.sleep(5)
 sys.exit(0)
def get(f):return {r.id.removeprefix('_R_'):str(r.seq).upper() for r in SeqIO.parse(f,'fasta')}
if not band:
 with open(d/'aligned.fa','w') as h:subprocess.run(['mafft','--quiet','--thread','2','--adjustdirection','--retree','1','--maxiterate','0',str(d/'target.fa')],stdout=h,check=True)
al=get(source/'aligned.fa');ref=max(al,key=lambda k:len(al[k].replace('-','')));cols=[i for i,b in enumerate(al[ref]) if b!='-'];N=len(al);sites=[{},{}];iup={frozenset(v):k for k,v in I.items() if k != 'X'}
for rev in (0,1):
 for a in range(0,len(cols)-24,2 if band else 3):
  for le in (18,20,22,24):
   seg=[s[cols[a]:cols[a+le-1]+1].replace('-','') for s in al.values()]
   seg=[(str(Seq(w).reverse_complement()) if rev else w) if len(w)==le and set(w)<=set('ACGT') else '' for w in seg];freq=Counter(w for w in seg if w)
   if sum(freq.values())<max(6,N*.12):continue
   consensus=''.join(Counter(w[j] for w in seg if w).most_common(1)[0][0] for j in range(le))
   freq_seeds=list(freq.most_common(6))
   if consensus not in dict(freq_seeds):freq_seeds.append((consensus,0))
   opts=[]
   for q,_ in freq_seeds:
    if not terminal(q):continue
    ds=list(q);variants=[q]
    for _,j,b in sorted([(sum(v for w,v in freq.items() if w[j]==b),j,b) for j in range(le//2,le) for b in 'ACGT' if q[j]!=b],reverse=True):
     if _<2:break
     prev=ds[j];ds[j]=iup[frozenset(I[prev]+b)]
     cand=''.join(ds)
     if deg(cand)>8 or not terminal(cand):ds[j]=prev;continue
     variants.append(cand)
    for v in variants:
     byword={w:[j for j,b in enumerate(w) if b not in I[v[j]]] for w in freq}
     hit=sum(1<<i for i,w in enumerate(seg) if w and len(byword[w])<=2 and all(j<le-3 for j in byword[w]))
     ex=sum(1<<i for i,w in enumerate(seg) if w and not byword[w])
     opts.append((v,hit,ex))
   if opts:
    unique={v[0]:v for v in opts}
    sites[rev][a,le]=sum([sorted([v for v in unique.values() if (deg(v[0])>1)==dg],key=lambda v:(v[1].bit_count(),v[2].bit_count()),reverse=True)[:4] for dg in [False,True]],[])
pools=defaultdict(list)
for (a,fl),fs in sites[0].items():
 for (b,rl),rs in sites[1].items():
  span=b+rl-a
  if not low<=span<=high:continue
  ordinary=[];degen=[]
  for f,fh,fe in fs:
   for r,rh,re in rs:
    z=((fh&rh).bit_count(),(fe&re).bit_count(),f,r)
    (degen if deg(f)*deg(r)>1 else ordinary).append(z)
  if ordinary and degen:
   window=(a+1,a+fl,b+1,b+rl)
   pools[window]=[max(ordinary),max(degen)]
ranked=sorted(pools.items(),key=lambda t:(sum(x[0] for x in t[1]),sum(x[1] for x in t[1])),reverse=True)
print(gene,band,'candidate windows',len(ranked),flush=True)
if band:
 # Interleave coordinate neighborhoods so conserved but nonspecific sites cannot monopolize the shortlist.
 groups=defaultdict(list)
 for item in ranked:
  win=item[0];groups[win[0]//25,win[2]//25].append(item)
 diverse=[v[i] for i in range(min(8,max(map(len,groups.values()),default=0))) for v in groups.values() if len(v)>i]
 ranked=ranked[:300]+diverse+ranked
sel={}
for win,(o,g) in ranked:
 if win in sel:continue
 if len(sel)>=600:break
 if not all(thermal(x[2],x[3])[0] and min(min(v['F_tm_min'],v['R_tm_min']) for v in thermal(x[2],x[3])[1])>=45 for x in [o,g]):continue
 sel[win]=[o,g]
print(gene,band,'selected windows',len(sel),flush=True)
while not (source/'background.csv').exists():time.sleep(5)
for k in ['target','background']:shutil.copyfile(source/(k+'.fa'),d/(k+'.fasta'))
pairs={}
for i,(win,tw) in enumerate(sel.items()):
 for j,x in enumerate(tw):pairs[f'p{i}'+('D' if j else '')]=(win,*x[2:])
(d/'independent_primers.tsv').write_text(''.join(f'{k}\t{f}\t{r}\n' for k,(_,f,r) in pairs.items()))
subprocess.run([str(root/'scripts'/'scan'),str(d)],stderr=subprocess.DEVNULL,check=True)
bg_metadata={r['accession']:r for r in csv.DictReader((source/'background.csv').open())}
near_n=sum(r['group']=='other Phyllodocida' for r in bg_metadata.values())
stats=defaultdict(Counter)
for x in csv.DictReader(open(d/'independent_hits.tsv'),delimiter='\t'):
 s=stats[x['name']]
 if x['role']=='target':
  s['covered']+=int(x['le2last3']);s['exact']+=int(x['exact'])
 else:
  m=x['model']
  if bg_metadata[x['accession']]['group']=='other Phyllodocida':s['near'+m]+=int(x['le2last3'])
  s[m]+=int(x['le2last3']);s[m+'loose']+=int(x['le2'] if m=='GT' else x['le3'])
background_n=len(get(source/'background.fa'))
def passed(k):
 s=stats[k];n=background_n;return s['WC']<=int(n*.01) and s['GT']<=int(n*.03) and s['WCloose']<=int(n*.05) and s['GTloose']<=int(n*.10) and max(s['nearWC'],s['nearGT'])<=max(2,int(near_n*.03))
order=sorted([k for k in pairs if not k.endswith('D')],key=lambda k:(passed(k) and passed(k+'D'),min(stats[k]['covered'],stats[k+'D']['covered']),sum(stats[n]['exact'] for n in [k,k+'D'])),reverse=True)
with open(d/'screen.csv','w') as h:
 w=csv.writer(h);w.writerow(['id','reference','F_start','F_end','R_start','R_end','F','R','target_species','covered','exact','background_species','WC','GT','WCloose','GTloose','pass'])
 for k in order:
  for n in (k,k+'D'):
   win,f,r=pairs[n];s=stats[n];w.writerow([n,ref,*win,f,r,N,s['covered'],s['exact'],background_n,*[s[key] for key in ['WC','GT','WCloose','GTloose']],passed(n)])

(d/'complete.flag').write_text('complete\n')
(d/'running.lock').rmdir()
