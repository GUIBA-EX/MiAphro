from pathlib import Path
from collections import Counter,defaultdict
import csv,subprocess,shutil,sys,time
from Bio import SeqIO
from Bio.Seq import Seq
from chemistry_round2 import I,thermal,terminal,deg,oligo_thermal
def eligible(f,r):
 ft,rt=oligo_thermal(f),oligo_thermal(r)
 if any(min(a[0],b[0])<45 or max(a[2],b[2])>=55 or abs((a[0]+a[1]-b[0]-b[1])/2)>5 for a,b in zip(ft,rt)):return False
 return thermal(f,r)[0]
root=Path(__file__).resolve().parents[1]
gene=sys.argv[1];source=root/gene
band=sys.argv[2] if len(sys.argv)>2 else None
limits={'S':(100,200),'M':(201,300),'L':(301,480)}
low,high=limits[band] if band else (135,330)
d=root/'Expanded3'/gene/band if band else source
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
import pickle
cache=source/'expanded3_sites.pkl'
if cache.exists():sites=pickle.loads(cache.read_bytes())
else:
 for rev in (0,1):
  for a in range(0,len(cols)-18):
   for le in range(18,31):
    if a+le>len(cols):continue
    seg=[s[cols[a]:cols[a+le-1]+1].replace('-','') for s in al.values()]
    seg=[(str(Seq(w).reverse_complement()) if rev else w) if len(w)==le and set(w)<=set('ACGT') else '' for w in seg];freq=Counter(w for w in seg if w)
    if sum(freq.values())<max(4,N*.06):continue
    consensus=''.join(Counter(w[j] for w in seg if w).most_common(1)[0][0] for j in range(le))
    seeds=[w for w,n in freq.most_common(16)]
    if consensus not in seeds:seeds.append(consensus)
    masks=[{base:sum(1<<i for i,w in enumerate(seg) if w and w[j]==base) for base in 'ACGT'} for j in range(le)]
    valid=sum(1<<i for i,w in enumerate(seg) if w)
    base_counts=[Counter({b:sum(n for w,n in freq.items() if w[j]==b) for b in 'ACGT'}) for j in range(le)]
    from functools import lru_cache
    @lru_cache(None)
    def bits(v):
     z=o=t=valid
     for j,c in enumerate(v):
      good=0
      for base in I[c]:good|=masks[j][base]
      bad=valid^good
      if j>=le-3:z&=good;o&=good;t&=good
      else:t=(t&good)|(o&bad);o=(o&good)|(z&bad);z&=good
     return t,z
    def quality(v):
     t,z=bits(v);return t.bit_count(),z.bit_count(),v
    def mixed(vs,n):
     return list(dict.fromkeys(sorted(vs,key=quality,reverse=True)[:n]+sorted(vs,key=lambda v:(quality(v)[1],quality(v)[0],v),reverse=True)[:n]))
    variants=set()
    for q in seeds:
     if not terminal(q):continue
     beam=[q];variants.add(q)
     for depth in range(3):
      expanded=set()
      for v in beam:
       for j in range(le//2,le):
        for b in 'ACGT':
         if b in I[v[j]] or base_counts[j][b]<2:continue
         ds=list(v);ds[j]=iup[frozenset(I[v[j]]+b)];cand=''.join(ds)
         if deg(cand)<=8 and terminal(cand):expanded.add(cand)
      if not expanded:break
      beam=mixed(expanded,2);variants.update(beam)
    if variants:
     keep=[]
     for dg in [False,True]:keep+=mixed([v for v in variants if (deg(v)>1)==dg],3)
     sites[rev][a,le]=[(v,*bits(v)) for v in keep]
  print('Sites',rev,len(sites[rev]),flush=True)
 cache.write_bytes(pickle.dumps(sites))
# Stage 1: full coordinate sweep using coverage- and exact-oriented site representatives.
# Stage 2: evaluate all retained oligo combinations in the shortlisted windows.
import heapq
Hcov=[];Hex=[];grid={};candidate_count=0
F={key:(max(v,key=lambda t:(t[1].bit_count(),t[2].bit_count())),max(v,key=lambda t:(t[2].bit_count(),t[1].bit_count()))) for key,v in sites[0].items()}
R={key:(max(v,key=lambda t:(t[1].bit_count(),t[2].bit_count())),max(v,key=lambda t:(t[2].bit_count(),t[1].bit_count()))) for key,v in sites[1].items()}
def retain(heap,item,limit=6000):
 if len(heap)<limit:heapq.heappush(heap,item)
 elif item>heap[0]:heapq.heapreplace(heap,item)
for (a,fl),(fc,fe) in F.items():
 for (b,rl),(rc,re) in R.items():
  span=b+rl-a
  if not low<=span<=high:continue
  cov=(fc[1]&rc[1]).bit_count();ex=(fe[2]&re[2]).bit_count();win=(a+1,a+fl,b+1,b+rl);candidate_count+=1
  retain(Hcov,(cov,ex,win));retain(Hex,(ex,cov,win))
  key=(a//12,b//12)
  if key not in grid or (cov,ex,win)>grid[key]:grid[key]=(cov,ex,win)
wins={t[2] for t in Hcov+Hex+list(grid.values())}
# Always re-evaluate earlier screened windows that are representable in this reference.
old=root/'Round2'/gene/band/'screen.csv'
for row in csv.DictReader(old.open()):
 if row['reference']==ref:wins.add(tuple(int(row[k]) for k in ['F_start','F_end','R_start','R_end']))
pools={}
for win in sorted(wins):
 a,ae,b,be=win;fs=sites[0].get((a-1,ae-a+1),[]);rs=sites[1].get((b-1,be-b+1),[]);ordinary=[];degen=[]
 for f,fh,fe in fs:
  for r,rh,re in rs:
   z=((fh&rh).bit_count(),(fe&re).bit_count(),f,r)
   (degen if deg(f)*deg(r)>1 else ordinary).append(z)
 if ordinary and degen:
  def choices(v):return list(dict.fromkeys(sorted(v,reverse=True)[:6]+sorted(v,key=lambda t:(t[1],t[0],t[2],t[3]),reverse=True)[:6]))
  pools[win]=[choices(ordinary),choices(degen)]
print(gene,band,'coordinate windows',candidate_count,'full combinations windows',len(pools),flush=True)
covrank=sorted(pools.items(),key=lambda t:(max(x[0] for x in t[1][1]),max(x[1] for x in t[1][1])),reverse=True)
exrank=sorted(pools.items(),key=lambda t:(max(x[1] for x in t[1][1]),max(x[0] for x in t[1][1])),reverse=True)
groups=defaultdict(list)
for item in covrank:groups[item[0][0]//12,item[0][2]//12].append(item)
diverse=[v[i] for i in range(4) for v in groups.values() if len(v)>i]
ranked=[]
for i in range(max(len(covrank),len(exrank),len(diverse))):
 for seq in (covrank,exrank,diverse):
  if i<len(seq):ranked.append(seq[i])
sel={}
for win,(os,gs) in ranked:
 if win in sel:continue
 if len(sel)>=3000:break
 chosen=[]
 for choices in [os,gs]:
  good=next((x for x in choices if eligible(x[2],x[3])),None)
  if good:chosen.append(good)
 if len(chosen)==2:sel[win]=chosen
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
