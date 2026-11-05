from pathlib import Path
import sys,csv,pickle,subprocess,shutil,time,concurrent.futures
from collections import Counter,defaultdict
from functools import lru_cache
from Bio import SeqIO
from Bio.Seq import Seq
R=Path(__file__).resolve().parent;sys.path.insert(0,str(R/'Syllidae/scripts'))
from chemistry_round2 import I,deg,thermal,terminal
read=lambda p:list(csv.DictReader(p.open()))
iup={frozenset(v):k for k,v in I.items() if k!='X'}
def options(seg,le):
 freq=Counter(w for w in seg if len(w)==le and set(w)<=set('ACGT'))
 if not freq:return []
 valid=sum(1<<i for i,w in enumerate(seg) if len(w)==le and set(w)<=set('ACGT'))
 masks=[{b:sum(1<<i for i,w in enumerate(seg) if len(w)==le and set(w)<=set('ACGT') and w[j]==b) for b in 'ACGT'} for j in range(le)]
 counts=[Counter({b:sum(n for w,n in freq.items() if w[j]==b) for b in 'ACGT'}) for j in range(le)]
 @lru_cache(None)
 def bits(v):
  z=o=t=valid
  for j,c in enumerate(v):
   good=0
   for b in I[c]:good|=masks[j][b]
   bad=valid^good
   if j>=le-3:z&=good;o&=good;t&=good
   else:t=(t&good)|(o&bad);o=(o&good)|(z&bad);z&=good
  return t,z
 def rank(q):
  c,e=bits(q);return c.bit_count(),e.bit_count(),q
 def mixed(v,n):return list(dict.fromkeys(sorted(v,key=rank,reverse=True)[:n]+sorted(v,key=lambda q:(rank(q)[1],rank(q)[0],q),reverse=True)[:n]))
 seeds=list(dict.fromkeys([w for w,n in freq.most_common(32)]+[''.join(c.most_common(1)[0][0] for c in counts)]))
 # Single-base consensus-guided improvements of ordinary seeds.
 ordinary=set(q for q in seeds if terminal(q))
 for q in seeds[:8]:
  for _ in range(2):
   neighbors={q}
   for j in range(le):
    b=counts[j].most_common(1)[0][0];v=q[:j]+b+q[j+1:]
    if terminal(v):neighbors.add(v)
   best=max(neighbors,key=rank)
   if best==q:break
   q=best;ordinary.add(q)
 variants=set(ordinary)
 for q in mixed(ordinary,12):
  beam=[q]
  for _ in range(3):
   expanded=set()
   for v in beam:
    for j in range(le//2,le):
     for b in 'ACGT':
      if b in I[v[j]] or counts[j][b]<2:continue
      w=v[:j]+iup[frozenset(I[v[j]]+b)]+v[j+1:]
      if deg(w)<=8 and terminal(w):expanded.add(w)
   if not expanded:break
   beam=mixed(expanded,4);variants.update(beam)
 keep=[]
 for dg in [False,True]:keep+=mixed([q for q in variants if (deg(q)>1)==dg],6)
 return [(q,*bits(q)) for q in keep]
def run(job):
 fam,gene=job;root=R/fam;src=root/gene;target=read(src/'target.csv');bg={x['accession']:x for x in read(src/'background.csv')};N=len(target);n=len(bg);near=sum(x['group']=='other Phyllodocida' for x in bg.values())
 for band in 'SML':
  while not (root/'Expanded3'/gene/band/'complete.flag').exists():time.sleep(5)
  out=root/'Refine3'/gene/band;out.mkdir(parents=True,exist_ok=True)
  if (out/'complete.flag').exists():continue
  candidates=[]
  for source in ['Expanded3','Round2','Controls']:
   candidates+=read(root/source/gene/band/'screen.csv')
  rows=[x for x in candidates if x['id'].endswith('D') and x['pass']=='True'];old=[x for x in read(root/'Round2/results.csv') if x['name']==f'{fam}-{gene}-{band}D'];ordered=old[:]
  rankings=[sorted(rows,key=lambda x:(int(x['covered']),int(x['exact'])),reverse=True),sorted(rows,key=lambda x:(int(x['exact']),int(x['covered'])),reverse=True)]
  for i in range(max(map(len,rankings),default=0)):
   for rs in rankings:
    if i<len(rs):ordered.append(rs[i])
  selected={}
  for x in ordered:
   key=tuple(x[k] for k in ['reference','F_start','F_end','R_start','R_end'])
   selected.setdefault(key,x)
   if len(selected)>=24:break
  alignments=[]
  for p in [src/'aligned.fa',root/'Controls'/gene/'coordinate_reference_alignment.fa']:
   alignments.append({x.id.removeprefix('_R_'):str(x.seq).upper() for x in SeqIO.parse(p,'fasta')})
  pairs={};refseq={};index=0
  for key,x in selected.items():
   ref=key[0];al=next(a for a in alignments if ref in a);cols=[i for i,c in enumerate(al[ref]) if c!='-'];opts=[]
   for end in 'FR':
    a,b=int(x[end+'_start']),int(x[end+'_end']);le=b-a+1
    seg=[al[z['accession']][cols[a-1]:cols[b-1]+1].replace('-','') for z in target]
    if end=='R':seg=[str(Seq(w).reverse_complement()) for w in seg]
    opts.append(options(seg,le))
   ordinary=[];degen=[]
   for f,fh,fe in opts[0]:
    for r,rh,re in opts[1]:
     z=((fh&rh).bit_count(),(fe&re).bit_count(),f,r)
     (degen if deg(f)*deg(r)>1 else ordinary).append(z)
   for dg,rs in [(False,ordinary),(True,degen)]:
    chosen=list(dict.fromkeys(sorted(rs,reverse=True)[:20]+sorted(rs,key=lambda z:(z[1],z[0],z[2],z[3]),reverse=True)[:20]));accepted=[]
    for objective in ['coverage','exact']:
     ordered_pairs=sorted(chosen,key=(lambda z:(z[0],z[1],z[2],z[3])) if objective=='coverage' else (lambda z:(z[1],z[0],z[2],z[3])),reverse=True);ct=0
     for z in ordered_pairs:
      ok,tr=thermal(z[2],z[3])
      if ok and min(min(t['F_tm_min'],t['R_tm_min']) for t in tr)>=45:
       if z not in accepted:accepted.append(z)
       ct+=1
       if ct>=5:break
    for z in accepted:
     name=f'r{index}'+('D' if dg else '');index+=1;pairs[name]=dict(reference=ref,**{k:x[k] for k in ['F_start','F_end','R_start','R_end']},F=z[2],R=z[3])
  print('Refine scan',fam,gene,band,len(selected),'windows',len(pairs),'pairs',flush=True)
  for role in ['target','background']:shutil.copyfile(src/(role+'.fa'),out/(role+'.fasta'))
  (out/'independent_primers.tsv').write_text(''.join(f"{k}\t{x['F']}\t{x['R']}\n" for k,x in pairs.items()))
  subprocess.run([str(root/'scripts/scan'),str(out)],stderr=subprocess.DEVNULL,check=True);stats=defaultdict(Counter)
  with (out/'independent_hits.tsv').open() as h:
   for z in csv.DictReader(h,delimiter='\t'):
    s=stats[z['name']]
    if z['role']=='target':s['covered']+=int(z['le2last3']);s['exact']+=int(z['exact'])
    else:
     m=z['model'];s[m]+=int(z['le2last3']);s[m+'loose']+=int(z['le2'] if m=='GT' else z['le3'])
     if bg[z['accession']]['group']=='other Phyllodocida':s['near'+m]+=int(z['le2last3'])
  output=[]
  for name,x in pairs.items():
   s=stats[name];passed=s['WC']<=int(n*.01) and s['GT']<=int(n*.03) and s['WCloose']<=int(n*.05) and s['GTloose']<=int(n*.10) and max(s['nearWC'],s['nearGT'])<=max(2,int(near*.03))
   output.append(dict(id=name,**x,target_species=N,background_species=n,**{k:s[k] for k in ['covered','exact','WC','GT','WCloose','GTloose']},**{'pass':passed}))
  with (out/'screen.csv').open('w') as h:w=csv.DictWriter(h,output[0]);w.writeheader();w.writerows(output)
  (out/'complete.flag').write_text('complete\n');print('Refine complete',fam,gene,band,flush=True)
with concurrent.futures.ThreadPoolExecutor(max_workers=2) as p:list(p.map(run,[('Syllidae','16S'),('Nereididae','28S')]))
