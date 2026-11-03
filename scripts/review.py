"""Review the published panels and screen candidate consensus corrections."""
import csv
from collections import Counter
from itertools import combinations
from Bio import SeqIO
from Bio.Seq import Seq
from chemistry import I,deg,thermal,terminal
from independent import flags
from template_matching import ROOT,render

def rows(p):return list(csv.DictReader(p.open()))
def seqs(p):return {r.id.removeprefix('_R_'):str(r.seq).upper() for r in SeqIO.parse(p,'fasta')}
def write(name,rows):
 with (ROOT/'audit'/name).open('w') as f:
  w=csv.DictWriter(f,rows[0],lineterminator='\n');w.writeheader();w.writerows(rows)
support=[];trials=[];audit=[]
for base in [ROOT/'baseline',ROOT/'SML']:
 pairs=rows(base/'results.csv');oligos=seqs(base/'primers.fasta')
 for p in pairs:
  name=p['name'];gene=name.split('-')[1][:3];d=ROOT/gene
  al=seqs(d/'aligned.fa');targets=seqs(d/'target.fa');bg=seqs(d/'background.fa')
  cols=[i for i,b in enumerate(al[p['reference']]) if b!='-'];edits=[]
  for side in 'FR':
   q=p[side];a,b=[int(p[side+'_'+k]) for k in ['start','end']]
   assert q==oligos[name+'_'+side] and len(q)==b-a+1 and deg(q)<=8 and terminal(q) and set(q[:len(q)//2])<=set('ACGT')
   words=[]
   for s in al.values():
    site=s[cols[a-1]:cols[b-1]+1].replace('-','')
    word=''.join(s[i] for i in cols[a-1:b])
    if word==site and set(word)<=set('ACGT'):words.append(str(Seq(word).reverse_complement()) if side=='R' else word)
   for j,basecall in enumerate(q):
    counts=Counter(w[j] for w in words);major,count=counts.most_common(1)[0]
    support.append(dict(name=name,side=side,position=j+1,base=basecall,n=len(words),**{b:counts[b] for b in 'ACGT'}))
    if major not in I[basecall] and count>max(counts[b] for b in I[basecall]):edits.append((side,j,major))
  twin=next(r for r in pairs if r['name']==(name[:-1] if name.endswith('D') else name+'D'))
  assert all(p[k]==twin[k] for k in ('reference','F_start','F_end','R_start','R_end'))
  def evaluate(f,r):
   raw={k:flags(s,f,r,True,False) for k,s in targets.items()}
   h={k for k,v in raw.items() if v[1]};e={k for k,v in raw.items() if v[0]}
   back=[set() for _ in range(4)]
   for k,s in bg.items():
    for m in (0,1):
     _,v,w=flags(s,f,r,False,m)
     if v:back[m].add(k)
     if w:back[m+2].add(k)
   return h,e,back
  h,e,back=evaluate(p['F'],p['R'])
  assert list(map(len,[h,e,*back]))==[int(p[k]) for k in ['covered','exact','WC','GT','WCloose','GTloose']]
  valid,tm=thermal(p['F'],p['R']);assert valid and all(min(x['F_tm_min'],x['R_tm_min'])>=45 for x in tm)
  audit.append(dict(name=name,degF=deg(p['F']),degR=deg(p['R']),consensus_flags=len(edits),checks=len(targets)+2*len(bg)))
  for n in range(1,len(edits)+1):
   for combo in combinations(edits,n):
    seq={s:list(p[s]) for s in 'FR'}
    for s,j,b in combo:seq[s][j]=b
    f,r=[''.join(seq[s]) for s in 'FR'];hh,ee,bb=evaluate(f,r)
    trials.append(dict(name=name,edits=' '.join(f'{s}{j+1}{p[s][j]}>{b}' for s,j,b in combo),F=f,R=r,covered=len(hh),exact=len(ee),lost=len(h-hh),new_background=sum(len(v-u) for u,v in zip(back,bb)),thermal=thermal(f,r)[0],terminal=terminal(f) and terminal(r)))
  if base.name=='SML':
   detail=rows(d/'target.csv')
   for x in detail:x.update(covered=str(int(x['accession'] in h)),exact=str(int(x['accession'] in e)))
   render(gene,p,detail,al,cols)
   if name in ['MiAphro-16S-SD','MiAphro-28S-MD']:render(gene,p,detail,al,cols,True)
write('pair_review.csv',audit);write('position_support.csv',support)
if trials:write('consensus_trials.csv',trials)
