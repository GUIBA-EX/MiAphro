from Bio.Data.IUPACData import ambiguous_dna_values as I
from functools import lru_cache
import re
comp=str.maketrans('ACGT','TGCA')
@lru_cache(None)
def pattern(q,reverse,gt,k):
 sets=[set(I[b]) for b in q]
 for chars in sets:
  if gt:
   if 'G' in chars:chars.add('A')
   if 'T' in chars:chars.add('C')
 if reverse:sets=[set(''.join(c).translate(comp)) for c in sets[::-1]]
 seeds=[]
 for j in range(k+1):
  a=j*len(q)//(k+1);b=(j+1)*len(q)//(k+1)
  seeds.append((a,re.compile('(?='+''.join('['+''.join(sorted(c))+']' for c in sets[a:b])+')')))
 return sets,seeds
def sites(s,q,rev,gt,k):
 sets,seeds=pattern(q,rev,gt,k)
 starts={m.start()-a for a,rg in seeds for m in rg.finditer(s)}
 for i in starts:
  if i<0 or i+len(q)>len(s):continue
  t=s[i:i+len(q)]
  if set(t)-set('ACGT'):continue
  bad=[j for j,c in enumerate(t) if c not in sets[j]]
  if len(bad)<=k:yield i,len(bad),all(j>=3 if rev else j<len(q)-3 for j in bad)
def flags(s,f,r,target,gt):
 k=2 if target or gt else 3;ans=[False]*3
 for strand in (s,s.translate(comp)[::-1]):
  fs=list(sites(strand,f,False,gt,k))
  if not fs:continue
  for j,m,t in sites(strand,r,True,gt,k):
   for i,n,u in fs:
    gap=j-i-len(f);span=j+len(r)-i
    if (gap<50 or not 100<=span<=500) if target else not 0<=gap<=1000:continue
    ans[0]|=max(m,n)==0;ans[1]|=max(m,n)<=2 and t and u;ans[2]=True
 return ans
