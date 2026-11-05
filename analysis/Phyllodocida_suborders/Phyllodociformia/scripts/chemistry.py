import primer3,itertools,functools,math
from Bio.Data.IUPACData import ambiguous_dna_values as I
kw=dict(mv_conc=50,dv_conc=1.5,dntp_conc=.6)
@functools.lru_cache(None)
def expand(q):return list(map(''.join,itertools.product(*(I[b] for b in q))))
def deg(q):return math.prod(len(I[c]) for c in q)
def terminal(q):
 if set.intersection(*(set(I[c]) for c in q[-3:])):return False
 if not any(set(I[c])<=set('GC') for c in q[-3:]):return False
 if not 1<=sum(set(I[c])<=set('GC') for c in q[-5:]):return False
 if sum(bool(set(I[c])&set('GC')) for c in q[-5:])>4:return False
 if any(set.intersection(*(set(I[c]) for c in q[i:i+5])) for i in range(len(q)-4)):return False
 return True
@functools.lru_cache(None)
def single(q,c):return primer3.calc_tm(q,dna_conc=c,**kw),max(primer3.calc_hairpin(q,dna_conc=c,**kw).tm,primer3.calc_homodimer(q,dna_conc=c,**kw).tm)
@functools.lru_cache(None)
def hetero(a,b,c):return primer3.calc_heterodimer(a,b,dna_conc=c,**kw).tm
@functools.lru_cache(None)
def oligo_thermal(q):
 vs=expand(q);reports=[]
 for c in [50,200/len(vs)]:
  v=[single(s,c) for s in vs];reports.append((min(x[0] for x in v),max(x[0] for x in v),max(x[1] for x in v)))
 return reports
@functools.lru_cache(None)
def thermal(F,R):
 fs,rs=expand(F),expand(R);report=[]
 for scenario in ['50_per_variant','200_per_direction']:
  fc=50 if scenario=='50_per_variant' else 200/len(fs);rc=50 if scenario=='50_per_variant' else 200/len(rs);ft=[single(s,fc) for s in fs];rt=[single(s,rc) for s in rs];st=max(x[1] for x in ft+rt)
  for a,b in itertools.combinations(fs,2):st=max(st,hetero(a,b,fc))
  for a,b in itertools.combinations(rs,2):st=max(st,hetero(a,b,rc))
  for a,b in itertools.product(fs,rs):st=max(st,hetero(a,b,max(fc,rc)))
  report.append(dict(scenario=scenario,F_tm_min=min(x[0] for x in ft),F_tm_max=max(x[0] for x in ft),R_tm_min=min(x[0] for x in rt),R_tm_max=max(x[0] for x in rt),structure_max=st))
 return all(x['structure_max']<55 and abs((x['F_tm_min']+x['F_tm_max']-x['R_tm_min']-x['R_tm_max'])/2)<=5 for x in report),report
