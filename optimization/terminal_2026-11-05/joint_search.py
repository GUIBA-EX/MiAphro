from pathlib import Path
import sys,csv,json,itertools,math,subprocess,shutil
from collections import Counter
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.Data.IUPACData import ambiguous_dna_values as I
R=Path(__file__).resolve().parents[2];O=Path(__file__).resolve().parent
sys.path.insert(0,str(R/'scripts'));from independent import sites,flags
from chemistry_rules import thermal,terminal,oligo_thermal
I={k:v for k,v in I.items() if k!='X'};inv={frozenset(v):k for k,v in I.items()}
read=lambda p:list(csv.DictReader(p.open()))
sel={x['name'] for x in read(O/'terminal_audit_baseline.csv') if x['variant']=='D' and int(x['terminal_excluded'])>0}
rows=[x for x in read(O/'baseline_catalog.csv') if x['name'] in sel];rows.sort(key=lambda x:x['target_taxon']!='Syllidae')
def deg(q):return math.prod(len(I[c]) for c in q)
def optimize(q,seqs):
 words={}
 for si,s in enumerate(seqs):
  for z in [s,str(Seq(s).reverse_complement())]:
   for start,_,_ in sites(z,q[:len(q)//2],False,False,2):
    w=z[start:start+len(q)]
    if len(w)==len(q) and set(w)<=set('ACGT'):words[w]=words.get(w,0)|(1<<si)
 ws=list(words);n=len(ws);valid=(1<<n)-1
 masks=[{b:sum(1<<i for i,w in enumerate(ws) if w[j]==b) for b in 'ACGT'} for j in range(len(q))]
 species_masks=[sum(1<<i for i,w in enumerate(ws) if words[w]>>si&1) for si in range(len(seqs))]
 def quality(v):
  a=b=c=valid
  for j,ch in enumerate(v):
   good=0
   for base in I[ch]:good|=masks[j][base]
   if j>=len(q)-3:a&=good;b&=good;c&=good
   else:c=(c&good)|(b&(valid^good));b=(b&good)|(a&(valid^good));a&=good
  cov=sum(1<<i for i,m in enumerate(species_masks) if m&c);ex=sum(1<<i for i,m in enumerate(species_masks) if m&a)
  return cov,ex
 seeds=set(''.join(x) for x in itertools.product(*(I[c] for c in q)))
 seeds.update(q[:len(q)//2]+w[len(q)//2:] for w in sorted(ws,key=lambda w:words[w].bit_count(),reverse=True)[:16])
 candidates={q}
 for seed in seeds:
  options=[]
  for j in range(len(q)//2,len(q)):
   bs=Counter({b:sum(words[w].bit_count() for w in ws if w[j]==b) for b in 'ACGT'})
   # Include binary and ternary choices, all observed terminal bases.
   opts=[(j,k) for k,v in I.items() if seed[j] in v and 1<len(v)<=4 and all(bs[b]>0 for b in v)]
   options.append(opts)
  candidates.add(seed)
  for count in [1,2,3]:
   for inds in itertools.combinations(range(len(options)),count):
    for changes in itertools.product(*(options[i] for i in inds)):
     if math.prod(len(I[c]) for _,c in changes)>8:continue
     z=list(seed)
     for j,c in changes:z[j]=c
     candidates.add(''.join(z))
 # Also preserve existing nonterminal ambiguity while exhaustively reallocating terminal masks.
 for seed in [q]+list(seeds):
  for trip in itertools.product(I,repeat=3):
   z=seed[:-3]+''.join(trip)
   if deg(z)<=8:candidates.add(z)
 result=[]
 for v in candidates:
  if not terminal(v):continue
  cov,ex=quality(v);result.append((v,cov,ex))
 key=lambda z:(z[1].bit_count(),z[2].bit_count(),-deg(z[0]),z[0])
 top=sorted(result,key=key,reverse=True)[:64]+sorted(result,key=lambda z:(z[2].bit_count(),z[1].bit_count(),-deg(z[0]),z[0]),reverse=True)[:32]
 top.append((q,*quality(q)));top=list({z[0]:z for z in top}.values())
 return top,len(candidates)
for x in rows:
 out=O/x['name'];p=R/'panels'/x['target_taxon']/'references'/x['gene']
 if (out/'joint_pairs.json').exists():continue
 seqs=[str(z.seq).upper() for z in SeqIO.parse(p/'target.fa','fasta')]
 fs,nf=optimize(x['F'],seqs);rs,nr=optimize(x['R'],seqs)
 print(x['name'],'oligos searched',nf,nr,'shortlist',len(fs),len(rs),flush=True)
 pairpool=[( (f[1]&r[1]).bit_count(),(f[2]&r[2]).bit_count(),f[0],r[0]) for f in fs for r in rs]
 choices=sorted(pairpool,reverse=True)[:100]+sorted(pairpool,key=lambda z:(z[1],z[0],z[2:]),reverse=True)[:60]
 # Include original local search candidates and baseline.
 pairs={(z[2],z[3]) for z in choices};pairs.update((z['F'],z['R']) for z in json.loads((out/'candidates.json').read_text()));pairs.add((x['F'],x['R']))
 results=[]
 for f,r in sorted(pairs):
  vv=[flags(s,f,r,True,False) for s in seqs];c=sum(z[1] for z in vv);e=sum(z[0] for z in vv)
  if (f,r)!=(x['F'],x['R']) and (c<int(x['covered']) or e<int(x['exact']) or (c,e)==(int(x['covered']),int(x['exact']))):continue
  ft,rt=oligo_thermal(f),oligo_thermal(r)
  ok=all(min(a[0],b[0])>=45 for a,b in zip(ft,rt)) and thermal(f,r)[0]
  results.append(dict(F=f,R=r,covered=c,exact=e,thermal_pass=ok,baseline=(f,r)==(x['F'],x['R'])))
 (out/'joint_pairs.json').write_text(json.dumps(results,indent=2))
 (out/'search_size.json').write_text(json.dumps(dict(F=nf,R=nr,pairs_screened=len(pairs))))
 print(x['name'],'target improvements',len(results)-1,'thermal',sum(z['thermal_pass'] and not z['baseline'] for z in results),flush=True)
print('JOINT SEARCH COMPLETE',flush=True)
