from pathlib import Path
import csv,json,io,subprocess,collections
import numpy as np
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.Data.IUPACData import ambiguous_dna_values as I
O=Path(__file__).resolve().parent;R=O.parents[1];root=R/'Syllidae';source=root/'16S'
read=lambda p:list(csv.DictReader(p.open()))
x=next(x for x in read(R/'Optimization3/results.csv') if x['objective']=='strict_profile' and x['name']=='Syllidae-16S-SD');F=x['F'];Q=x['R'];meta={x['accession']:x for x in read(source/'background.csv')};seq={x.id:str(x.seq).upper() for x in SeqIO.parse(source/'background.fa','fasta')}
old=root/x['source_pool']/'16S/S/independent_hits.tsv'
with old.open() as h:header=h.readline()
r=subprocess.run(['rg','--no-filename','^'+x['id']+r'\t',str(old)],capture_output=True,text=True,check=True)
prior={(z['model'],z['accession']):z for z in csv.DictReader(io.StringIO(header+r.stdout),delimiter='\t') if z['role']=='background'}
def sites(s,q,reverse,gt=False,k=3):
 # Exhaustive vectorized comparison at every sequence start, independent of seed and bit-parallel scanners.
 allowed=[set(I[c]) for c in q]
 if gt:
  for a in allowed:
   if 'G' in a:a.add('A')
   if 'T' in a:a.add('C')
 if reverse:allowed=[{str(Seq(b).complement()) for b in a} for a in allowed[::-1]]
 ar=np.frombuffer(s.encode(),dtype='S1');L=len(q)
 if len(ar)<L:return []
 windows=np.lib.stride_tricks.sliding_window_view(ar,L);mismatch=np.ones(windows.shape,dtype=bool)
 for j,a in enumerate(allowed):mismatch[:,j]=~np.isin(windows[:,j],np.array(sorted(a),dtype='S1'))
 valid=np.isin(windows,np.array(list('ACGT'),dtype='S1')).all(axis=1);inds=np.flatnonzero(valid&(mismatch.sum(axis=1)<=k));out=[]
 for i in inds:
  bad=np.flatnonzero(mismatch[i]);bad=sorted(L-1-int(j) for j in bad) if reverse else [int(j) for j in bad]
  binding=s[int(i):int(i)+L];binding=str(Seq(binding).reverse_complement()) if reverse else binding
  out.append((int(i),bad,binding))
 return out
products=[];records=[];model_counts=collections.Counter()
for acc,s in seq.items():
 modelpairs={}
 for gt in [False,True]:
  model='GT' if gt else 'WC';pairs=[]
  for orientation,strand in [('+',s),('-',str(Seq(s).reverse_complement()))]:
   fs=sites(strand,F,False,gt,2 if gt else 3);rs=sites(strand,Q,True,gt,2 if gt else 3)
   for i,fb,ft in fs:
    for j,rb,rt in rs:
     gap=j-i-len(F);span=j+len(Q)-i
     if not 0<=gap<=1000:continue
     f3=[len(F)-b for b in fb];r3=[len(Q)-b for b in rb]
     p=dict(accession=acc,model=model,orientation=orientation,product_start=i+1,product_end=j+len(Q),span=span,F_start=i+1,R_start=j+1,F_mm=len(fb),R_mm=len(rb),F_bad5=';'.join(str(b+1) for b in fb),R_bad5=';'.join(str(b+1) for b in rb),F_bad3=';'.join(map(str,sorted(f3))),R_bad3=';'.join(map(str,sorted(r3))),F_binding=ft,R_binding=rt,F_last1_clear=1 not in f3,R_last1_clear=1 not in r3,F_last3_clear=not any(b<=3 for b in f3),R_last3_clear=not any(b<=3 for b in r3),F_last5_clear=not any(b<=5 for b in f3),R_last5_clear=not any(b<=5 for b in r3))
     p['strict']=max(len(fb),len(rb))<=2 and p['F_last3_clear'] and p['R_last3_clear'];p['exact']=len(fb)==len(rb)==0
     pairs.append(p)
  modelpairs[model]=pairs
  e=any(p['exact'] for p in pairs);strict=any(p['strict'] for p in pairs);loose=bool(pairs);priorrow=prior.get((model,acc),{})
  assert [e,strict,loose]==[priorrow.get(k)=='1' for k in ['exact','le2last3','le2' if gt else 'le3']],(acc,model)
  model_counts[model+'_loose']+=loose;model_counts[model+'_strict']+=strict;model_counts[model+'_exact']+=e
 ps=modelpairs['WC']
 if not ps:continue
 products+=ps
 flags={f'any_both_last{k}_clear':any(p[f'F_last{k}_clear'] and p[f'R_last{k}_clear'] for p in ps) for k in [1,3,5]}
 strict=any(p['strict'] for p in ps)
 category='A_strict' if strict else ('B_last3_clear_three_mismatch' if flags['any_both_last3_clear'] else ('C_last1_clear_internal_last3_mismatch' if flags['any_both_last1_clear'] else 'D_terminal_mismatch_in_every_product'))
 # Representative most permissive product prioritizes intact 3prime ends, then fewer mismatches.
 best=min(ps,key=lambda p:(not(p['F_last1_clear'] and p['R_last1_clear']),not(p['F_last3_clear'] and p['R_last3_clear']),p['F_mm']+p['R_mm'],max(p['F_mm'],p['R_mm']),abs(p['span']-156)))
 records.append(dict(**meta[acc],category=category,products=len(ps),strict=strict,**flags,any_both_le2=any(max(p['F_mm'],p['R_mm'])<=2 for p in ps),any_span100_200=any(100<=p['span']<=200 for p in ps),any_span100_500=any(100<=p['span']<=500 for p in ps),any_last3_clear_span100_200=any(p['F_last3_clear'] and p['R_last3_clear'] and 100<=p['span']<=200 for p in ps),any_last3_clear_span100_500=any(p['F_last3_clear'] and p['R_last3_clear'] and 100<=p['span']<=500 for p in ps),min_span=min(p['span'] for p in ps),max_span=max(p['span'] for p in ps),**{'representative_'+k:best[k] for k in ['span','F_mm','R_mm','F_bad3','R_bad3','F_binding','R_binding']}))
def write(name,rows):
 with (O/name).open('w') as h:w=csv.DictWriter(h,rows[0]);w.writeheader();w.writerows(rows)
write('hits_1007.csv',records);write('all_WC_products.csv',products)
groups=[]
for level in ['group','family']:
 totals=collections.Counter(z[level] for z in meta.values())
 for taxon,n in collections.Counter(z[level] for z in records).most_common():
  subset=[z for z in records if z[level]==taxon];groups.append(dict(level=level,taxon=taxon,background_species=totals[taxon],loose_hits=n,strict_hits=sum(z['strict'] for z in subset),last1_clear=sum(z['any_both_last1_clear'] for z in subset),last3_clear=sum(z['any_both_last3_clear'] for z in subset),last5_clear=sum(z['any_both_last5_clear'] for z in subset),S_span=sum(z['any_span100_200'] for z in subset)))
write('taxonomy.csv',groups)
positions=[]
for end,q in [('F',F),('R',Q)]:
 for j in range(1,len(q)+1):
  accs={p['accession'] for p in products if str(j) in p[end+'_bad3'].split(';')};rep=sum(str(j) in z['representative_'+end+'_bad3'].split(';') for z in records)
  positions.append(dict(primer=end,position_from_3prime=j,position_from_5prime=len(q)-j+1,primer_code=q[-j],records_any_product_mismatch=len(accs),representative_product_mismatch=rep))
write('mismatch_positions.csv',positions)
summary=dict(primers={'F':F,'R':Q},background_species=len(meta),counts=dict(model_counts),WC_products=len(products),categories=dict(collections.Counter(z['category'] for z in records)),record_flags={k:sum(z[k] for z in records) for k in records[0] if k.startswith('any_')},all_products_have_last3_mismatch=sum(not z['any_both_last3_clear'] for z in records),validation='All 3675 records WC and GT match prior independent flags')
(O/'summary.json').write_text(json.dumps(summary,indent=2));print(json.dumps(summary,indent=2));print(groups[:12])
