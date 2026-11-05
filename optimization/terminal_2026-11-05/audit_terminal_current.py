from pathlib import Path
import csv,sys,json,html
from collections import Counter
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.Data.IUPACData import ambiguous_dna_values as I
ROOT=Path(__file__).resolve().parents[2];OUT=Path(__file__).resolve().parent/'current_terminal'
sys.path.insert(0,str(ROOT/'scripts'))
from independent import sites
read=lambda p:list(csv.DictReader(p.open()))
rows=read(ROOT/'catalog.csv');summaries=[];details=[];positions=[];cache={};recorded={}
for x in read(ROOT/'species.csv'):recorded[x['name'],x['accession']]=int(x['covered'])
def load(t,g,ref):
 key=t,g,ref
 if key in cache:return cache[key]
 p=ROOT/'panels'/t/'references'/g
 a={x.id:str(x.seq).upper() for x in SeqIO.parse(p/'aligned.fa','fasta')}
 if ref not in a:a={x.id:str(x.seq).upper() for x in SeqIO.parse(p/'coordinate_reference_alignment.fa','fasta')}
 assert ref in a
 raw={x.id:str(x.seq).upper() for x in SeqIO.parse(p/'target.fa','fasta')};meta=read(p/'target.csv')
 for m in meta:
  z=a[m['accession']].replace('-','');s=raw[m['accession']]
  assert z in (s,str(Seq(s).reverse_complement())),(key,m['accession'],'alignment identity')
 anchor=[i for i,b in enumerate(a[ref]) if b!='-'];cache[key]=(a,raw,meta,anchor);return cache[key]
def terminal(seq,cols,side,primer):
 # Reconstruct the complete homologous binding window before taking its 3' triplet.
 # Per-column triplets are misleading when compensating alignment gaps shift bases.
 z=seq[min(cols):max(cols)+1].replace('-','')
 if len(z)!=len(primer):return 'unresolved',''
 if side=='R':z=str(Seq(z).reverse_complement())
 b=z[-3:]
 if any(c not in 'ACGT' for c in b):return 'unresolved',b
 return ('match' if all(a in I[q] for a,q in zip(b,primer[-3:])) else 'mismatch'),b
def witnesses(s,f,r):
 valid=[];loose=[]
 for sign,strand in [('+',s),('-',str(Seq(s).reverse_complement()))]:
  fs=list(sites(strand,f,False,False,2));rs=list(sites(strand,r,True,False,2))
  for i,n,u in fs:
   for j,m,t in rs:
    if j-i-len(f)<50 or not 100<=j+len(r)-i<=500:continue
    fw=strand[i:i+len(f)];rv=str(Seq(strand[j:j+len(r)]).reverse_complement())
    fbad=sum(c not in I[q] for c,q in zip(fw[-3:],f[-3:]));rbad=sum(c not in I[q] for c,q in zip(rv[-3:],r[-3:]))
    assert u==(fbad==0) and t==(rbad==0)
    item=(n+m,fbad+rbad,sign,i+1,j+len(r),fw[-3:],rv[-3:])
    loose.append(item)
    if fbad==rbad==0:valid.append(item)
 return sorted(valid),sorted(loose)
for x in rows:
 a,raw,meta,anchor=load(x['target_taxon'],x['gene'],x['reference']);cols={e:[anchor[i-1] for i in range(int(x[e+'_start']),int(x[e+'_end'])+1)] for e in 'FR'}
 counts=Counter();basecounts={e:[Counter() for _ in range(3)] for e in 'FR'}
 for m in meta:
  acc=m['accession'];res={e:terminal(a[acc],cols[e],e,x[e]) for e in 'FR'}
  for e,(status,b) in res.items():
   counts[e+'_'+status]+=1
   if status in ['match','mismatch']:
    for i,c in enumerate(b):basecounts[e][i][c]+=1
  joint='match' if all(res[e][0]=='match' for e in 'FR') else ('unresolved' if any(res[e][0]=='unresolved' for e in 'FR') else 'incompatible')
  counts['both_'+joint]+=1
  valid,loose=witnesses(raw[acc],x['F'],x['R']);covered=bool(valid);assert covered==bool(recorded[x['name'],acc]),(x['name'],acc,'coverage discrepancy')
  counts['covered']+=covered;counts['relaxed_terminal']+=bool(loose);counts['terminal_excluded']+=bool(loose) and not covered
  w=(valid or loose or [None])[0]
  details.append(dict(name=x['name'],species=m['species'],accession=acc,family=m['family'],F_status=res['F'][0],F_target_3prime=res['F'][1],R_status=res['R'][0],R_target_3prime=res['R'][1],homologous_joint=joint,covered=int(covered),pair_with_terminal_constraint_removed=int(bool(loose)),witness_strand=w[2] if w else '',witness_start=w[3] if w else '',witness_end=w[4] if w else '',witness_F_terminal=w[5] if w else '',witness_R_terminal=w[6] if w else ''))
 assert counts['covered']==int(x['covered'])
 summaries.append(dict(name=x['name'],target=x['target_taxon'],gene=x['gene'],variant=x['variant'],N=len(meta),F_3prime=x['F'][-3:],R_3prime=x['R'][-3:],F_match=counts['F_match'],F_mismatch=counts['F_mismatch'],F_unresolved=counts['F_unresolved'],F_indel=counts['F_indel'],R_match=counts['R_match'],R_mismatch=counts['R_mismatch'],R_unresolved=counts['R_unresolved'],R_indel=counts['R_indel'],both_match=counts['both_match'],both_incompatible=counts['both_incompatible'],both_unresolved=counts['both_unresolved'],covered=counts['covered'],relaxed_terminal=counts['relaxed_terminal'],terminal_excluded=counts['terminal_excluded'],exact=x['exact']))
 for e in 'FR':
  for i,c in enumerate(basecounts[e]):positions.append(dict(name=x['name'],end=e,position_from_3prime=3-i,primer_base=x[e][-3+i],callable=sum(c.values()),**{b:c[b] for b in 'ACGT'},incompatible=sum(n for b,n in c.items() if b not in I[x[e][-3+i]])))
 print('PASS',x['name'],f"covered {counts['covered']}/{len(meta)}; homologous terminals {counts['both_match']}/{len(meta)}; excluded by terminal {counts['terminal_excluded']}",flush=True)
for filename,rs in [('summary.csv',summaries),('species.csv',details),('terminal_positions.csv',positions)]:
 with (OUT/filename).open('w') as h:
  w=csv.DictWriter(h,fieldnames=list(rs[0]),lineterminator='\n');w.writeheader();w.writerows(rs)
(OUT/'audit.json').write_text(json.dumps(dict(status='PASS',pairs=len(rows),oligos=2*len(rows),species_pair_evaluations=len(details),coverage_discrepancies=0,covered_terminal_mismatches=0,method='Homologous anchor projection plus independent unaligned paired-site scan; reverse primer checked in primer orientation; IUPAC compatible; no GT relaxation.'),indent=2))
