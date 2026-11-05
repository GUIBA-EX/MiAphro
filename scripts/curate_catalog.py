"""Reproducible editorial pruning; original sequences and scores are unchanged."""
from pathlib import Path
import csv,json
R=Path(__file__).resolve().parents[1]
read=lambda p:list(csv.DictReader(p.open()))
def write(p,rows,fields):
 with p.open('w') as h:
  w=csv.DictWriter(h,fieldnames=fields,lineterminator='\n');w.writeheader();w.writerows(rows)
def fasta(p,rs):p.write_text(''.join(f">{x['name']}_{e}\n{x[e]}\n" for x in rs for e in 'FR'))
# Always start with the preserved original export. The source builder must reproduce it.
A=R/'archive/2026-11-05-unfiltered';full=read(A/'catalog.csv');alts=read(A/'conditional/results.csv')
keep=set();decisions=[]
for x in full+alts:
 if x['variant']!='D':continue
 reasons=[]
 if x['target_taxon'] in ['Phyllodocida','Nereidiformia']:reasons.append('Broad-target unit withdrawn')
 if int(x['covered'])/int(x['target_species'])<.5:reasons.append('D coverage below 50%')
 if int(x['exact'])==0:reasons.append('No D exact matches')
 if x['profile']=='conditional':reasons.append('Relaxed loose-background gate; not in core catalogue')
 if not reasons:keep.update([x['name'],x['name'][:-1]])
 decisions.append(dict(set=x['name'][:-1],target_taxon=x['target_taxon'],gene=x['gene'],band=x['band'],profile=x['profile'],covered=x['covered'],exact=x['exact'],target_species=x['target_species'],status='archived' if reasons else 'retained',reason='; '.join(reasons) if reasons else 'D coverage >=50%, exact >0, original full background profile; ordinary retained as same-coordinate comparator'))
rows=[x for x in full if x['name'] in keep];taxa={x['target_taxon'] for x in rows}
write(R/'catalog.csv',rows,list(full[0]));fasta(R/'primers.fasta',rows)
sp=read(A/'species.csv');write(R/'species.csv',[x for x in sp if x['name'] in keep],list(sp[0]))
ps=read(A/'panels.csv');write(R/'panels.csv',[x for x in ps if x['target_taxon'] in taxa],list(ps[0]))
mp=read(A/'provenance.csv');write(R/'provenance.csv',[x for x in mp if x['name'] in keep],list(mp[0]))
for t in {x['target_taxon'] for x in full}:
 p=R/'panels'/t;rs=[x for x in rows if x['target_taxon']==t]
 write(p/'results.csv',rs,list(full[0]));write(p/'species.csv',[x for x in sp if x['name'] in keep and x['name'].startswith('MiPhyllo-'+t+'-')],list(sp[0]));fasta(p/'primers.fasta',rs)
write(R/'conditional/results.csv',[],list(alts[0]));write(R/'conditional/species.csv',[],list(sp[0]));fasta(R/'conditional/primers.fasta',[])
write(R/'archive/selection_decisions.csv',decisions,list(decisions[0]))
removed=[x for x in full+alts if x['name'] not in keep];write(R/'archive/removed_candidates.csv',removed,list(full[0]));fasta(R/'archive/removed_primers.fasta',removed)
# Preserve the complete independent scan; record exact identity of the retained subset.
verified=json.loads((A/'audit/catalog_rescan.json').read_text());assert keep <= {x['name'] for x in verified['records']}
(R/'audit/curation.json').write_text(json.dumps(dict(status='PASS',date='2026-11-05',original_pairs=len(full),original_conditional_pairs=len(alts),current_pairs=len(rows),archived_pairs=len(removed),sets=len(rows)//2,units=sorted(taxa),sequence_changes=0,rule='Retain standard ordinary/D sets with D coverage >=50% and exact >0; exclude Phyllodocida and Nereidiformia; archive all relaxed-background alternatives.',retained_names=sorted(keep),independent_verification='Unchanged subset of audit/catalog_rescan.json'),indent=2)+'\n')
print(f'Curated: {len(rows)} pairs, {len(removed)} archived, {len(taxa)} units')
