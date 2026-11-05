"""Build the unified catalogue from frozen, independently checked source results."""
from pathlib import Path
import csv,json,math,hashlib,shutil
from collections import Counter
from Bio import SeqIO
from Bio.Data.IUPACData import ambiguous_dna_values as I
R=Path(__file__).resolve().parents[1];A=R/'analysis';P=R/'panels';P.mkdir(exist_ok=True)
read=lambda p:list(csv.DictReader(p.open()))
def write(p,rows):
 p.parent.mkdir(parents=True,exist_ok=True)
 with p.open('w') as h:
  w=csv.DictWriter(h,rows[0],lineterminator='\n');w.writeheader();w.writerows(rows)
units=[('Aphroditiformia','suborder',R/'studies/Aphroditiformia','SML'),('Phyllodocida','order',A/'Phyllodocida','SML')]+[(n,'suborder',A/'Phyllodocida_suborders'/n,'Round2' if n=='Nereidiformia' else 'Expanded') for n in ['Glyceriformia','Nereidiformia','Phyllodociformia']]+[(n,'family',A/'Phyllodocida_families'/n,'Round2') for n in ['Syllidae','Nereididae','Hesionidae']]
opt=A/'Phyllodocida_families/Optimization3';new=read(opt/'results.csv');newdetails=read(opt/'species.csv');catalog=[];conditional=[];allsp=[];condsp=[];inventory=[];mapping=[]
fields=['name','target_taxon','rank','gene','band','variant','profile','reference','F_start','F_end','R_start','R_end','span','F','R','deg_F','deg_R','target_species','covered','exact','background_species','WC','GT','WCloose','GTloose','source_name','source_results','source_pool','source_id']
def normalized(x,taxon,rank,profile,path):
 gene,band=x['name'].split('-')[-2:];d=band.endswith('D');length=band.removesuffix('D');slot=length+('2' if profile=='conditional' else '')+('D' if d else '');name=f'MiPhyllo-{taxon}-{gene}-{slot}'
 z=dict(name=name,target_taxon=taxon,rank=rank,gene=gene,band=length,variant='D' if d else 'ordinary',profile=profile,reference=x['reference'],span=int(x['R_end'])-int(x['F_start'])+1,source_name=x['name'],source_results=str(path.relative_to(R)),source_pool=x.get('source_pool','SML' if taxon in ['Aphroditiformia','Phyllodocida'] else 'Expanded'),source_id=x['id'])
 for k in ['F_start','F_end','R_start','R_end','F','R','target_species','covered','exact','background_species','WC','GT','WCloose','GTloose']:z[k]=x[k]
 for end in 'FR':z['deg_'+end]=math.prod(len(I[b]) for b in x[end])
 return {k:z[k] for k in fields}
for taxon,rank,root,folder in units:
 basepath=root/folder/'results.csv';base=read(basepath);details=read(root/folder/'species.csv');selected={x['name']:(x,basepath,details) for x in base}
 if rank=='family':
  for x in new:
   if x['objective']=='balanced' and x['name'].startswith(taxon+'-'):selected[x['name']]=(x,opt/'results.csv',[s for s in newdetails if s['objective']=='balanced'])
 rows=[];sps=[]
 for oldname,(x,path,ds) in selected.items():
  z=normalized(x,taxon,rank,'standard',path);rows.append(z);mapping.append(dict(profile='standard',source_name=oldname,name=z['name'],source_results=z['source_results']))
  sps.extend(dict(name=z['name'],species=s['species'],accession=s['accession'],family=s['family'],covered=s['covered'],exact=s['exact']) for s in ds if s['name']==oldname)
 rows.sort(key=lambda x:(x['gene'],['S','M','L'].index(x['band']),x['variant']=='D'))
 d=P/taxon;d.mkdir(exist_ok=True);write(d/'results.csv',rows);write(d/'species.csv',sps);(d/'primers.fasta').write_text(''.join(f">{x['name']}_{e}\n{x[e]}\n" for x in rows for e in 'FR'));catalog+=rows;allsp+=sps
 for gene in ['16S','28S']:
  dest=d/'references'/gene;dest.mkdir(parents=True,exist_ok=True)
  for filename in ['target.fa','target.csv','background.fa','background.csv','aligned.fa']:shutil.copy2(root/gene/filename,dest/filename)
  # Preserve coordinate anchors outside a focal-family target panel.
  control=root/'Controls'/gene/'coordinate_reference_alignment.fa'
  if control.exists():shutil.copy2(control,dest/'coordinate_reference_alignment.fa')
  inventory.append(dict(target_taxon=taxon,rank=rank,gene=gene,target_species=len(read(dest/'target.csv')),background_species=len(read(dest/'background.csv')),source_root=str(root.relative_to(R))))
 if rank=='family':
  for band in ['S','M','L']:
   for gene in ['16S','28S']:
    alts=[x for x in new if x['objective']=='strict_profile' and x['name'] in [f'{taxon}-{gene}-{band}',f'{taxon}-{gene}-{band}D']]
    if len(alts)!=2:continue
    if all(x['F']==selected[x['name']][0]['F'] and x['R']==selected[x['name']][0]['R'] for x in alts):continue
    for x in alts:
     z=normalized(x,taxon,rank,'conditional',opt/'results.csv');conditional.append(z);mapping.append(dict(profile='conditional',source_name=x['name'],name=z['name'],source_results=z['source_results']))
     condsp.extend(dict(name=z['name'],species=s['species'],accession=s['accession'],family=s['family'],covered=s['covered'],exact=s['exact']) for s in newdetails if s['objective']=='strict_profile' and s['name']==x['name'])
write(R/'catalog.csv',catalog);write(R/'species.csv',allsp);write(R/'panels.csv',inventory);write(R/'provenance.csv',mapping)
(R/'primers.fasta').write_text(''.join(f">{x['name']}_{e}\n{x[e]}\n" for x in catalog for e in 'FR'))
write(R/'conditional/results.csv',conditional);write(R/'conditional/species.csv',condsp);(R/'conditional/primers.fasta').write_text(''.join(f">{x['name']}_{e}\n{x[e]}\n" for x in conditional for e in 'FR'))
print('Catalogue',len(catalog),'pairs;',len(conditional),'conditional pairs;',len(units),'target units')
