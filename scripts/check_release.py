"""Structural release gate for the integrated MiPhyllo catalogue."""
from pathlib import Path
import csv,json,re,math,hashlib,sys
from collections import Counter
from Bio import SeqIO
from Bio.Data.IUPACData import ambiguous_dna_values as I
R=Path(__file__).resolve().parents[1]
read=lambda p:list(csv.DictReader(p.open()))
rows=read(R/'catalog.csv');alt=read(R/'conditional/results.csv');assert len(rows)==48 and len(alt)==0
assert len({x['name'] for x in rows+alt})==48
verified=json.loads((R/'audit/catalog_rescan.json').read_text());assert verified['status']=='PASS' and verified['pairs']==102
assert {(x['name'],x['profile']) for x in rows+alt}<={(x['name'],x['profile']) for x in verified['records']}
for rs,folder in [(rows,R),(alt,R/'conditional')]:
 seq={x.id:str(x.seq) for x in SeqIO.parse(folder/'primers.fasta','fasta')};assert seq=={x['name']+'_'+side:x[side] for x in rs for side in 'FR'}
 details=read(folder/'species.csv');lookup={x['name']:x for x in rs}
 for x in rs:
  expected=f"MiPhyllo-{x['target_taxon']}-{x['gene']}-{x['band']}"+('2' if x['profile']=='conditional' else '')+('D' if x['variant']=='D' else '')
  assert x['name']==expected
  for side in 'FR':
   assert len(x[side])==int(x[side+'_end'])-int(x[side+'_start'])+1
   assert math.prod(len(I[b]) for b in x[side])==int(x['deg_'+side])<=8
   assert set(x[side][:len(x[side])//2])<=set('ACGT')
   if x['variant']=='ordinary':assert set(x[side])<=set('ACGT')
  if x['variant']=='D':
   a=lookup[x['name'][:-1]];assert all(a[k]==x[k] for k in ['reference','F_start','F_end','R_start','R_end'])
  ss=[s for s in details if s['name']==x['name']];assert len(ss)==int(x['target_species'])==len({s['species'] for s in ss})
  for k in ['covered','exact']:assert sum(int(s[k]) for s in ss)==int(x[k])
  docroot=R/'conditional' if x['profile']=='conditional' else R/'panels'/x['target_taxon']
  for lang in ['README.md','README.en.md']:
   line=next(t for t in (docroot/lang).read_text().splitlines() if t.startswith('|'+x['name']+'|') and '`' in t)
   assert all('`'+x[k]+'`' in line for k in ['F','R']) and f"{x['covered']}/{x['target_species']}" in line and f"{x['exact']}/{x['target_species']}" in line
for x in read(R/'panels.csv'):
 p=R/'panels'/x['target_taxon']/'references'/x['gene'];species={}
 for role in ['target','background']:
  metadata=read(p/(role+'.csv'));seq=list(SeqIO.parse(p/(role+'.fa'),'fasta'));assert len(seq)==len(metadata)==len({z['species'] for z in metadata})==int(x[role+'_species']);assert {z.id for z in seq}=={z['accession'] for z in metadata};species[role]={z['species'] for z in metadata}
 assert not species['target']&species['background']
for p in R.rglob('*.md'):
 if '.git' in p.parts or '.recheck' in p.parts:continue
 for link in re.findall(r'\]\(([^)]+)\)',p.read_text()):
  if not link.startswith(('https://','http://','#')):assert (p.parent/link.split('#')[0]).is_file(),(p,link)
for x in json.loads((R/'analysis/source_manifest.json').read_text()):assert hashlib.sha256((R/x['path']).read_bytes()).hexdigest()==x['sha256'],x['path']
assert (R/'figures/catalog_coverage.svg').exists()
if (R/'SHA256SUMS').exists():
 for line in (R/'SHA256SUMS').read_text().splitlines():
  digest,path=line.split('  ',1);assert hashlib.sha256((R/path).read_bytes()).hexdigest()==digest,path
print('PASS: 48 current pairs; 54 archived pairs; unique names; FASTA; same endpoints; degeneracy; panels; bilingual tables; links; source hashes; independent rescans.')

# Selection identity and exclusion checks are independent of the export builder.
original=read(R/'archive/2026-11-05-unfiltered/catalog.csv')
original_alt=read(R/'archive/2026-11-05-unfiltered/conditional/results.csv')
expected=set()
for x in original:
 if x['variant']=='D' and x['target_taxon'] not in {'Phyllodocida','Nereidiformia'} and 2*int(x['covered'])>=int(x['target_species']) and int(x['exact'])>0:
  expected.update((x['name'],x['name'][:-1]))
assert {x['name'] for x in rows}==expected
updates={x['name']:x for x in read(R/'optimization/terminal_2026-11-05/selected_results.csv')}
assert len(updates)==4
assert rows==[updates.get(x['name'],x) for x in original if x['name'] in expected]
comparisons=json.loads((R/'optimization/terminal_2026-11-05/comparison.json').read_text())
assert len(comparisons)==9 and all(x['independent_status']=='PASS' for x in comparisons)
for x in comparisons:
 if x['status']=='replace':
  assert x['after']['thermal_pass'] and x['after']['background_pass']
  assert x['after']['covered']>=x['before']['covered'] and x['after']['exact']>=x['before']['exact']
  assert any(x['after'][k]>x['before'][k] for k in ['covered','exact'])
  assert all(str(x['after'][k])==updates[x['name']][k] for k in ['F','R','covered','exact','WC','GT','WCloose','GTloose'])
removed=read(R/'archive/removed_candidates.csv')
assert {x['name'] for x in removed}=={x['name'] for x in original+original_alt}-expected
assert len(removed)==54
for p in (R/'panels').iterdir():
 if not p.is_dir():continue
 expected_rows=[x for x in rows if x['target_taxon']==p.name]
 assert read(p/'results.csv')==expected_rows
 assert {x.id:str(x.seq) for x in SeqIO.parse(p/'primers.fasta','fasta')}=={x['name']+'_'+e:x[e] for x in expected_rows for e in 'FR'}
manifest=json.loads((R/'publications/MiPhyllo_manifest.json').read_text())
assert manifest['current_pairs']==48 and manifest['conditional_pairs']==0 and len(manifest['units'])==6
qa=json.loads((R/'audit/integrated_preprint.json').read_text())
assert qa['status']=='PASS' and qa['pairs']==48 and qa['oligos_checked']==96
assert qa['pdf_sha256']==hashlib.sha256((R/'publications/MiPhyllo_preprint.pdf').read_bytes()).hexdigest()
print('PASS: editorial criteria, archive partition, panel exports and preprint identity.')
