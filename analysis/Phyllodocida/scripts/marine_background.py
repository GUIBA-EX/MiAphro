"""Assemble traceable marine-priority background; verify genus habitat with WoRMS."""
from pathlib import Path
from Bio import SeqIO
import hashlib
import csv,json,urllib.request,urllib.parse,concurrent.futures,time,re,gzip,io
root=Path(__file__).resolve().parents[1];outputs=root.parent;cache=root/'marine_taxonomy';cache.mkdir(exist_ok=True)
items=[]
def add(gene,seq,meta=None,group=None):
 md={r['accession']:r for r in meta} if meta else {}
 for r in SeqIO.parse(seq,'fasta'):
  m=md.get(r.id,{});sp=m.get('organism',m.get('species',r.description[len(r.id)+1:])).strip();s=str(r.seq).upper()
  if sp.startswith(r.id+' '):sp=sp[len(r.id)+1:]
  if not re.fullmatch('[A-Z][a-z-]+ [a-z][a-z-]+',sp) or len(s)<150 or sum(c not in 'ACGT' for c in s)>.03*len(s):continue
  items.append(dict(gene=gene,accession=r.id,species=sp,sequence=s,source=str(seq.relative_to(outputs)),source_group=group or m.get('clade','')))
add('16S',outputs/'MiLily/references/background.fasta',list(csv.DictReader((outputs/'MiLily/references/metadata.tsv').open(),delimiter='\t')))
add('16S',outputs/'MiLily/references/target.fasta',list(csv.DictReader((outputs/'MiLily/references/metadata.tsv').open(),delimiter='\t')))
for role in ['target','background']:add('16S',outputs/f'Actiniaria/16S/{role}.fasta',list(csv.DictReader((outputs/f'Actiniaria/16S/{role}_panel.csv').open())))
for f in (outputs/'MiUrchin_28S_trial').glob('*.fa'):
 if f.stem not in ['fragments','ref','Polychaeta']:add('28S',f,group=f.stem)
for role in ['target','background']:add('28S',outputs/f'MiAnemone_28S/{role}.fasta',list(csv.DictReader((outputs/'MiAnemone_28S/panel.csv').open())))
# Newly retrieved non-target annelids: exclude terrestrial/freshwater clitellates by WoRMS marine status later.
for gene in ['16S','28S']:
 add(gene,root/gene/'background.fa',list(csv.DictReader((root/gene/'background.csv').open())),group='Annelida')
existing={}
for p in cache.glob('*.json'):
 d=json.loads(p.read_text())
 for n,v in zip(d['names'],d['records']):existing[n]=v
supp=root/'supplement/16S/background.fa'
if supp.exists():add('16S',supp,list(csv.DictReader(supp.with_suffix('.csv').open())),group='intertidal supplement')
genera=sorted({r['species'].split()[0] for r in items}-set(existing))
def fetch(start):
 names=genera[start:start+40];p=cache/('batch_'+hashlib.sha256('|'.join(names).encode()).hexdigest()[:12]+'.json');q=urllib.parse.urlencode([('scientificnames[]',n) for n in names]+[('marine_only','false')]);url='https://www.marinespecies.org/rest/AphiaRecordsByNames?'+q
 if not p.exists():
  for k in range(5):
   try:
    with urllib.request.urlopen(url,timeout=100) as h:res=json.load(h)
    p.write_text(json.dumps(dict(names=names,records=res),indent=2));break
   except Exception:
    if k==4:raise
    time.sleep(3)
 d=json.loads(p.read_text());assert d['names']==names
 print('WoRMS',start,len(genera),flush=True);return d
with concurrent.futures.ThreadPoolExecutor(max_workers=2) as pool:data=list(pool.map(fetch,range(0,len(genera),40)))
data.append(dict(names=list(existing),records=list(existing.values())))
lookup={}
for d in data:
 for name,hits in zip(d['names'],d['records']):
  hits=[r for r in (hits or []) if r.get('rank')=='Genus' and r.get('scientificname')==name and r.get('kingdom')=='Animalia' and r.get('isMarine')==1]
  if len(hits)==1:lookup[name]=hits[0]
excluded=[]
for gene in ['16S','28S']:
 best={};target={r['species'] for r in csv.DictReader((root/gene/'target.csv').open())}
 for r in items:
  if r['gene']!=gene or r['species'] in target:continue
  tax=lookup.get(r['species'].split()[0])
  if not tax or tax.get('order')=='Phyllodocida':excluded.append(dict(gene=gene,species=r['species'],accession=r['accession'],reason='marine genus not uniquely confirmed or target order'));continue
  r=dict(r,group=tax['phylum'],family=tax['family'] or 'unknown',worms_genus_id=tax['AphiaID'],worms_url=tax['url']);sp=r['species']
  if sp not in best or sum(c in 'ACGT' for c in r['sequence'])>sum(c in 'ACGT' for c in best[sp]['sequence']):best[sp]=r
 d=root/gene;rows=sorted(best.values(),key=lambda r:(r['group'],r['species']));(d/'marine_background.fa').write_text(''.join('>'+r['accession']+'\n'+r['sequence']+'\n' for r in rows))
 with (d/'marine_background.csv').open('w') as h:
  w=csv.DictWriter(h,[k for k in rows[0] if k!='sequence']);w.writeheader();w.writerows({k:v for k,v in r.items() if k!='sequence'} for r in rows)
 print(gene,'marine background',len(rows),flush=True)
with (cache/'excluded.csv').open('w') as h:w=csv.DictWriter(h,excluded[0]);w.writeheader();w.writerows(excluded)
