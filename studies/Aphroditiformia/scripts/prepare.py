from pathlib import Path
import json,time,urllib.request,urllib.parse,gzip,io,csv
from Bio import SeqIO
root=Path(__file__).resolve().parents[1]
def api(method,**kw):
 for i in range(5):
  try:
   time.sleep(.4)
   return urllib.request.urlopen('https://eutils.ncbi.nlm.nih.gov/entrez/eutils/'+method+'.fcgi?'+urllib.parse.urlencode(dict(db='nuccore',**kw)),timeout=80).read()
  except Exception:
   if i==4:raise
   time.sleep(2)
for gene in ['16S','28S']:
 for role,taxa in [('target','Aphroditiformia Polynoidae Sigalionidae Aphroditidae Acoetidae Iphionidae Eulepethidae Pholoidae Pisionidae'),('background','Nereididae Phyllodocidae Glyceridae Nephtyidae Hesionidae Syllidae')]:
  d=root/gene;d.mkdir(exist_ok=True);q='('+' OR '.join(x+'[Organism]' for x in taxa.split())+') AND '+gene+'[Title] AND 150:6000[SLEN]';p=d/(role+'.json')
  if not p.exists():p.write_bytes(api('esearch',term=q,retmode='json',retmax=10000))
  ids=json.loads(p.read_text())['esearchresult']['idlist'];print(gene,role,len(ids),flush=True)
  best={}
  for i in range(0,len(ids),150):
   f=d/f'{role}_{i}.gb.gz'
   if not f.exists():f.write_bytes(gzip.compress(api('efetch',rettype='gb',retmode='text',id=','.join(ids[i:i+150]))))
   for r in SeqIO.parse(io.StringIO(gzip.decompress(f.read_bytes()).decode()),'genbank'):
    sp=r.annotations['organism'];line=r.annotations['taxonomy']
    if len(sp.split())!=2 or not sp.split()[1].islower() or '.' in sp or any(z in sp for z in ['uncultured','unidentified']):continue
    feats=[x for x in r.features if x.type=='rRNA' and gene.lower() in ' '.join(x.qualifiers.get('product',[])).lower()]
    if not feats:continue
    seq=str(max(feats,key=lambda x:len(x.location)).extract(r.seq)).upper()
    if not 150<=len(seq)<=5500 or len(seq.strip('ACGT'))>len(seq)*.03 and sum(c not in 'ACGT' for c in seq)>len(seq)*.03:continue
    fam=next((z for z in reversed(line) if z.endswith('idae')),'unknown')
    item=(sum(c in 'ACGT' for c in seq),r.id,seq,fam)
    if sp not in best or item>best[sp]:best[sp]=item
  (d/(role+'.fa')).write_text(''.join(f'>{a}\n{s}\n' for _,a,s,_ in best.values()))
  with (d/(role+'.csv')).open('w') as h:
   w=csv.writer(h);w.writerow(['species','accession','family']);w.writerows((sp,a,f) for sp,(_,a,s,f) in best.items())
