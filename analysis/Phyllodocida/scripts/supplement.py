from pathlib import Path
import sys
import concurrent.futures
import json,time,urllib.request,urllib.parse,gzip,io,csv
from Bio import SeqIO
root=Path(__file__).resolve().parents[1]/'supplement';root.mkdir(exist_ok=True)
def api(method,**kw):
 for i in range(5):
  try:
   time.sleep(.4)
   return urllib.request.urlopen('https://eutils.ncbi.nlm.nih.gov/entrez/eutils/'+method+'.fcgi?'+urllib.parse.urlencode(dict(db='nuccore',**kw)),timeout=80).read()
  except Exception:
   if i==4:raise
   time.sleep(2)
for gene in (sys.argv[1:] or ['16S','28S']):
 for role,taxa in [('background','(Carcinus[Organism] OR Portunus[Organism] OR Charybdis[Organism] OR Scylla[Organism] OR Cancer[Organism] OR Pagurus[Organism] OR Palaemon[Organism] OR Alpheus[Organism] OR Penaeus[Organism] OR Litopenaeus[Organism] OR Balanus[Organism] OR Chthamalus[Organism] OR Amphibalanus[Organism] OR Gammarus[Organism] OR Caprella[Organism] OR Crassostrea[Organism] OR Magallana[Organism] OR Mytilus[Organism] OR Perna[Organism] OR Ruditapes[Organism] OR Meretrix[Organism] OR Patella[Organism] OR Littorina[Organism] OR Nucella[Organism] OR Haliotis[Organism] OR Aplysia[Organism] OR Octopus[Organism] OR Sepia[Organism] OR Loligo[Organism] OR Haliclona[Organism] OR Halichondria[Organism] OR Suberites[Organism] OR Cliona[Organism] OR Hymeniacidon[Organism] OR Bugula[Organism] OR Bugulina[Organism] OR Membranipora[Organism] OR Watersipora[Organism] OR Ciona[Organism] OR Botryllus[Organism])')]:
  d=root/gene;d.mkdir(exist_ok=True);q='('+taxa+') AND ('+gene+'[Title] OR '+gene+'[Gene Name]) AND 150:6000[SLEN]';p=d/(role+'.json')
  if not p.exists():p.write_bytes(api('esearch',term=q,retmode='json',retmax=10000))
  result=json.loads(p.read_text())['esearchresult'];ids=result['idlist'];
  for start in range(len(ids),int(result['count']),10000):
   pg=d/f'{role}_search_{start}.json'
   if not pg.exists():pg.write_bytes(api('esearch',term=q,retmode='json',retmax=10000,retstart=start))
   ids+=json.loads(pg.read_text())['esearchresult']['idlist']
  assert len(ids)==int(result['count'])
  (d/(role+'_query.txt')).write_text(q+'\nRetrieved 2026-11-05\n');print(gene,role,len(ids),flush=True)
  def fetch(i):
   f=d/f'{role}_{i}.gb.gz'
   if not f.exists():f.write_bytes(gzip.compress(api('efetch',rettype='gb',retmode='text',id=','.join(ids[i:i+150]))))
  with concurrent.futures.ThreadPoolExecutor(max_workers=2) as pool:list(pool.map(fetch,range(0,len(ids),150)))
  best={}
  for i in range(0,len(ids),150):
   f=d/f'{role}_{i}.gb.gz'
   if not f.exists():f.write_bytes(gzip.compress(api('efetch',rettype='gb',retmode='text',id=','.join(ids[i:i+150]))))
   for r in SeqIO.parse(io.StringIO(gzip.decompress(f.read_bytes()).decode()),'genbank'):
    sp=r.annotations['organism'];line=r.annotations['taxonomy']
    if len(sp.split())!=2 or not sp.split()[0][0].isupper() or not sp.split()[1].islower() or '.' in sp or any(z in sp for z in ['uncultured','unidentified']):continue
    feats=[x for x in r.features if x.type=='rRNA' and gene.lower() in ' '.join(x.qualifiers.get('product',[])).lower()]
    if not feats:continue
    seq=str(max(feats,key=lambda x:len(x.location)).extract(r.seq)).upper()
    if not 150<=len(seq)<=5500 or sum(c not in 'ACGT' for c in seq)>len(seq)*.03:continue
    fam=next((z for z in reversed(line) if z.endswith('idae')),'unknown')
    item=(sum(c in 'ACGT' for c in seq),r.id,seq,fam)
    if sp not in best or item>best[sp]:best[sp]=item
  (d/(role+'.fa')).write_text(''.join(f'>{a}\n{s}\n' for _,a,s,_ in best.values()))
  with (d/(role+'.csv')).open('w') as h:
   w=csv.writer(h);w.writerow(['species','accession','family']);w.writerows((sp,a,f) for sp,(_,a,s,f) in best.items())
