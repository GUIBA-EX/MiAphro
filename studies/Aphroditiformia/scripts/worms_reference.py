"""Cache exact WoRMS genus queries and resolve accepted genus records."""
from pathlib import Path
import csv,json,time,urllib.request,urllib.parse
root=Path(__file__).resolve().parents[1];out=root/'taxonomy';out.mkdir(exist_ok=True)
rows=[dict(gene=g,**r) for g in ['16S','28S'] for r in csv.DictReader((root/g/'target.csv').open())]
def get(endpoint):
 p=out/(urllib.parse.quote(endpoint,safe='')+'.json')
 if p.exists():return json.loads(p.read_text())
 for n in range(4):
  try:
   with urllib.request.urlopen('https://www.marinespecies.org/rest/'+endpoint,timeout=40) as h: b=h.read().decode()
   v=json.loads(b) if b else [];p.write_text(json.dumps(v,indent=2));return v
  except Exception:
   if n==3:raise
   time.sleep(2)
lookup={}
for genus in sorted({r['species'].split()[0] for r in rows}):
 hits=get('AphiaRecordsByName/'+genus+'?like=false&marine_only=true');hits=[x for x in hits if x['rank']=='Genus' and x['phylum']=='Annelida' and x['scientificname']==genus]
 if len(hits)!=1:raise ValueError((genus,hits))
 a=hits[0];b=get('AphiaRecordByAphiaID/'+str(a['valid_AphiaID'])) if a['AphiaID']!=a['valid_AphiaID'] else a
 lookup[genus]=b;print(genus,b['family'],flush=True)
for r in rows:
 b=lookup[r['species'].split()[0]];r.update(worms_genus=b['scientificname'],worms_family=b['family'],aphia_id=b['AphiaID'],worms_url=b['url'],accessed='2026-11-05')
with (out/'references.csv').open('w') as h:
 w=csv.DictWriter(h,rows[0]);w.writeheader();w.writerows(rows)
