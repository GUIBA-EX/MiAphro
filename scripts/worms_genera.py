"""Enumerate accepted genera recursively, with complete 50-record pagination."""
from pathlib import Path
import json,csv,urllib.request,urllib.parse,time
root=Path(__file__).resolve().parents[1];out=root/'taxonomy';cache=out/'genus_census';cache.mkdir(exist_ok=True)
def get(id,offset):
 p=cache/f'{id}_{offset}.json'
 if p.exists():return json.loads(p.read_text())
 url=f'https://www.marinespecies.org/rest/AphiaChildrenByAphiaID/{id}?marine_only=false&offset={offset}'
 for k in range(5):
  try:
   with urllib.request.urlopen(url,timeout=50) as r:b=r.read().decode()
   v=json.loads(b) if b else [];p.write_text(json.dumps(v,indent=2));return v
  except Exception:
   if k==4:raise
   time.sleep(3)
def walk(id,fam,seen):
 if id in seen:return []
 seen.add(id);ans=[];offset=1
 while True:
  rows=get(id,offset)
  for r in rows:
   if r['status']!='accepted':continue
   if r['rank']=='Genus':ans.append(dict(family=fam,genus=r['scientificname'],aphia_id=r['AphiaID'],isMarine=r['isMarine'],isExtinct=r['isExtinct'],url=r['url']))
   elif r['taxonRankID']<180:ans+=walk(r['AphiaID'],fam,seen)
  if len(rows)<50:break
  offset+=50
 return ans
allrows=[]
for r in json.load((out/'aphroditiformia_children.json').open()):
 if r['rank']=='Family' and r['status']=='accepted':
  rows=walk(r['AphiaID'],r['scientificname'],set());allrows+=rows;print(r['scientificname'],len(rows),flush=True)
assert len({r['aphia_id'] for r in allrows})==len(allrows)
with (out/'worms_accepted_genera.csv').open('w') as h:
 w=csv.DictWriter(h,allrows[0],lineterminator='\n');w.writeheader();w.writerows(allrows)
