"""Cross-check both READMEs, frozen panels, figures and exported primer tables."""
import csv
import re
from collections import Counter
from pathlib import Path
from xml.etree import ElementTree as ET
from Bio import SeqIO

root=Path(__file__).resolve().parents[1]
def rows(path):return list(csv.DictReader((root/path).open()))
def fasta(path):return {r.id:str(r.seq) for r in SeqIO.parse(root/path,'fasta')}
docs=[(root/p).read_text() for p in ('README.md','README.en.md')]
for file in root.rglob('*.md'):
    for url in re.findall(r'\]\(([^)]+)\)',file.read_text()):
        if not url.startswith(('https://','http://','#')):assert (file.parent/url).is_file(),(file,url)
current=rows('SML/results.csv');baseline=rows('baseline/results.csv')
assert len(current)==12 and len(baseline)==4
for sub,rs in [('baseline',baseline),('SML',current)]:
    assert fasta(f'{sub}/primers.fasta')=={r['name']+'_'+side:r[side] for r in rs for side in 'FR'}
for g in ('16S','28S'):
    for role in ('target','background'):
        meta=rows(f'{g}/{role}.csv');lib=fasta(f'{g}/{role}.fa')
        assert len(lib)==len(meta)==len({r['species'] for r in meta})
        assert set(lib)=={r['accession'] for r in meta}
    assert not set(fasta(f'{g}/target.fa')) & set(fasta(f'{g}/background.fa'))
species=rows('SML/species.csv');families=rows('SML/families.csv')
for r in current:
    name=r['name'];g,band=name.split('-')[1:]
    path=f'figures/template_matching/{name}.svg'; labels=list(ET.parse(root/path).getroot().itertext())
    assert r['F'] in labels and r['R'] in labels
    for doc in docs:
        line=next(l for l in doc.splitlines() if l.startswith('|') and path in l)
        assert [v.strip() for v in line.split('|')[1:8]]==[band,r['span'],'`'+r['F']+'`','`'+r['R']+'`',r['covered']+'/'+r['target_species'],r['exact']+'/'+r['target_species'],r['WC']+'/'+r['GT']]
    ss=[x for x in species if x['name']==name];assert len(ss)==int(r['target_species'])==len({x['accession'] for x in ss})
    for x in ss:
        idx=labels.index(x['species']+' · '+x['accession'])
        assert labels[idx+1]==('Exact' if x['exact']=='1' else 'Hit' if x['covered']=='1' else 'Miss')
    ff=[x for x in families if x['name']==name]
    for k in ('covered','exact'):
        assert sum(int(x[k]) for x in ss)==sum(int(x[k]) for x in ff)==int(r[k])
    pos=list(ET.parse(root/f'figures/{g}_positions.svg').getroot().itertext())
    for side in 'FR':assert f"{side} {r[side+'_start']}–{r[side+'_end']}" in pos
print('PASS: Markdown links, both result tables, all current/baseline FASTA sequences, panel membership, binding coordinates, species-level figure outcomes and family totals.')
