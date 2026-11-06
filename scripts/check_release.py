"""Structural release gate for the integrated MiPhyllo catalogue."""
from pathlib import Path
import csv, json, re, math, hashlib, sys
from collections import Counter
from Bio import SeqIO
from Bio.Data.IUPACData import ambiguous_dna_values as I
R = Path(__file__).resolve().parents[1]
read = lambda p: list(csv.DictReader(p.open()))
rows = read(R / 'catalog.csv')
alt = read(R / 'conditional/results.csv')
if not (len(rows) == 48 and len(alt) == 0):
    raise ValueError('Validation failed: len(rows) == 48 and len(alt) == 0')
if not len({x['name'] for x in rows + alt}) == 48:
    raise ValueError("Validation failed: len({x['name'] for x in rows + alt}) == 48")
verified = json.loads((R / 'audit/catalog_rescan.json').read_text())
if not (verified['status'] == 'PASS' and verified['pairs'] == 102):
    raise ValueError("Validation failed: verified['status'] == 'PASS' and verified['pairs'] == 102")
if not {(x['name'], x['profile']) for x in rows + alt} <= {(x['name'], x['profile']) for x in verified['records']}:
    raise ValueError("Validation failed: {(x['name'], x['profile']) for x in rows + alt} <= {(x['name'], x['profile']) for x in verified['records']}")
for rs, folder in [(rows, R), (alt, R / 'conditional')]:
    seq = {x.id: str(x.seq) for x in SeqIO.parse(folder / 'primers.fasta', 'fasta')}
    if not seq == {x['name'] + '_' + side: x[side] for x in rs for side in 'FR'}:
        raise ValueError("Validation failed: seq == {x['name'] + '_' + side: x[side] for x in rs for side in 'FR'}")
    details = read(folder / 'species.csv')
    lookup = {x['name']: x for x in rs}
    for x in rs:
        expected = f'MiPhyllo-{x['target_taxon']}-{x['gene']}-{x['band']}' + ('2' if x['profile'] == 'conditional' else '') + ('D' if x['variant'] == 'D' else '')
        if not x['name'] == expected:
            raise ValueError("Validation failed: x['name'] == expected")
        for side in 'FR':
            if not len(x[side]) == int(x[side + '_end']) - int(x[side + '_start']) + 1:
                raise ValueError("Validation failed: len(x[side]) == int(x[side + '_end']) - int(x[side + '_start']) + 1")
            if not math.prod((len(I[b]) for b in x[side])) == int(x['deg_' + side]) <= 8:
                raise ValueError("Validation failed: math.prod((len(I[b]) for b in x[side])) == int(x['deg_' + side]) <= 8")
            if not set(x[side][:len(x[side]) // 2]) <= set('ACGT'):
                raise ValueError("Validation failed: set(x[side][:len(x[side]) // 2]) <= set('ACGT')")
            if x['variant'] == 'ordinary':
                if not set(x[side]) <= set('ACGT'):
                    raise ValueError("Validation failed: set(x[side]) <= set('ACGT')")
        if x['variant'] == 'D':
            a = lookup[x['name'][:-1]]
            if not all((a[k] == x[k] for k in ['reference', 'F_start', 'F_end', 'R_start', 'R_end'])):
                raise ValueError("Validation failed: all((a[k] == x[k] for k in ['reference', 'F_start', 'F_end', 'R_start', 'R_end']))")
        ss = [s for s in details if s['name'] == x['name']]
        if not len(ss) == int(x['target_species']) == len({s['species'] for s in ss}):
            raise ValueError("Validation failed: len(ss) == int(x['target_species']) == len({s['species'] for s in ss})")
        for k in ['covered', 'exact']:
            if not sum((int(s[k]) for s in ss)) == int(x[k]):
                raise ValueError('Validation failed: sum((int(s[k]) for s in ss)) == int(x[k])')
        docroot = R / 'conditional' if x['profile'] == 'conditional' else R / 'panels' / x['target_taxon']
        for lang in ['README.md', 'README.en.md']:
            line = next((t for t in (docroot / lang).read_text().splitlines() if t.startswith('|' + x['name'] + '|') and '`' in t))
            if not (all(('`' + x[k] + '`' in line for k in ['F', 'R'])) and f'{x['covered']}/{x['target_species']}' in line and (f'{x['exact']}/{x['target_species']}' in line)):
                raise ValueError("Validation failed: all(('`' + x[k] + '`' in line for k in ['F', 'R'])) and f'{x['covered']}/{x['target_species']}' in line and (f'{x['exact']}/{x['target_species']}' in line)")
for x in read(R / 'panels.csv'):
    p = R / 'panels' / x['target_taxon'] / 'references' / x['gene']
    species = {}
    for role in ['target', 'background']:
        metadata = read(p / (role + '.csv'))
        seq = list(SeqIO.parse(p / (role + '.fa'), 'fasta'))
        if not len(seq) == len(metadata) == len({z['species'] for z in metadata}) == int(x[role + '_species']):
            raise ValueError("Validation failed: len(seq) == len(metadata) == len({z['species'] for z in metadata}) == int(x[role + '_species'])")
        if not {z.id for z in seq} == {z['accession'] for z in metadata}:
            raise ValueError("Validation failed: {z.id for z in seq} == {z['accession'] for z in metadata}")
        species[role] = {z['species'] for z in metadata}
    if not not species['target'] & species['background']:
        raise ValueError("Validation failed: not species['target'] & species['background']")
for p in R.rglob('*.md'):
    if '.git' in p.parts or '.recheck' in p.parts:
        continue
    for link in re.findall('\\]\\(([^)]+)\\)', p.read_text()):
        if not link.startswith(('https://', 'http://', '#')):
            if not (p.parent / link.split('#')[0]).is_file():
                raise ValueError("Validation failed: (p.parent / link.split('#')[0]).is_file()", (p, link))
for x in json.loads((R / 'analysis/source_manifest.json').read_text()):
    if not hashlib.sha256((R / x['path']).read_bytes()).hexdigest() == x['sha256']:
        raise ValueError("Validation failed: hashlib.sha256((R / x['path']).read_bytes()).hexdigest() == x['sha256']", x['path'])
if not (R / 'figures/catalog_coverage.svg').exists():
    raise ValueError("Validation failed: (R / 'figures/catalog_coverage.svg').exists()")
if (R / 'SHA256SUMS').exists():
    for line in (R / 'SHA256SUMS').read_text().splitlines():
        digest, path = line.split('  ', 1)
        if not hashlib.sha256((R / path).read_bytes()).hexdigest() == digest:
            raise ValueError('Validation failed: hashlib.sha256((R / path).read_bytes()).hexdigest() == digest', path)
print('PASS: 48 current pairs; 54 archived pairs; unique names; FASTA; same endpoints; degeneracy; panels; bilingual tables; links; source hashes; independent rescans.')
original = read(R / 'archive/2026-11-05-unfiltered/catalog.csv')
original_alt = read(R / 'archive/2026-11-05-unfiltered/conditional/results.csv')
expected = set()
for x in original:
    if x['variant'] == 'D' and x['target_taxon'] not in {'Phyllodocida', 'Nereidiformia'} and (2 * int(x['covered']) >= int(x['target_species'])) and (int(x['exact']) > 0):
        expected.update((x['name'], x['name'][:-1]))
if not {x['name'] for x in rows} == expected:
    raise ValueError("Validation failed: {x['name'] for x in rows} == expected")
updates = {x['name']: x for x in read(R / 'optimization/terminal_2026-11-05/selected_results.csv')}
if not len(updates) == 4:
    raise ValueError('Validation failed: len(updates) == 4')
if not rows == [updates.get(x['name'], x) for x in original if x['name'] in expected]:
    raise ValueError("Validation failed: rows == [updates.get(x['name'], x) for x in original if x['name'] in expected]")
comparisons = json.loads((R / 'optimization/terminal_2026-11-05/comparison.json').read_text())
if not (len(comparisons) == 9 and all((x['independent_status'] == 'PASS' for x in comparisons))):
    raise ValueError("Validation failed: len(comparisons) == 9 and all((x['independent_status'] == 'PASS' for x in comparisons))")
for x in comparisons:
    if x['status'] == 'replace':
        if not (x['after']['thermal_pass'] and x['after']['background_pass']):
            raise ValueError("Validation failed: x['after']['thermal_pass'] and x['after']['background_pass']")
        if not (x['after']['covered'] >= x['before']['covered'] and x['after']['exact'] >= x['before']['exact']):
            raise ValueError("Validation failed: x['after']['covered'] >= x['before']['covered'] and x['after']['exact'] >= x['before']['exact']")
        if not any((x['after'][k] > x['before'][k] for k in ['covered', 'exact'])):
            raise ValueError("Validation failed: any((x['after'][k] > x['before'][k] for k in ['covered', 'exact']))")
        if not all((str(x['after'][k]) == updates[x['name']][k] for k in ['F', 'R', 'covered', 'exact', 'WC', 'GT', 'WCloose', 'GTloose'])):
            raise ValueError("Validation failed: all((str(x['after'][k]) == updates[x['name']][k] for k in ['F', 'R', 'covered', 'exact', 'WC', 'GT', 'WCloose', 'GTloose']))")
removed = read(R / 'archive/removed_candidates.csv')
if not {x['name'] for x in removed} == {x['name'] for x in original + original_alt} - expected:
    raise ValueError("Validation failed: {x['name'] for x in removed} == {x['name'] for x in original + original_alt} - expected")
if not len(removed) == 54:
    raise ValueError('Validation failed: len(removed) == 54')
for p in (R / 'panels').iterdir():
    if not p.is_dir():
        continue
    expected_rows = [x for x in rows if x['target_taxon'] == p.name]
    if not read(p / 'results.csv') == expected_rows:
        raise ValueError("Validation failed: read(p / 'results.csv') == expected_rows")
    if not {x.id: str(x.seq) for x in SeqIO.parse(p / 'primers.fasta', 'fasta')} == {x['name'] + '_' + e: x[e] for x in expected_rows for e in 'FR'}:
        raise ValueError("Validation failed: {x.id: str(x.seq) for x in SeqIO.parse(p / 'primers.fasta', 'fasta')} == {x['name'] + '_' + e: x[e] for x in expected_rows for e in 'FR'}")
manifest = json.loads((R / 'publications/MiPhyllo_manifest.json').read_text())
if not (manifest['current_pairs'] == 48 and manifest['conditional_pairs'] == 0 and (len(manifest['units']) == 6)):
    raise ValueError("Validation failed: manifest['current_pairs'] == 48 and manifest['conditional_pairs'] == 0 and (len(manifest['units']) == 6)")
qa = json.loads((R / 'audit/integrated_preprint.json').read_text())
if not (qa['status'] == 'PASS' and qa['pairs'] == 48 and (qa['oligos_checked'] == 96)):
    raise ValueError("Validation failed: qa['status'] == 'PASS' and qa['pairs'] == 48 and (qa['oligos_checked'] == 96)")
if not qa['pdf_sha256'] == hashlib.sha256((R / 'publications/MiPhyllo_preprint.pdf').read_bytes()).hexdigest():
    raise ValueError("Validation failed: qa['pdf_sha256'] == hashlib.sha256((R / 'publications/MiPhyllo_preprint.pdf').read_bytes()).hexdigest()")
print('PASS: editorial criteria, archive partition, panel exports and preprint identity.')
