"""Re-scan every current and conditional pair against its own frozen panel."""
from pathlib import Path
import csv, json, subprocess, shutil, concurrent.futures
from collections import Counter
from Bio import SeqIO
from independent import flags
R = Path(__file__).resolve().parents[1]
read = lambda p: list(csv.DictReader(p.open()))

def run(job):
    taxon, gene = job
    rows = [x for x in read(R / 'catalog.csv') + read(R / 'conditional/results.csv') if x['target_taxon'] == taxon and x['gene'] == gene]
    d = R / '.recheck' / taxon / gene
    d.mkdir(parents=True, exist_ok=True)
    source = R / 'panels' / taxon / 'references' / gene
    for role in ['target', 'background']:
        shutil.copyfile(source / (role + '.fa'), d / (role + '.fasta'))
    ids = {f'p{i}': x for i, x in enumerate(rows)}
    (d / 'independent_primers.tsv').write_text(''.join((f'{k}\t{x['F']}\t{x['R']}\n' for k, x in ids.items())))
    subprocess.run([str(R / 'scripts/scan'), str(d)], stderr=subprocess.DEVNULL, check=True)
    hits = {(x['name'], x['role'], x['model'], x['accession']): x for x in csv.DictReader((d / 'independent_hits.tsv').open(), delimiter='\t')}
    libs = {role: {x.id: str(x.seq).upper() for x in SeqIO.parse(source / (role + '.fa'), 'fasta')} for role in ['target', 'background']}
    answer = []
    for key, x in ids.items():
        count = Counter()
        for role, seqs in libs.items():
            for acc, s in seqs.items():
                for gt in [False] if role == 'target' else [False, True]:
                    model = 'GT' if gt else 'WC'
                    e, c, l = flags(s, x['F'], x['R'], role == 'target', gt)
                    y = hits.get((key, role, model, acc), {})
                    if not [e, c, l] == [y.get(k) == '1' for k in ['exact', 'le2last3', 'le3' if role == 'background' and (not gt) else 'le2']]:
                        raise ValueError("Validation failed: [e, c, l] == [y.get(k) == '1' for k in ['exact', 'le2last3', 'le3' if role == 'background' and (not gt) else 'le2']]", (x['name'], role, acc))
                    if role == 'target':
                        count['covered'] += c
                        count['exact'] += e
                    else:
                        count[model] += c
                        count[model + 'loose'] += l
        if not all((count[k] == int(x[k]) for k in ['covered', 'exact', 'WC', 'GT', 'WCloose', 'GTloose'])):
            raise ValueError("Validation failed: all((count[k] == int(x[k]) for k in ['covered', 'exact', 'WC', 'GT', 'WCloose', 'GTloose']))", (x['name'], count))
        answer.append(dict(name=x['name'], profile=x['profile'], status='PASS', target_species=len(libs['target']), background_species=len(libs['background'])))
    print('PASS', taxon, gene, len(answer), 'pairs', flush=True)
    return answer
if __name__ == '__main__':
    subprocess.run(['c++', '-std=c++17', '-O2', str(R / 'scripts/raw_scan.cpp'), '-o', str(R / 'scripts/scan')], check=True)
    jobs = list(dict.fromkeys(((x['target_taxon'], x['gene']) for x in read(R / 'catalog.csv'))))
    with concurrent.futures.ProcessPoolExecutor(max_workers=2) as p:
        result = [x for batch in p.map(run, jobs) for x in batch]
    (R / 'audit').mkdir(exist_ok=True)
    (R / 'audit/current_rescan.json').write_text(json.dumps(dict(status='PASS', pairs=len(result), records=result), indent=2))
    print('ALL PASS', len(result), flush=True)
