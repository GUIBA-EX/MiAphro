"""Summarize exact insert identity within same-assay recovered references.

Run after recheck_catalog.py; this is not full-library resolution.
"""
import collections
import csv
from pathlib import Path

ROOT=Path(__file__).resolve().parents[1]
read=lambda path:list(csv.DictReader(path.open()))

def main():
    rows=read(ROOT/'catalog.csv'); summary=[]; details=[]
    for taxon,gene in dict.fromkeys((r['target_taxon'],r['gene']) for r in rows):
        pairs=[r for r in rows if r['target_taxon']==taxon and r['gene']==gene]
        meta={r['accession']:r['species'] for r in read(ROOT/'panels'/taxon/'references'/gene/'target.csv')}
        found=collections.defaultdict(set)
        path=ROOT/'.recheck'/taxon/gene/'target_barcodes.tsv'
        with path.open() as handle:
            for hit in csv.DictReader(handle,delimiter='\t'):
                found[(hit['name'],hit['accession'])].add(hit['sequence'])
        for i,row in enumerate(pairs):
            key=f'p{i}'
            records={acc:next(iter(seq)) for (name,acc),seq in found.items() if name==key and len(seq)==1}
            multi=sum(name==key and len(seq)>1 for (name,acc),seq in found.items())
            names=collections.defaultdict(set)
            for acc,seq in records.items():names[seq].add(meta[acc])
            unique=sum(len(names[seq])==1 for seq in records.values())
            summary.append(dict(name=row['name'],target_records=len(meta),pair_covered=row['covered'],single_clean_insert_records=len(records),multiple_clean_insert_records=multi,identity_unique_among_recovered=unique,scope='unique against other recovered species labels only'))
            for acc,seq in records.items():details.append(dict(name=row['name'],accession=acc,species=meta[acc],insert=seq,unique_among_recovered=int(len(names[seq])==1)))
    out=ROOT/'audit/review_20261107'
    for filename,data in [('recovered_insert_identity.csv',summary),('recovered_inserts.csv',details)]:
        with (out/filename).open('w') as h:
            w=csv.DictWriter(h,fieldnames=list(data[0]),lineterminator='\n');w.writeheader();w.writerows(data)
    print('pairs',len(summary),'single-insert records',len(details),'multiple',sum(r['multiple_clean_insert_records'] for r in summary))

if __name__=='__main__':main()
