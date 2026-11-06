"""Audit homologous projected 3-prime 1/3/5 windows on frozen alignments."""
import csv
from pathlib import Path
from collections import Counter
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.Data.IUPACData import ambiguous_dna_values as I
ROOT=Path(__file__).resolve().parents[1]

def read(p):
    with p.open() as f:return list(csv.DictReader(f))

def main():
    details=[]; summaries=[]; cache={}
    for row in read(ROOT/'catalog.csv'):
        key=row['target_taxon'],row['gene'],row['reference']
        if key not in cache:
            directory=ROOT/'panels'/key[0]/'references'/key[1]
            alignment={r.id:str(r.seq).upper() for r in SeqIO.parse(directory/'aligned.fa','fasta')}
            if key[2] not in alignment:alignment={r.id:str(r.seq).upper() for r in SeqIO.parse(directory/'coordinate_reference_alignment.fa','fasta')}
            if key[2] not in alignment:raise ValueError(f'Missing anchor {key}')
            cache[key]=(alignment,read(directory/'target.csv'),[i for i,b in enumerate(alignment[key[2]]) if b!='-'])
        alignment,meta,anchor=cache[key];counts=Counter()
        for rec in meta:
            seq=alignment[rec['accession']]; item=dict(name=row['name'],accession=rec['accession'],species=rec['species'])
            for end in ('F','R'):
                cols=[anchor[i-1] for i in range(int(row[end+'_start']),int(row[end+'_end'])+1)]
                window=seq[min(cols):max(cols)+1].replace('-','')
                if len(window)!=len(row[end]):
                    status={n:'indel_or_incomplete' for n in (1,3,5)}
                else:
                    if end=='R':window=str(Seq(window).reverse_complement())
                    status={}
                    for n in (1,3,5):
                        frag=window[-n:];primer=row[end][-n:]
                        status[n]='ambiguous' if set(frag)-set('ACGT') else ('match' if all(base in I[code] for base,code in zip(frag,primer)) else 'mismatch')
                for n,v in status.items():item[f'{end}_last{n}']=v;counts[(end,n,v)]+=1
            for n in (1,3,5):
                pair=(item[f'F_last{n}'],item[f'R_last{n}'])
                item[f'paired_last{n}']='match' if pair==('match','match') else ('unresolved' if any(x in ('ambiguous','indel_or_incomplete') for x in pair) else 'mismatch')
                counts[('pair',n,item[f'paired_last{n}'])]+=1
            details.append(item)
        result=dict(name=row['name'],target=row['target_taxon'],gene=row['gene'],variant=row['variant'],reference_records=len(meta))
        for end in ('F','R','pair'):
            for n in (1,3,5):
                for status in (('match','mismatch','unresolved') if end=='pair' else ('match','mismatch','ambiguous','indel_or_incomplete')):
                    result[f'{end}_last{n}_{status}']=counts[(end,n,status)]
        summaries.append(result)
    out=ROOT/'audit/review_20261107'
    for name,rows in [('terminal_windows_records.csv',details),('terminal_windows_summary.csv',summaries)]:
        with (out/name).open('w') as f:
            w=csv.DictWriter(f,fieldnames=list(rows[0]),lineterminator='\n');w.writeheader();w.writerows(rows)
    print('pairs',len(summaries),'records',len(details))
if __name__=='__main__':main()
