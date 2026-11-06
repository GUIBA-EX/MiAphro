"""Recompute primary-pair structures across all concrete IUPAC variants."""
import csv
import itertools
import json
import sys
from pathlib import Path
import primer3

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT/'optimization/terminal_2026-11-05'))
from chemistry_rules import expand, single, hetero


def main():
    records = list(csv.DictReader((ROOT/'catalog.csv').open()))
    result = []
    for row in records:
        forward, reverse = expand(row['F']), expand(row['R'])
        for label in ('50_per_variant', '200_per_direction'):
            fc = 50 if label == '50_per_variant' else 200/len(forward)
            rc = 50 if label == '50_per_variant' else 200/len(reverse)
            tm = [(single(s, c)[0], end, s) for end, variants, c in [('F',forward,fc),('R',reverse,rc)] for s in variants]
            cases = []
            for end, variants, c in [('F',forward,fc),('R',reverse,rc)]:
                for seq in variants:
                    kw = dict(mv_conc=50,dv_conc=1.5,dntp_conc=.6,dna_conc=c)
                    cases += [(primer3.calc_hairpin(seq,**kw).tm, end+'_hairpin',seq,''), (primer3.calc_homodimer(seq,**kw).tm,end+'_homodimer',seq,seq)]
                cases += [(hetero(a,b,c),end+'_variant_heterodimer',a,b) for a,b in itertools.combinations(variants,2)]
            cases += [(hetero(a,b,max(fc,rc)),'F_R_heterodimer',a,b) for a,b in itertools.product(forward,reverse)]
            worst = max(cases); lo=min(x[0] for x in tm); hi=max(x[0] for x in tm)
            fts=[x[0] for x in tm if x[1]=='F']; rts=[x[0] for x in tm if x[1]=='R']
            delta=abs((min(fts)+max(fts)-min(rts)-max(rts))/2)
            result.append(dict(name=row['name'],scenario=label,F_nM=fc,R_nM=rc,primer_min_tm=lo,primer_max_tm=hi,structure_tm=worst[0],structure_type=worst[1],variant_1=worst[2],variant_2=worst[3],margin=lo-worst[0],midpoint_difference=delta,basic_pass=int(worst[0]<55 and lo>=45 and hi<=75 and delta<=5)))
    out=ROOT/'audit/review_20261107';out.mkdir(parents=True,exist_ok=True)
    with (out/'all_pair_structures.csv').open('w') as f:
        writer=csv.DictWriter(f,fieldnames=list(result[0]),lineterminator='\n');writer.writeheader();writer.writerows(result)
    (out/'structure_metadata.json').write_text(json.dumps({'primer3_py':primer3.__version__,'pairs':len(records),'scenarios':2,'model':'50 mM monovalent, 1.5 mM Mg, 0.6 mM dNTP; no adapters; F/R heterodimers use higher participating concentration','status':'PASS' if all(x['basic_pass'] for x in result) else 'BASIC_GATE_FAILURE'},indent=2)+'\n')
    print('pairs',len(records),'low-margin',len({r['name'] for r in result if r['margin']<5}),'basic failures',sum(not r['basic_pass'] for r in result))

if __name__=='__main__':
    main()
