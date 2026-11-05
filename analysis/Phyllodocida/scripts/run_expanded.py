from pathlib import Path
import subprocess,sys,concurrent.futures
r=Path(__file__).resolve().parents[1]
def run(gb):
 g,b=gb
 with (r/f'expanded_{g}_{b}.log').open('w') as h:subprocess.run([sys.executable,str(r/'scripts/design_expanded.py'),g,b],stdout=h,stderr=h,check=True)
 print('Complete',g,b,flush=True)
with concurrent.futures.ThreadPoolExecutor(max_workers=2) as p:list(p.map(run,[(g,b) for g in ['16S','28S'] for b in 'SML']))
subprocess.run([sys.executable,str(r/'scripts/finish_expanded.py')],check=True)
