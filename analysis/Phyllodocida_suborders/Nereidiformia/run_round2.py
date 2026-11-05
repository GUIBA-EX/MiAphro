from pathlib import Path
import concurrent.futures,subprocess,sys
r=Path(__file__).resolve().parent
def run(gene):
 for band in 'SML':
  with (r/f'round2_{gene}_{band}.log').open('w') as h:subprocess.run([sys.executable,str(r/'scripts/design_round2.py'),gene,band],stdout=h,stderr=h,check=True)
  print('Complete',gene,band,flush=True)
with concurrent.futures.ThreadPoolExecutor(max_workers=2) as p:list(p.map(run,['16S','28S']))
subprocess.run([sys.executable,str(r/'scripts/finish_round2.py')],check=True)
