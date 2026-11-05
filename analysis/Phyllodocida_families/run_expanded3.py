from pathlib import Path
import subprocess,concurrent.futures,sys
R=Path(__file__).resolve().parent
def run(job):
 f,g=job
 for b in 'SML':
  with (R/f/f'expanded3_{g}_{b}.log').open('w') as h:subprocess.run([sys.executable,str(R/f/'scripts/design_expanded3.py'),g,b],stdout=h,stderr=h,check=True)
  print('Complete',f,g,b,flush=True)
with concurrent.futures.ThreadPoolExecutor(max_workers=2) as p:list(p.map(run,[('Syllidae','16S'),('Nereididae','28S')]))
