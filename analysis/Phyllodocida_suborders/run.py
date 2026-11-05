from pathlib import Path
import subprocess,concurrent.futures,sys
r=Path(__file__).resolve().parent
names=['Glyceriformia','Nereidiformia','Phyllodociformia']
def run(job):
 name,g,b=job;d=r/name
 with (d/f'{g}_{b}.log').open('w') as h:subprocess.run([sys.executable,str(d/'scripts/design_expanded.py'),g,b],stdout=h,stderr=h,check=True)
 print('Complete',name,g,b,flush=True)
with concurrent.futures.ThreadPoolExecutor(max_workers=2) as pool:list(pool.map(run,[(n,g,b) for g in ['16S','28S'] for b in 'SML' for n in names]))
for n in names:subprocess.run([sys.executable,str(r/n/'scripts/finish_expanded.py')],check=True)
