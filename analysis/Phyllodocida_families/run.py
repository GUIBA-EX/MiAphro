from pathlib import Path
import concurrent.futures,subprocess,sys
root=Path(__file__).resolve().parent
families=['Syllidae','Nereididae','Hesionidae']
subprocess.run([sys.executable,str(root/'controls.py')],check=True)
def run(job):
 fam,gene=job;d=root/fam
 if not (d/gene/'aligned.fa').exists():
  with (d/gene/'aligned.fa').open('w') as out,(d/gene/'alignment.log').open('w') as err:
   subprocess.run(['mafft','--quiet','--thread','2','--adjustdirection','--retree','1','--maxiterate','0',str(d/gene/'target.fa')],stdout=out,stderr=err,check=True)
 for band in 'SML':
  with (d/f'{gene}_{band}.log').open('w') as h:subprocess.run([sys.executable,str(d/'scripts/design_round2.py'),gene,band],stdout=h,stderr=h,check=True)
  print('Complete',fam,gene,band,flush=True)
with concurrent.futures.ThreadPoolExecutor(max_workers=2) as pool:list(pool.map(run,[(f,g) for f in families for g in ['16S','28S']]))
for f in families:
 subprocess.run([sys.executable,str(root/f/'scripts/finish_round2.py')],check=True)
 print('Independent check complete',f,flush=True)
