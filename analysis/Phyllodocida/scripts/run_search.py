from pathlib import Path
import subprocess,sys,concurrent.futures
root=Path(__file__).resolve().parents[1]
for gene in ['16S','28S']:
 with (root/gene/'aligned.fa').open('w') as h, (root/gene/'mafft.log').open('w') as e:
  subprocess.run(['mafft','--thread','2','--adjustdirection','--retree','1','--maxiterate','0',str(root/gene/'target.fa')],stdout=h,stderr=e,check=True)
def run(item):
 gene,band=item;log=root/f'{gene}_{band}.log'
 with log.open('w') as h:subprocess.run([sys.executable,str(root/'scripts/design.py'),gene,band],stdout=h,stderr=h,check=True)
 print('Finished',gene,band,flush=True)
with concurrent.futures.ThreadPoolExecutor(max_workers=2) as pool:list(pool.map(run,[(g,b) for g in ['16S','28S'] for b in 'SML']))
subprocess.run([sys.executable,str(root/'scripts/finish_sml.py')],check=True)
