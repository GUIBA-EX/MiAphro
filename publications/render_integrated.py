"""Build the integrated manuscript and print its HTML to PDF with Chrome."""
from pathlib import Path
import subprocess,sys,os,tempfile
P=Path(__file__).resolve().parent
subprocess.run([sys.executable,str(P/'build_integrated.py')],check=True)
chrome=os.environ.get('CHROME','/Applications/Google Chrome.app/Contents/MacOS/Google Chrome')
with tempfile.TemporaryDirectory() as t:
 subprocess.run([chrome,'--headless','--disable-gpu','--allow-file-access-from-files','--no-pdf-header-footer','--user-data-dir='+t,'--print-to-pdf='+str(P/'MiPhyllo_preprint.pdf'),(P/'MiPhyllo_preprint.html').as_uri()],check=True)
