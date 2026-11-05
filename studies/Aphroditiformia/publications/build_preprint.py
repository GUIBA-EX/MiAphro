"""Render the manuscript using pandoc and Chrome/Chromium (set CHROME if needed)."""
import os,json
from pathlib import Path
import re
import shutil
import subprocess
import tempfile

out=Path(__file__).resolve().parent
c={'title': 'MiAphro: in silico design and evaluation of mitochondrial and nuclear rRNA primers for scale-worm (Aphroditiformia) eDNA metabarcoding', 'version': '1.0'}
article=subprocess.check_output(['pandoc','--from=gfm','--to=html5',str(out/'MiAphro_preprint.md')],text=True)
article=re.sub(r'href="(?!https?:|#)([^"]+)"',lambda m:'href="https://github.com/GUIBA-EX/MiAphro/blob/main/'+str((out/m[1]).resolve().relative_to(out.parent))+'"',article)
article=re.sub(r'(<p><img\b[^>]*></p>|<figure>.*?</figure>)\s*(<p><strong>Figure.*?</p>)',r'<div class=figure-block>\1\2</div>',article,flags=re.S)
article=article.replace('′','&#39;')
css="""@page{size:A4;margin:17mm 18mm 18mm;@bottom-left{content:'MiAphro · Fei Xia · Computational preprint';font:8pt 'Avenir Next';color:#54635c}@bottom-right{content:counter(page);font:8pt 'Avenir Next'}}.figure-block{break-inside:avoid}body{font:10.5pt/1.38 'Avenir Next',sans-serif;color:#172b27}h1{font-size:21pt;text-align:center;line-height:1.18}h2{font-size:15pt}h3{font-size:12pt}h2,h3{break-after:avoid}p{margin:6pt 0 9pt;orphans:3;widows:3}img{width:100%;display:block}img[src$="_positions.svg"]{width:86%;margin:auto}figure{margin:8pt 0;break-inside:avoid;break-after:avoid}figcaption{display:none}table{break-inside:avoid;border-collapse:collapse;width:100%;font:8.7pt/1.4 'Avenir Next'}td,th{padding:3pt;border-bottom:1px solid #cfdfda}th{background:#deede8}tr{break-inside:avoid}a{color:#216459}h2[id^=table]{break-before:page}code{font-family:'Avenir Next',sans-serif}p:has(>strong:first-child){break-before:avoid}h2#references~p{font-size:9pt;line-height:1.25;margin:5pt 0}h2[id^="table-1"]+p+table{font-size:8pt;line-height:1.1}h2[id^="table-1"]+p+table td,h2[id^="table-1"]+p+table th{padding:1.5pt}"""
chrome=os.environ.get('CHROME') or shutil.which('chromium') or '/Applications/Google Chrome.app/Contents/MacOS/Google Chrome'
with tempfile.TemporaryDirectory() as t:
    doc=Path(t)/'paper.html'
    doc.write_text(f'<!doctype html><meta charset=utf-8><base href="{out.as_uri()}/"><title>MiAphro · Fei Xia</title><style>{css}</style>{article}')
    subprocess.run([chrome,'--headless','--disable-gpu','--allow-file-access-from-files','--no-pdf-header-footer',f'--user-data-dir={t}/profile',f'--print-to-pdf={out}/MiAphro_preprint.pdf',doc.as_uri()],check=True)

