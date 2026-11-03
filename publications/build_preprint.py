"""Render the manuscript using pandoc and Chrome/Chromium (set CHROME if needed)."""
import os
from pathlib import Path
import re
import shutil
import subprocess
import tempfile

out=Path(__file__).resolve().parent
article=subprocess.check_output(['pandoc','--from=gfm','--to=html5',str(out/'MiAphro_preprint.md')],text=True)
article=re.sub(r'href="(?!https?:|#)([^"]+)"',lambda m:'href="https://github.com/GUIBA-EX/MiAphro/blob/main/'+str((out/m[1]).resolve().relative_to(out.parent))+'"',article)
css="""@page{size:A4;margin:17mm 18mm 18mm;@bottom-left{content:'MiAphro · Fei Xia · Computational preprint';font:8pt Arial;color:#54635c}@bottom-right{content:counter(page);font:8pt Arial}}body{font:11pt/1.46 Georgia,serif;color:#172b27}h1{font-size:21pt;text-align:center;line-height:1.18}h2{font-size:15pt}h3{font-size:12pt}h2,h3{break-after:avoid}p{margin:6pt 0 9pt;orphans:3;widows:3}img{width:100%;display:block}img[src$="_positions.svg"]{width:86%;margin:auto}figure{margin:8pt 0;break-inside:avoid;break-after:avoid}figcaption{display:none}table{border-collapse:collapse;width:100%;font:8.7pt/1.4 Arial}td,th{padding:3pt;border-bottom:1px solid #cfdfda}th{background:#deede8}tr{break-inside:avoid}a{color:#216459}h2#materials-and-methods,h2[id^=results],h2#discussion,h2#primer-sequences,h2#primertemplate-comparisons{break-before:page}p:has(>strong:first-child){break-before:avoid}"""
chrome=os.environ.get('CHROME') or shutil.which('chromium') or '/Applications/Google Chrome.app/Contents/MacOS/Google Chrome'
with tempfile.TemporaryDirectory() as t:
    doc=Path(t)/'paper.html'
    doc.write_text(f'<!doctype html><meta charset=utf-8><base href="{out.as_uri()}/"><title>MiAphro · Fei Xia</title><style>{css}</style>{article}')
    subprocess.run([chrome,'--headless','--disable-gpu','--allow-file-access-from-files','--no-pdf-header-footer',f'--user-data-dir={t}/profile',f'--print-to-pdf={out}/MiAphro_preprint.pdf',doc.as_uri()],check=True)
