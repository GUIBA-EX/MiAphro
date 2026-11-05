from pathlib import Path
import hashlib
R=Path(__file__).resolve().parents[1]
files=[p for p in R.rglob('*') if p.is_file() and not any(x in {'.git','.recheck','__pycache__'} for x in p.relative_to(R).parts) and p.name not in {'SHA256SUMS','.DS_Store','scan','target_barcodes.tsv'} and p.suffix!='.pyc']
(R/'SHA256SUMS').write_text(''.join(hashlib.sha256(p.read_bytes()).hexdigest()+'  '+str(p.relative_to(R))+'\n' for p in sorted(files)))
print('Checksums',len(files),'files')
