#!/usr/bin/env python3
"""Refresh local package SHA256SUMS.txt after adding analysis notebooks and docs."""
import hashlib
from pathlib import Path

HERE=Path(__file__).resolve().parents[2]
TARGET=HERE/'reproducibility/manifests/SHA256SUMS.txt'
EXCLUDE_DIRS={'.git','.venv','.venv-analysis','__pycache__','.ipynb_checkpoints','analysis_smoke','build'}


def sha(path):
    h=hashlib.sha256()
    with path.open('rb') as f:
        for block in iter(lambda:f.read(1<<20),b''):
            h.update(block)
    return h.hexdigest()


def main():
    lines=[]
    for p in sorted(HERE.rglob('*')):
        if not p.is_file() or p==TARGET or any(part in EXCLUDE_DIRS for part in p.relative_to(HERE).parts):
            continue
        if p.is_symlink():
            raise ValueError(f'Symlink in repository package: {p}')
        if p.suffix.lower() in ('.hdf5','.h5','.nio','.so','.pyc'):
            raise ValueError(f'Generated simulation result mistakenly included: {p}')
        lines.append(f'{sha(p)}  {p.relative_to(HERE).as_posix()}')
    TARGET.parent.mkdir(parents=True,exist_ok=True)
    TARGET.write_text('\n'.join(lines)+'\n',encoding='utf-8')
    print('SHA256SUMS refreshed:',len(lines),'files')
    print(TARGET)

if __name__=='__main__':
    main()
