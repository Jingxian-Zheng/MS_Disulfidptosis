"""Offline checks of the complete upload tree, with no optional imports."""
from __future__ import annotations
import ast
import csv
import hashlib
import json
from pathlib import Path
import re
import sys
import tomllib

ROOT=Path(__file__).resolve().parents[1]

def check(root=ROOT):
    errors=[]
    for path in root.rglob('*'):
        if not path.is_file() or any(p in {'.git','work','__pycache__'} for p in path.relative_to(root).parts):
            continue
        rel=path.relative_to(root).as_posix()
        if path.stat().st_size >=25*1024*1024:
            errors.append(rel+': exceeds browser file limit')
        if path.suffix=='.py':
            try:ast.parse(path.read_text(encoding='utf-8-sig'),filename=rel)
            except SyntaxError as exc:errors.append(f'{rel}:{exc.lineno}: {exc.msg}')
        if path.suffix=='.toml':
            try:tomllib.loads(path.read_text(encoding='utf-8-sig'))
            except tomllib.TOMLDecodeError as exc:errors.append(rel+': '+str(exc))
        if path.suffix in {'.py','.R','.r','.jsx','.ps1','.md','.json','.toml','.yml','.yaml','.csv','.tsv','.txt'}:
            text=path.read_text(encoding='utf-8-sig',errors='replace')
            patterns=[r'[A-Za-z]:[\\/]Users[\\/][^\s\"\']+',r'/data/' + r'fulei/',
                      r'gh[pousr]_[A-Za-z0-9]{30,}',r'github_pat_[A-Za-z0-9_]{30,}',
                      r'AKIA[A-Z0-9]{16}',r'-----BEGIN (?:RSA |EC |OPENSSH )?PRIVATE KEY-----']
            if any(re.search(pattern,text) for pattern in patterns):
                errors.append(rel+': workstation path or credential-like content')
    manifest=json.loads((root/'provenance/release_files.json').read_text(encoding='utf-8'))
    for row in manifest:
        f=root/row['path']
        if not f.is_file() or hashlib.sha256(f.read_bytes()).hexdigest()!=row['sha256']:
            errors.append(row['path']+': release checksum mismatch')
    sources=list(csv.DictReader((root/'provenance/source_inventory.csv').open(encoding='utf-8-sig',newline='')))
    for row in sources:
        if not (root/row['release_path']).is_file():errors.append('Unmapped source: '+row['source_path'])
    status={'status':'passed' if not errors else 'failed','source_files_accounted_for':len(sources),
            'public_source_records':sum(r.get('distribution')=='public_code' for r in sources),
            'local_manuscript_source_records':sum(r.get('distribution')=='local_manuscript_tools' for r in sources),
            'release_manifest_files':len(manifest),'errors':errors}
    print(json.dumps(status,ensure_ascii=False,indent=2))
    return 0 if not errors else 1

if __name__=='__main__':
    raise SystemExit(check())
