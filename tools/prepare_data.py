"""Verify the study data archive before extracting into a new workspace."""
from __future__ import annotations
import argparse
import hashlib
import json
from pathlib import Path, PurePosixPath
import shutil
import stat
import tempfile
import zipfile

ROOT = Path(__file__).resolve().parents[1]

def sha256(path):
    h = hashlib.sha256()
    with Path(path).open('rb') as handle:
        for chunk in iter(lambda:handle.read(4*1024*1024), b''):
            h.update(chunk)
    return h.hexdigest()

def safe_member(info):
    name = info.filename.replace('\\', '/')
    p = PurePosixPath(name)
    if p.is_absolute() or '..' in p.parts or any(':' in x for x in p.parts):
        raise ValueError('Unsafe archive member: '+name)
    if stat.S_ISLNK(info.external_attr >> 16):
        raise ValueError('Symlinks are not accepted in the data archive')
    return name

def prepare(archive, workspace):
    archive, workspace = Path(archive).resolve(), Path(workspace).resolve()
    expected = json.loads((ROOT/'config/data_archive.json').read_text(encoding='utf-8'))
    if sha256(archive) != expected['sha256']:
        raise ValueError('Archive SHA-256 differs from the recorded revision data package.')
    if workspace.exists():
        raise FileExistsError('Choose a new workspace; existing data will not be overwritten.')
    workspace.parent.mkdir(parents=True, exist_ok=True)
    with tempfile.TemporaryDirectory(prefix='ms-data-', dir=workspace.parent) as temporary:
        stage = Path(temporary)
        with zipfile.ZipFile(archive) as z:
            for info in z.infolist():
                name = safe_member(info)
                target = stage/name
                if info.is_dir():
                    target.mkdir(parents=True, exist_ok=True)
                else:
                    target.parent.mkdir(parents=True, exist_ok=True)
                    with z.open(info) as src, target.open('wb') as dst:
                        shutil.copyfileobj(src, dst)
        manifest = json.loads((stage/'S1_Data_manifest.json').read_text(encoding='utf-8'))
        for row in manifest:
            file = stage/row['archive_path']
            if not file.is_file() or file.stat().st_size!=row['bytes'] or sha256(file)!=row['sha256']:
                raise ValueError('Data manifest mismatch: '+row['archive_path'])
        for folder in ['audit','single_nucleus','machine_learning/original_expression_sensitivity',
                       'immune','drug','spatial/source_h5ad','figures/source','figures/illustrator','reports']:
            (stage/folder).mkdir(parents=True,exist_ok=True)
        for filename in ['original_R4_target_checkpoints.txt']:
            shutil.copy2(ROOT/'resources'/filename,stage/'audit'/filename)
        stage.rename(workspace)
    print(f'Prepared {len(manifest)} verified data files in {workspace}')
    return 0

def main():
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument('--archive',type=Path,required=True)
    p.add_argument('--workspace',type=Path,required=True)
    a=p.parse_args()
    return prepare(a.archive,a.workspace)

if __name__=='__main__':
    raise SystemExit(main())
