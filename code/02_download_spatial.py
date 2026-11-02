import os
from pathlib import Path
import json, urllib.request, hashlib, concurrent.futures, time
out = Path(os.environ.get('REVISION_ROOT',str(Path.cwd())))
files = [f for h in json.loads((out / 'audit/hca_files.txt').read_text())['hits'] for f in h['files'] if f['name'].startswith('visium_') or f['name'] == 'ctype_SC.h5ad']
dest = out / 'spatial/source_h5ad'
dest.mkdir(parents=True, exist_ok=True)

def fetch(f):
    p = dest / f['name']
    if p.exists() and p.stat().st_size == f['size'] and hashlib.sha256(p.read_bytes()).hexdigest() == f['sha256']:
        return {'file': f['name'], 'status': 'present', 'bytes': p.stat().st_size}
    for attempt in range(3):
        try:
            req = urllib.request.Request(f['azul_url'], headers={'User-Agent': 'research-revision'})
            with urllib.request.urlopen(req, timeout=90) as r:
                if 'json' in r.headers.get('Content-Type', ''):
                    info = json.load(r)
                    if 'Location' in info:
                        r2 = urllib.request.urlopen(info['Location'], timeout=90)
                        r = r2
                    else:
                        raise RuntimeError(str(info)[:300])
                h = hashlib.sha256()
                n = 0
                with p.with_suffix('.partial').open('wb') as w:
                    while True:
                        b = r.read(4 * 1024 * 1024)
                        if not b:
                            break
                        w.write(b)
                        h.update(b)
                        n += len(b)
            if n != f['size'] or h.hexdigest() != f['sha256']:
                raise ValueError(f'Checksum or size mismatch {n}')
            p.with_suffix('.partial').replace(p)
            result = {'file': f['name'], 'status': 'verified', 'bytes': n, 'sha256': h.hexdigest(), 'source': f['azul_url']}
            print(json.dumps(result), flush=True)
            return result
        except Exception as e:
            err = str(e)
            print(f['name'], 'attempt', attempt + 1, err[:180], flush=True)
    return {'file': f['name'], 'status': 'failed', 'error': err}
with concurrent.futures.ThreadPoolExecutor(max_workers=3) as pool:
    results = list(pool.map(fetch, files))
(out / 'audit/public_spatial_downloads.json').write_text(json.dumps(results, indent=2), encoding='utf-8')
print('DOWNLOAD_COMPLETE', sum((r['status'] != 'failed' for r in results)), '/', len(results), flush=True)

if any(r['status'] == 'failed' for r in results):
    raise SystemExit('One or more spatial downloads failed checksum verification.')
