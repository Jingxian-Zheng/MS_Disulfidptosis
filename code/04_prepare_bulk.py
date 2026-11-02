import os
from pathlib import Path
import pandas as pd, numpy as np, csv, gzip, json, io
out = Path(os.environ.get('REVISION_ROOT', str(Path.cwd())))
base = Path(os.environ.get('BULK_INPUT_ROOT', 'public_bulk_probe_inputs'))

def matrix_meta(p):
    op = gzip.open if p.suffix == '.gz' else open
    meta = {}
    buf = []
    inside = False
    with op(p, 'rt', encoding='utf-8') as f:
        for line in f:
            if line.startswith('!series_matrix_table_begin'):
                inside = True
                continue
            if line.startswith('!series_matrix_table_end'):
                break
            if inside:
                buf.append(line)
            elif line.startswith('!Sample_'):
                row = next(csv.reader([line], delimiter='\t'))
                meta.setdefault(row[0], []).append(row[1:])
    return (pd.read_csv(io.StringIO(''.join(buf)), sep='\t', index_col=0), meta)
train = pd.read_csv(base / 'GSE17048_matrix_exp.txt', sep='\t', index_col=0)
brain = pd.read_csv(base / 'GSE32915_matrix_exp.txt', sep='\t', index_col=0)
ext, emeta = matrix_meta(out / 'runtime/GSE203241_series_matrix.txt.gz')
_, tmeta = matrix_meta(base / 'GSE17048_series_matrix.txt')
_, bmeta = matrix_meta(base / 'GSE32915_series_matrix.txt')
with gzip.open(out / 'runtime/GPL571.annot.gz', 'rt') as f:
    rows = []
    inside = False
    for line in f:
        if line.startswith('!platform_table_begin'):
            inside = True
            continue
        if line.startswith('!platform_table_end'):
            break
        if inside:
            rows.append(line)
annot = pd.read_csv(io.StringIO(''.join(rows)), sep='\t', dtype=str)
maps = {'GSE17048': pd.read_csv(base / 'GSE17048_GPL6947.txt', sep='\t', dtype=str), 'GSE32915': pd.read_csv(base / 'GSE32915_GPL4133.txt', sep='\t', dtype=str), 'GSE203241': annot[['ID', 'Gene symbol']].rename(columns={'Gene symbol': 'gene_name'})}
mats = {}
log = []
for ds, x in [('GSE17048', train), ('GSE32915', brain), ('GSE203241', ext)]:
    x.index = x.index.astype(str)
    mp = maps[ds].dropna().drop_duplicates('ID')
    mp = mp[mp.gene_name.str.match('^[A-Za-z0-9_.-]+$')]
    mp = mp[mp.ID.isin(x.index)]
    xx = x.loc[mp.ID].astype(float)
    xx.index = mp.gene_name
    q99 = np.nanquantile(xx.to_numpy(), 0.99)
    linear = bool(q99 > 50)
    n_floor = int((xx < 1).sum().sum()) if linear else 0
    if linear:
        xx = np.log2(xx.clip(lower=1))
    mats[ds] = xx.groupby(level=0, sort=True).mean()
    mats[ds].to_csv(out / 'machine_learning' / (ds + '_gene_log2.tsv'), sep='\t')
    log.append({'dataset': ds, 'n_samples': x.shape[1], 'n_probes': x.shape[0], 'n_genes': mats[ds].shape[0], 'q99': float(q99), 'log2_transformed': linear, 'values_floored_to_1_before_log2': n_floor})
mrna = set(pd.read_csv(Path(os.environ.get('LEGACY_STUDY_ROOT', str(out / 'legacy_inputs'))) / 'Data/GSE17048_mRNA_exp.tsv', sep='\t', index_col=0).index)
universe = sorted(set.intersection(*(set(m.index) for m in mats.values())) & mrna)
(out / 'machine_learning/frozen_rank_gene_universe.txt').write_text('\n'.join(universe), encoding='utf-8')
features = (Path(os.environ.get('LEGACY_STUDY_ROOT', str(out / 'legacy_inputs'))) / 'Result1/disulfidptosis_markers.txt').read_text().split()
for ds, x in mats.items():
    missing = set(features) - set(x.index)
    print(ds, 'missing', ','.join(sorted(missing)), 'universe', len(universe))
    rank = x.loc[universe].rank(axis=0, method='average', pct=True)
    feat = rank.reindex(features).T
    feat.index.name = 'sample'
    feat.to_csv(out / 'machine_learning' / (ds + '_features.csv'))

def sample_meta(meta):
    ids = meta['!Sample_geo_accession'][0]
    titles = meta['!Sample_title'][0]
    d = pd.DataFrame({'sample': ids, 'title': titles})
    for j, vals in enumerate(meta.get('!Sample_characteristics_ch1', [])):
        if len(vals) == len(d):
            d['characteristic_' + str(j + 1)] = vals
    return d
for ds, meta in [('GSE17048', tmeta), ('GSE32915', bmeta), ('GSE203241', emeta)]:
    d = sample_meta(meta)
    if ds == 'GSE203241':
        d['group'] = np.where(d.title.str.startswith(('POMS', 'AOMS')), 'MS', 'Control')
        d['age_stratum'] = np.where(d.title.str.startswith(('POMS', 'PHC')), 'Pediatric', 'Adult')
        d['donor'] = d.title.str.extract('\\[([^\\]]+)\\]')[0]
        for c in d.columns:
            if c.startswith('characteristic') and d[c].str.startswith('age (years):').any():
                d['age'] = pd.to_numeric(d[c].str.split(':').str[-1].str.strip())
    else:
        orig = pd.read_csv(Path(os.environ.get('LEGACY_STUDY_ROOT', str(out / 'legacy_inputs'))) / 'Data' / (ds + '_sample_type.txt'), sep='\t')
        d = d.merge(orig, on='sample', validate='one_to_one')
        if ds == 'GSE32915':
            d['donor'] = d.title.str.replace('(_new|_\\d+)$', '', regex=True).str.replace(' ', '_', regex=False)
        else:
            d['donor'] = d['sample']
            d['clinical_subtype'] = d.title.str.extract('(PPMS|RRMS|SPMS)')[0]
    d.to_csv(out / 'machine_learning' / (ds + '_metadata.csv'), index=False)
(out / 'audit/bulk_preprocessing.json').write_text(json.dumps(log, indent=2), encoding='utf-8')
print(json.dumps(log))
print('PREPROCESSING_COMPLETE')
