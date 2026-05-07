# -*- coding: utf-8 -*-
"""
Inspect KCYPS dta files: shape, ID, variable labels.
.dta files have variable labels (Korean) and value labels embedded - safer than xlsx.
"""
import pandas as pd
import numpy as np
import os, json, sys

ROOT = r'D:/2026/SCI/청소년패널조사'
OUT = os.path.join(ROOT, 'output', 'codebook')
os.makedirs(OUT, exist_ok=True)

# Use pyreadstat for richer metadata; fallback to pandas read_stata
try:
    import pyreadstat
    HAVE_PYREAD = True
    print('pyreadstat available')
except Exception as e:
    HAVE_PYREAD = False
    print('pyreadstat not available, using pandas:', e)

waves = list(range(1, 8))
files_y = [os.path.join(ROOT, f'KCYPS2018m1Yw{w}.dta') for w in waves]
files_p = [os.path.join(ROOT, f'KCYPS2018m1Pw{w}.dta') for w in waves]

summary = {'youth': {}, 'parent': {}}
all_var_labels = {}  # var -> {wave: (label, type, n_unique)}

def inspect_file(path, kind, wave):
    print(f'  reading {os.path.basename(path)} ...')
    if HAVE_PYREAD:
        df, meta = pyreadstat.read_dta(path)
        var_labels = meta.column_names_to_labels
        val_labels = meta.variable_value_labels
    else:
        df = pd.read_stata(path, convert_categoricals=False)
        var_labels = {}
        val_labels = {}
    return df, var_labels, val_labels

# Youth files
print('\n=== YOUTH FILES ===')
for w, p in zip(waves, files_y):
    if not os.path.exists(p):
        print(f'  missing: {p}')
        continue
    df, vl, valL = inspect_file(p, 'youth', w)
    summary['youth'][w] = {
        'n_rows': int(df.shape[0]),
        'n_cols': int(df.shape[1]),
        'columns_first15': df.columns.tolist()[:15],
        'has_id': 'ID' in df.columns or 'id' in df.columns,
    }
    for c in df.columns:
        all_var_labels.setdefault(c, {})[f'Y{w}'] = (vl.get(c, ''), str(df[c].dtype), int(df[c].nunique(dropna=True)))
    # Save var label list for this wave
    pd.DataFrame({
        'var': df.columns,
        'label': [vl.get(c, '') for c in df.columns],
        'dtype': [str(df[c].dtype) for c in df.columns],
        'n_unique': [df[c].nunique(dropna=True) for c in df.columns],
        'n_missing': [int(df[c].isna().sum()) for c in df.columns],
    }).to_csv(os.path.join(OUT, f'youth_w{w}_vars.csv'), index=False, encoding='utf-8-sig')

print('\n=== PARENT FILES ===')
for w, p in zip(waves, files_p):
    if not os.path.exists(p):
        print(f'  missing: {p}')
        continue
    df, vl, valL = inspect_file(p, 'parent', w)
    summary['parent'][w] = {
        'n_rows': int(df.shape[0]),
        'n_cols': int(df.shape[1]),
        'columns_first15': df.columns.tolist()[:15],
    }
    pd.DataFrame({
        'var': df.columns,
        'label': [vl.get(c, '') for c in df.columns],
        'dtype': [str(df[c].dtype) for c in df.columns],
        'n_unique': [df[c].nunique(dropna=True) for c in df.columns],
        'n_missing': [int(df[c].isna().sum()) for c in df.columns],
    }).to_csv(os.path.join(OUT, f'parent_w{w}_vars.csv'), index=False, encoding='utf-8-sig')

# Build wave coverage matrix
print('\n=== Building variable coverage matrix ...')
var_set = sorted(all_var_labels.keys())
rows = []
for v in var_set:
    info = all_var_labels[v]
    label = next((info[k][0] for k in info if info[k][0]), '')
    waves_present = [k for k in info]
    rows.append({
        'var': v,
        'label': label,
        'n_waves': len(waves_present),
        'waves': ';'.join(sorted(waves_present)),
    })
cov = pd.DataFrame(rows).sort_values(['n_waves', 'var'], ascending=[False, True])
cov.to_csv(os.path.join(OUT, 'variable_coverage_youth.csv'), index=False, encoding='utf-8-sig')
print(f'Coverage saved: {cov.shape[0]} unique vars across youth files')
print(f'Vars in all 7 waves: {(cov.n_waves == 7).sum()}')
print(f'Vars in >= 5 waves: {(cov.n_waves >= 5).sum()}')

with open(os.path.join(OUT, 'summary.json'), 'w', encoding='utf-8') as f:
    json.dump(summary, f, ensure_ascii=False, indent=2)
print('\nDONE')
