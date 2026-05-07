# -*- coding: utf-8 -*-
"""
Build clean codebook inventory: variable -> construct label.
Then merge with the actual variable names in dta files (which have ...w{N} suffix).
Finally, output a wave-coverage matrix with construct labels.
"""
import pandas as pd
import numpy as np
import os, re, json

ROOT = r'D:/2026/SCI/청소년패널조사'
OUT = os.path.join(ROOT, 'output', 'codebook')

def parse_cb(name):
    """Parse a CB sheet. Cols are unnamed but standard layout: [_, Q_block, sub1, sub2, var, value, value_label]."""
    df = pd.read_csv(os.path.join(OUT, f'{name}.csv'), header=2, encoding='utf-8-sig')
    # After header=2 the actual columns become row index 4 ('문항' etc.)
    # Easier: read raw and find columns by content
    raw = pd.read_csv(os.path.join(OUT, f'{name}.csv'), header=None, encoding='utf-8-sig')
    # The file has Korean header at row 3 (index 3 because of pre-rows). Let's find it.
    for i in range(min(10, len(raw))):
        row = raw.iloc[i].astype(str).fillna('').tolist()
        if '변수명' in row:
            hdr_idx = i
            break
    else:
        raise RuntimeError(f'No header found in {name}')
    df = pd.read_csv(os.path.join(OUT, f'{name}.csv'), header=hdr_idx, encoding='utf-8-sig')
    # Remove unnamed leading column if all NaN
    df.columns = [c.strip() if isinstance(c, str) else c for c in df.columns]
    return df

def variable_dict(name):
    df = parse_cb(name)
    # Identify columns
    var_col = next((c for c in df.columns if c == '변수명'), None)
    item_col = next((c for c in df.columns if c == '문항'), None)
    val_col = next((c for c in df.columns if c == '값'), None)
    lab_col = next((c for c in df.columns if c == '값설명'), None)
    # subgroup columns are often 2 unnamed cols between '문항' and '변수명'
    print(f'{name}: cols={df.columns.tolist()}')
    cols = df.columns.tolist()
    # Find any unnamed sub columns between item and var
    if var_col and item_col:
        i_item = cols.index(item_col)
        i_var = cols.index(var_col)
        sub_cols = [c for c in cols[i_item+1:i_var]]
    else:
        sub_cols = []
    # Forward fill item and sub columns within blocks
    if item_col:
        df[item_col] = df[item_col].ffill()
    for c in sub_cols:
        df[c] = df[c].ffill()
    # Build var -> construct dict
    rows = []
    cur_var = None
    for _, r in df.iterrows():
        v = r.get(var_col)
        if isinstance(v, str) and v.strip():
            cur_var = v.strip()
            item = r.get(item_col, '')
            sub = ' / '.join([str(r.get(c, '')) for c in sub_cols if isinstance(r.get(c, ''), str) and r.get(c, '')])
            rows.append({
                'var': cur_var,
                'construct': str(item).replace('\n', ' ').strip() if isinstance(item, str) else '',
                'sub': sub,
                'first_val': r.get(val_col, ''),
                'first_val_label': r.get(lab_col, ''),
            })
    out = pd.DataFrame(rows).drop_duplicates('var')
    return out

inventory = {}
for n in ['CB_m', 'CB_h', 'CB_p']:
    inv = variable_dict(n)
    inventory[n] = inv
    inv.to_csv(os.path.join(OUT, f'inv_{n}.csv'), index=False, encoding='utf-8-sig')
    print(f'{n}: {len(inv)} variables')

# Combine youth (m + h) since same person uses m form when in middle school, h form when in high school
inv_y = pd.concat([inventory['CB_m'].assign(form='m'), inventory['CB_h'].assign(form='h')], ignore_index=True)
inv_y.to_csv(os.path.join(OUT, 'inv_youth_all.csv'), index=False, encoding='utf-8-sig')

# Dedup by var, prefer first occurrence with non-empty construct
dedup_y = inv_y.sort_values('construct', key=lambda s: s.fillna('').str.len(), ascending=False).drop_duplicates('var', keep='first')
dedup_y.to_csv(os.path.join(OUT, 'inv_youth_dedup.csv'), index=False, encoding='utf-8-sig')
print(f'Combined youth dedup: {len(dedup_y)} unique vars')

# Now build cross-wave coverage: for each youth dta variable (with w{N}), strip suffix and look up construct
# Read variable lists from previously saved files
y_dfs = []
for w in range(1, 8):
    f = os.path.join(OUT, f'youth_w{w}_vars.csv')
    if not os.path.exists(f): continue
    d = pd.read_csv(f, encoding='utf-8-sig')
    d['wave'] = w
    y_dfs.append(d)
all_y = pd.concat(y_dfs)
# Strip wave suffix
def strip_wave(v):
    m = re.match(r'^(.*?)w\d+$', str(v))
    return m.group(1) if m else str(v)
all_y['root'] = all_y['var'].apply(strip_wave)

# Pivot: root vs wave
cov = all_y.pivot_table(index='root', columns='wave', values='var', aggfunc='count', fill_value=0)
cov.columns = [f'w{c}' for c in cov.columns]
cov['n_waves'] = cov.sum(axis=1)
# Merge construct labels
inv_lookup = dedup_y.set_index('var')
cov = cov.reset_index()
cov['construct'] = cov['root'].map(lambda v: inv_lookup['construct'].get(v, ''))
cov['sub'] = cov['root'].map(lambda v: inv_lookup['sub'].get(v, ''))
cov['form'] = cov['root'].map(lambda v: inv_lookup['form'].get(v, ''))
cov = cov.sort_values(['n_waves', 'root'], ascending=[False, True])
cov.to_csv(os.path.join(OUT, 'youth_root_coverage.csv'), index=False, encoding='utf-8-sig')

print(f'\nYouth root coverage: {len(cov)} root variables')
print(f'  in 7 waves: {(cov.n_waves == 7).sum()}')
print(f'  in 6 waves: {(cov.n_waves == 6).sum()}')
print(f'  in 5 waves: {(cov.n_waves >= 5).sum()}')
print(f'  in 4 waves: {(cov.n_waves >= 4).sum()}')

print('\n--- Roots in 6+ waves with construct ---')
print(cov[cov.n_waves >= 6][['root', 'n_waves', 'construct', 'sub', 'form']].head(60).to_string())
