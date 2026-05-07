# -*- coding: utf-8 -*-
"""
Build a HIGH-LEVEL construct inventory:
  prefix (e.g., YPSY1, YPRT2, YPST1) -> # items, # waves with item, sample construct label
This is the basis for designing analysis: pick constructs that have:
  (a) >=4 waves of coverage,
  (b) clear validated scale meaning,
  (c) appropriate level of measurement.
"""
import pandas as pd
import numpy as np
import os, re, json

ROOT = r'D:/2026/SCI/청소년패널조사'
OUT = os.path.join(ROOT, 'output', 'codebook')

cov = pd.read_csv(os.path.join(OUT, 'youth_root_coverage.csv'), encoding='utf-8-sig')

# Extract prefix - first 4 alphabet chars + first numeric (e.g., YPSY1)
def construct_prefix(root):
    m = re.match(r'^([A-Z]+\d+)', str(root))
    return m.group(1) if m else str(root)[:4]

cov['prefix'] = cov['root'].apply(construct_prefix)

# Aggregate
agg = cov.groupby('prefix').agg(
    n_items=('root', 'count'),
    n_in_all_7=('n_waves', lambda s: (s == 7).sum()),
    n_in_6plus=('n_waves', lambda s: (s >= 6).sum()),
    n_in_4plus=('n_waves', lambda s: (s >= 4).sum()),
    sample_construct=('construct', lambda s: s.dropna().iloc[0] if len(s.dropna()) else ''),
).reset_index()
agg = agg.sort_values('n_items', ascending=False)
agg.to_csv(os.path.join(OUT, 'construct_prefix_summary.csv'), index=False, encoding='utf-8-sig')
print(f'Total prefixes: {len(agg)}')
print('\n--- Top 40 prefixes by item count ---')
print(agg.head(40).to_string())

# Also: filter to constructs that look like Likert scales with multi-item
# (n_items >= 5, in_all_7 == n_items meaning fully balanced panel)
balanced = cov[cov.n_waves == 7]
balanced_agg = balanced.groupby('prefix').agg(
    n_items=('root', 'count'),
    items=('root', lambda s: ';'.join(sorted(s))),
    sample_construct=('construct', lambda s: s.dropna().iloc[0] if len(s.dropna()) else ''),
).reset_index().sort_values('n_items', ascending=False)
balanced_agg.to_csv(os.path.join(OUT, 'construct_prefix_balanced.csv'), index=False, encoding='utf-8-sig')
print(f'\n--- Top 40 BALANCED prefixes (all items in 7 waves) ---')
print(balanced_agg.head(40).to_string())
