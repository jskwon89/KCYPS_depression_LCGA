# -*- coding: utf-8 -*-
"""Check ICC manually since MixedLM gave singular variance."""
import pandas as pd
import numpy as np
import os

ROOT = r'D:/2026/SCI/청소년패널조사'
OUT = os.path.join(ROOT, 'output', 'desc')
long_df = pd.read_parquet(os.path.join(OUT, 'long_panel.parquet'))

def icc1(group_var, value):
    """ICC(1,1) from one-way random effects ANOVA: ICC = (MS_between - MS_within) / (MS_between + (k_bar-1)*MS_within)."""
    df = pd.DataFrame({'g': group_var, 'y': value}).dropna()
    ng = df.groupby('g').size()
    k_bar = (df.shape[0] - (ng**2).sum() / df.shape[0]) / (ng.shape[0] - 1)  # harmonic-ish avg
    # Manual one-way ANOVA
    grand = df['y'].mean()
    ssb = (ng * (df.groupby('g')['y'].mean() - grand)**2).sum()
    ssw = ((df['y'] - df.groupby('g')['y'].transform('mean'))**2).sum()
    df_between = ng.shape[0] - 1
    df_within = df.shape[0] - ng.shape[0]
    msb = ssb / df_between
    msw = ssw / df_within
    icc = (msb - msw) / (msb + (k_bar - 1) * msw)
    return icc, msb, msw, k_bar

results = []
for s in ['depression', 'life_sat', 'self_esteem', 'parent_warmth', 'parent_reject',
          'peer_attach', 'teacher_rel', 'sm_addiction', 'aggression', 'somatic',
          'withdrawal', 'attention', 'happy', 'grit', 'parent_autonomy', 'parent_coerce',
          'parent_structure', 'parent_chaos']:
    icc, msb, msw, k = icc1(long_df['ID'], long_df[s])
    results.append({'scale': s, 'ICC1': icc, 'MS_between': msb, 'MS_within': msw, 'k_bar': k})

icc_df = pd.DataFrame(results).sort_values('ICC1', ascending=False)
icc_df.to_csv(os.path.join(OUT, 'icc_manual.csv'), index=False, encoding='utf-8-sig')
print(icc_df.round(3).to_string(index=False))
