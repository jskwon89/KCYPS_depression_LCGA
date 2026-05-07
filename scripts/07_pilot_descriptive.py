# -*- coding: utf-8 -*-
"""
Pilot descriptive analysis:
  (1) attrition by baseline characteristics
  (2) cross-wave correlations of key DVs
  (3) ICC and within/between variance decomposition
  (4) item-level psychometrics for depression scale (key DV)
"""
import pandas as pd
import numpy as np
import os, json
from scipy import stats

ROOT = r'D:/2026/SCI/청소년패널조사'
OUT = os.path.join(ROOT, 'output', 'desc')

long_df = pd.read_parquet(os.path.join(OUT, 'long_panel.parquet'))
wide = pd.read_parquet(os.path.join(OUT, 'wide_panel.parquet'))

print('=== 1. RESPONSE PATTERN ===')
# Define "responded" if depression has valid value
resp = long_df.assign(responded=long_df.depression.notna()).pivot_table(index='ID', columns='wave', values='responded', fill_value=False)
resp.columns = [f'w{c}' for c in resp.columns]
resp['n_waves_responded'] = resp.sum(axis=1)
print('Distribution of N waves responded:')
print(resp.n_waves_responded.value_counts().sort_index().to_string())

# Attrition pattern
patt = resp[[f'w{w}' for w in range(1,8)]].astype(int).astype(str).agg(''.join, axis=1)
print('\nTop 10 response patterns:')
print(patt.value_counts().head(10).to_string())

# Save
resp.to_csv(os.path.join(OUT, 'response_pattern.csv'), encoding='utf-8-sig')

print('\n=== 2. CROSS-WAVE CORRELATIONS (depression) ===')
dep_wide = wide[[f'depression_w{w}' for w in range(1,8)]]
print(dep_wide.corr(method='pearson').round(3).to_string())

print('\n=== 3. ICC (variance decomposition) ===')
# ICC = between-person var / total var
import statsmodels.formula.api as smf
import statsmodels.api as sm

icc_results = []
for scale in ['depression', 'life_sat', 'self_esteem', 'parent_warmth', 'parent_reject',
              'peer_attach', 'teacher_rel', 'sm_addiction', 'aggression', 'somatic',
              'withdrawal', 'attention', 'happy', 'grit']:
    ldf = long_df[['ID', 'wave', scale]].dropna()
    if len(ldf) < 100: continue
    try:
        model = sm.MixedLM(ldf[scale], np.ones((len(ldf), 1)), groups=ldf['ID'])
        res = model.fit(method='lbfgs', disp=False)
        var_re = float(res.cov_re.iloc[0,0])
        var_resid = float(res.scale)
        icc = var_re / (var_re + var_resid)
        icc_results.append({'scale': scale, 'between_var': var_re, 'within_var': var_resid, 'icc': icc, 'n_obs': len(ldf)})
    except Exception as e:
        print(f'  {scale}: ICC failed - {e}')
icc_df = pd.DataFrame(icc_results).sort_values('icc', ascending=False)
icc_df.to_csv(os.path.join(OUT, 'icc.csv'), index=False, encoding='utf-8-sig')
print(icc_df.round(3).to_string())

print('\n=== 4. ATTRITION ANALYSIS ===')
# Compare wave-1 baseline characteristics by 7-wave completion
w1_baseline = long_df[long_df.wave==1].set_index('ID')
completion = resp[['n_waves_responded']]
m = w1_baseline.join(completion, how='left')
m['completed_all_7'] = (m['n_waves_responded'] == 7).astype(int)
print(f'Completed all 7 waves: {m.completed_all_7.sum()} of {len(m)} = {m.completed_all_7.mean():.1%}')

attr_compare = []
for var in ['depression', 'life_sat', 'self_esteem', 'parent_warmth', 'parent_reject',
            'peer_attach', 'sm_addiction', 'YGENDER']:
    if var not in m.columns:
        continue
    g0 = m[m.completed_all_7 == 0][var].dropna()
    g1 = m[m.completed_all_7 == 1][var].dropna()
    if len(g0) < 5 or len(g1) < 5: continue
    t, p = stats.ttest_ind(g0, g1, equal_var=False)
    attr_compare.append({'var': var, 'mean_dropout': g0.mean(), 'mean_complete': g1.mean(), 't': t, 'p': p})
print(pd.DataFrame(attr_compare).round(3).to_string(index=False))

print('\n=== 5. DEPRESSION ITEM-LEVEL DESCRIPTIVES (KEY DV) ===')
# Re-load wave 1 file to get item-level
import pyreadstat
df1, _ = pyreadstat.read_dta(os.path.join(ROOT, 'KCYPS2018m1Yw1.dta'))
dep_items_w1 = [f'YPSY4E{i:02d}w1' for i in range(1,11)]
dep_w1 = df1[dep_items_w1].apply(lambda x: x.where(x>=1, np.nan)).where(lambda x: x<=4)
print('Wave 1 depression item descriptives (1=전혀 그렇지 않다 ~ 4=매우 그렇다):')
print(dep_w1.describe().T.round(3)[['count','mean','std','min','max']].to_string())
print('\nSuicidal ideation (E04 "죽고 싶은 생각이 든다") wave 1:')
print((df1['YPSY4E04w1'].where((df1['YPSY4E04w1']>=1)&(df1['YPSY4E04w1']<=4))).value_counts(dropna=False).sort_index().to_string())

# 자살사고 발생률 by 차수
print('\nSuicidal ideation (E04) endorsement rates per wave (>=2 = some endorsement):')
for w in range(1, 8):
    df, _ = pyreadstat.read_dta(os.path.join(ROOT, f'KCYPS2018m1Yw{w}.dta'))
    col = f'YPSY4E04w{w}'
    s = df[col].where((df[col]>=1)&(df[col]<=4))
    n = s.notna().sum()
    p_any = (s>=2).sum()
    p_strong = (s>=3).sum()
    print(f'  wave {w}: N={n}, any (>=2): {p_any} ({p_any/n*100:.1f}%), strong (>=3): {p_strong} ({p_strong/n*100:.1f}%)')

print('\n=== 6. AGE INFERENCE ===')
# Wave year - birth year ~ age. m1 cohort: 중1 = age 13~14 (Korean age) at 2018
df1, _ = pyreadstat.read_dta(os.path.join(ROOT, 'KCYPS2018m1Yw1.dta'))
print('YBRT1Aw1 (birth year) distribution:')
print(df1['YBRT1Aw1'].value_counts(dropna=False).sort_index().to_string())
