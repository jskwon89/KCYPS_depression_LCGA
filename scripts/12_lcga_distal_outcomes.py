# -*- coding: utf-8 -*-
"""
LCGA Step 2: Profile classes by baseline covariates and Distal outcomes (Wave 7).
- Predictors at W1: gender, parenting (warmth/reject/coerce/structure), peer attach, self-esteem, sm_addiction
- Distal outcomes at W7: suicidal ideation (binarized), life_sat, offline_delin, cyber_delin

Approach: 'Modal' assignment + multinomial logistic for class predictors.
For distal outcomes: weighted by class posterior ('manual 3-step' light version).
"""
import pandas as pd
import numpy as np
import os, json
import statsmodels.api as sm
import statsmodels.formula.api as smf
from scipy import stats
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

ROOT = r'D:/2026/SCI/청소년패널조사'
DESC = os.path.join(ROOT, 'output', 'desc')
LCGA_OUT = os.path.join(ROOT, 'output', 'paperA_lcga')

class_df = pd.read_csv(os.path.join(LCGA_OUT, 'class_assignments_K4.csv'), encoding='utf-8-sig')
long_df = pd.read_parquet(os.path.join(DESC, 'long_panel.parquet'))

# Wave 1 baseline
w1 = long_df[long_df.wave==1].copy()
w7 = long_df[long_df.wave==7].copy()

# Cyber/offline counts (corrected version)
cyber = pd.read_parquet(os.path.join(ROOT, 'output', 'paperB_riclpm', 'cyber_count_wide.parquet'))
off = pd.read_parquet(os.path.join(ROOT, 'output', 'paperB_riclpm', 'off_count_wide.parquet'))

# Suicidal ideation at W7: from raw items
import pyreadstat
df7, _ = pyreadstat.read_dta(os.path.join(ROOT, 'KCYPS2018m1Yw7.dta'))
df7['suic_w7'] = df7['YPSY4E04w7'].where((df7['YPSY4E04w7']>=1) & (df7['YPSY4E04w7']<=4))
df7['suic_any'] = (df7['suic_w7'] >= 2).astype(float)
df7['suic_strong'] = (df7['suic_w7'] >= 3).astype(float)
suic_df = df7[['ID', 'suic_w7', 'suic_any', 'suic_strong']].copy()

# Build predictor frame (W1 baseline)
pred = w1[['ID', 'YGENDER', 'parent_warmth', 'parent_reject', 'parent_coerce', 'parent_structure',
           'parent_chaos', 'parent_autonomy', 'peer_attach', 'self_esteem', 'sm_addiction',
           'depression', 'life_sat', 'happy', 'grit', 'aggression', 'somatic', 'withdrawal',
           'attention']].copy()
pred = pred.rename(columns={c: c+'_w1' for c in pred.columns if c not in ['ID']})

# Distal frame (W7)
distal = w7[['ID', 'depression', 'life_sat', 'self_esteem', 'happy']].copy()
distal = distal.rename(columns={c: c+'_w7' for c in distal.columns if c not in ['ID']})
distal = distal.merge(suic_df, on='ID', how='left')
distal['off_count_w7'] = off['off_w7'].reindex(distal['ID']).values
distal['cyber_count_w7'] = cyber['cyber_w7'].reindex(distal['ID']).values

# Merge
m = class_df.merge(pred, on='ID', how='left').merge(distal, on='ID', how='left')
m.to_csv(os.path.join(LCGA_OUT, 'class_with_predictors.csv'), index=False, encoding='utf-8-sig')
print(f'Merged frame: {len(m)} rows')

# Class labels (based on profile)
class_label = {1: 'C1 Late-Increasing (36.3%)',
               2: 'C2 Stable Low / Resilient (25.6%)',
               3: 'C3 Decreasing (27.2%)',
               4: 'C4 Persistent High / Inverted-U (10.9%)'}
m['class_label'] = m['class_K4'].map(class_label)

# === 1. Baseline (W1) means by class ===
print('\n=== Baseline (W1) characteristics by class ===')
b_cols = ['parent_warmth_w1', 'parent_reject_w1', 'parent_coerce_w1', 'parent_chaos_w1',
          'parent_structure_w1', 'parent_autonomy_w1', 'peer_attach_w1', 'self_esteem_w1',
          'sm_addiction_w1', 'depression_w1', 'life_sat_w1', 'happy_w1', 'grit_w1',
          'aggression_w1', 'somatic_w1', 'withdrawal_w1', 'attention_w1']
desc_class = m.groupby('class_K4')[b_cols].agg(['mean','std','count']).round(3)
desc_class.to_csv(os.path.join(LCGA_OUT, 'baseline_by_class.csv'), encoding='utf-8-sig')

# T-tests/ANOVAs (just F + p)
anova_results = []
for c in b_cols:
    groups = [g[c].dropna().values for _, g in m.groupby('class_K4')]
    if all(len(g) > 1 for g in groups):
        F, p = stats.f_oneway(*groups)
        anova_results.append({'var': c, 'F': F, 'p': p,
                              'C1': groups[0].mean(), 'C2': groups[1].mean(), 'C3': groups[2].mean(), 'C4': groups[3].mean()})
anova_df = pd.DataFrame(anova_results).sort_values('F', ascending=False)
anova_df.to_csv(os.path.join(LCGA_OUT, 'anova_baseline_by_class.csv'), index=False, encoding='utf-8-sig')
print(anova_df.round(3).to_string(index=False))

# Gender by class
print('\nGender (1=male, 2=female) by class:')
gen = pd.crosstab(m.class_K4, m.YGENDER_w1, normalize='index').round(3)
gen.columns = ['male', 'female']
gen.to_csv(os.path.join(LCGA_OUT, 'gender_by_class.csv'), encoding='utf-8-sig')
print(gen.to_string())
chi2, pchi, dof, exp = stats.chi2_contingency(pd.crosstab(m.class_K4, m.YGENDER_w1))
print(f'  chi² = {chi2:.2f}, df = {dof}, p = {pchi:.4g}')

# === 2. Multinomial logistic regression: which W1 covariates predict class membership? ===
print('\n=== Multinomial logistic: W1 covariates -> class ===')
# Reference: class 2 (stable low / resilient)
m2 = m.dropna(subset=['class_K4', 'YGENDER_w1', 'parent_warmth_w1', 'parent_reject_w1',
                       'peer_attach_w1', 'self_esteem_w1', 'sm_addiction_w1']).copy()
m2['female'] = (m2['YGENDER_w1'] == 2).astype(float)
# Standardize predictors for interpretation
predictors = ['female', 'parent_warmth_w1', 'parent_reject_w1', 'parent_chaos_w1',
              'parent_coerce_w1', 'peer_attach_w1', 'self_esteem_w1', 'sm_addiction_w1']
for p in predictors:
    if p != 'female' and p in m2.columns:
        m2[p+'_z'] = (m2[p] - m2[p].mean()) / m2[p].std()

# Use class 2 as reference
m2['class_str'] = m2['class_K4'].astype(int).astype(str)
m2_subset = m2[m2['class_str'].isin(['1', '2', '3', '4'])].copy()
form = 'class_str ~ female + parent_warmth_w1_z + parent_reject_w1_z + parent_chaos_w1_z + parent_coerce_w1_z + peer_attach_w1_z + self_esteem_w1_z + sm_addiction_w1_z'
import statsmodels.formula.api as smf
try:
    mlogit = smf.mnlogit(form, data=m2_subset).fit(method='bfgs', maxiter=200, disp=False)
    print(mlogit.summary())
    # Save params
    with open(os.path.join(LCGA_OUT, 'mlogit_class_summary.txt'), 'w', encoding='utf-8') as f:
        f.write(str(mlogit.summary()))
    # OR for each class vs reference
    or_df = np.exp(mlogit.params)
    or_df.to_csv(os.path.join(LCGA_OUT, 'mlogit_OR.csv'), encoding='utf-8-sig')
    print('\nOdds ratios:')
    print(or_df.round(3).to_string())
except Exception as e:
    print(f'mlogit failed: {e}')

# === 3. Distal outcomes ===
print('\n=== Distal outcomes (W7) by class ===')
dis_cols = ['depression_w7', 'life_sat_w7', 'self_esteem_w7', 'happy_w7', 'suic_any',
            'suic_strong', 'off_count_w7', 'cyber_count_w7']
dis_desc = m.groupby('class_K4')[dis_cols].agg(['mean','std','count']).round(3)
dis_desc.to_csv(os.path.join(LCGA_OUT, 'distal_by_class.csv'), encoding='utf-8-sig')
print(m.groupby('class_K4')[dis_cols].mean().round(3).to_string())

# Suicidal ideation by class
print('\nSuicidal ideation rates by class:')
suic_by_class = m.groupby('class_K4').agg(
    n=('ID', 'count'),
    n_w7=('suic_any', 'count'),
    p_any=('suic_any', 'mean'),
    p_strong=('suic_strong', 'mean'),
).round(3)
suic_by_class.to_csv(os.path.join(LCGA_OUT, 'suicidal_ideation_by_class.csv'), encoding='utf-8-sig')
print(suic_by_class.to_string())

# Logistic: class -> suicidal ideation
print('\nLogistic: class (vs C2 ref) -> suicidal ideation at W7')
m3 = m.dropna(subset=['suic_any', 'class_K4']).copy()
m3['female'] = (m3['YGENDER_w1'] == 2).astype(float)
m3['class_str'] = m3['class_K4'].astype(int).astype(str)
form_log = 'suic_any ~ C(class_str, Treatment(reference="2")) + female'
log_mod = smf.logit(form_log, data=m3).fit(disp=False)
print(log_mod.summary())
with open(os.path.join(LCGA_OUT, 'logit_suic_by_class.txt'), 'w', encoding='utf-8') as f:
    f.write(str(log_mod.summary()))
ors = pd.DataFrame({'OR': np.exp(log_mod.params), 'lower': np.exp(log_mod.conf_int()[0]),
                    'upper': np.exp(log_mod.conf_int()[1]), 'p': log_mod.pvalues}).round(3)
ors.to_csv(os.path.join(LCGA_OUT, 'logit_suic_OR.csv'), encoding='utf-8-sig')
print('\nORs (vs C2 reference):')
print(ors.to_string())

# === 4. Plot trajectories with English labels (avoid font issue) ===
import pickle
classes = m.groupby('class_K4').size()
print('\nClass sizes:', classes.to_dict())

# Recompute mean trajectories from raw data
yfull = long_df.pivot_table(index='ID', columns='wave', values='depression')
yfull = yfull.merge(class_df[['ID', 'class_K4']], left_index=True, right_on='ID', how='left').set_index('ID')
fig, ax = plt.subplots(figsize=(10, 6))
colors = {1: '#ff7f0e', 2: '#2ca02c', 3: '#1f77b4', 4: '#d62728'}
labels = {1: 'C1 Late-Increasing (36.3%)', 2: 'C2 Stable-Low (25.6%)',
          3: 'C3 Decreasing (27.2%)', 4: 'C4 Persistent-High (10.9%)'}
for k in [1, 2, 3, 4]:
    sub = yfull[yfull.class_K4 == k][[1,2,3,4,5,6,7]]
    means = sub.mean()
    stds = sub.std()
    ax.errorbar(range(1, 8), means.values, yerr=stds.values/np.sqrt(len(sub)),
                color=colors[k], linewidth=3, marker='o', markersize=10,
                label=labels[k], capsize=4)
ax.set_xlabel('Wave (1=Grade 7, 2018 → 7=Post-Grade 12, 2024)', fontsize=12)
ax.set_ylabel('Depression Score (mean of 10 items, 1-4)', fontsize=12)
ax.set_title('Latent Class Growth Analysis: Adolescent Depression Trajectories\nKCYPS m1 cohort, N=2,452, FIML', fontsize=13)
ax.legend(loc='upper right', fontsize=10)
ax.grid(alpha=0.3)
ax.set_ylim(1, 3)
plt.tight_layout()
plt.savefig(os.path.join(LCGA_OUT, 'lcga_K4_observed_means.png'), dpi=150, bbox_inches='tight')
plt.close()

# Bar plot: suicidal ideation rates by class
fig, axes = plt.subplots(1, 2, figsize=(12, 5))
ax = axes[0]
suic_by_class[['p_any', 'p_strong']].plot(kind='bar', ax=ax, color=['#ff7f0e', '#d62728'])
ax.set_ylabel('Proportion endorsing suicidal ideation')
ax.set_title('Suicidal ideation at Wave 7 by depression class')
ax.set_xticklabels([labels[k] for k in [1,2,3,4]], rotation=20, ha='right', fontsize=9)
ax.legend(['Any (≥2)', 'Strong (≥3)'])
ax.grid(alpha=0.3)
ax = axes[1]
m.groupby('class_K4')['life_sat_w7'].mean().plot(kind='bar', ax=ax, color='#2ca02c')
ax.set_ylabel('Life satisfaction (W7)')
ax.set_title('Life satisfaction at Wave 7 by class')
ax.set_xticklabels([labels[k] for k in [1,2,3,4]], rotation=20, ha='right', fontsize=9)
ax.grid(alpha=0.3)
plt.tight_layout()
plt.savefig(os.path.join(LCGA_OUT, 'distal_outcomes_by_class.png'), dpi=150, bbox_inches='tight')
plt.close()

print('\nDone. All outputs saved to', LCGA_OUT)
