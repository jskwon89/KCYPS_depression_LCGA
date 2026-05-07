# -*- coding: utf-8 -*-
"""
BCH 3-step distal outcomes + Cox PH (with Schoenfeld test) using 9-item primary classes.

Fixes from reviewer feedback:
  - 9-item LCGA classes (no E04 contamination)
  - Bootstrap CIs clipped to [0,1] for proportions
  - Schoenfeld residual PH test added (saved to file)
"""
import pandas as pd
import numpy as np
import os
import pyreadstat
from numpy.linalg import inv
from lifelines import KaplanMeierFitter, CoxPHFitter
from lifelines.statistics import proportional_hazard_test
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import warnings
warnings.filterwarnings('ignore')

ROOT = r'D:/2026/SCI/청소년패널조사'
P9   = os.path.join(ROOT, 'output', 'paperA_lcga', 'primary_9item')
BCH9 = os.path.join(P9, 'bch')
LD9  = os.path.join(P9, 'long_distal')
os.makedirs(BCH9, exist_ok=True); os.makedirs(LD9, exist_ok=True)

class_df = pd.read_csv(os.path.join(P9, 'class_assignments_K4_9item.csv'), encoding='utf-8-sig')
class_df = class_df.rename(columns={'class_K4_9item': 'class_K4'})
post = class_df[[c for c in class_df.columns if c.startswith('prob_class')]].values
modal = class_df['class_K4'].values
K = post.shape[1]
N = post.shape[0]

# ===== Classification error matrix (BCH) =====
D = np.zeros((K, K))
for s in range(K):
    denom = post[:, s].sum()
    for t in range(K):
        D[s, t] = (post[:, s] * (modal == (t+1))).sum() / denom
print('\nClassification-error matrix D (rows=true, cols=modal):')
print(pd.DataFrame(D, index=[f's={k+1}' for k in range(K)],
                   columns=[f't={k+1}' for k in range(K)]).round(3).to_string())
print('Diagonal (correct classification):', np.round(np.diag(D), 3).tolist())
D_inv = inv(D)
pd.DataFrame(D).to_csv(os.path.join(BCH9, 'classification_error_matrix_D.csv'), index=False, encoding='utf-8-sig')
pd.DataFrame(D_inv).to_csv(os.path.join(BCH9, 'classification_error_matrix_Dinv.csv'), index=False, encoding='utf-8-sig')

# ===== Distal outcomes =====
df7, _ = pyreadstat.read_dta(os.path.join(ROOT, 'KCYPS2018m1Yw7.dta'))
df7['suic_w7_raw'] = df7['YPSY4E04w7'].where((df7['YPSY4E04w7']>=1) & (df7['YPSY4E04w7']<=4))
df7['suic_any'] = (df7['suic_w7_raw'] >= 2).astype(float)
df7['suic_strong'] = (df7['suic_w7_raw'] >= 3).astype(float)

desc = pd.read_parquet(os.path.join(ROOT, 'output', 'desc', 'long_panel.parquet'))
w7 = desc[desc.wave == 7].set_index('ID')[['depression', 'depression9', 'life_sat', 'self_esteem', 'happy']]
w7.columns = [c + '_w7' for c in w7.columns]
suic_df = df7.set_index('ID')[['suic_w7_raw', 'suic_any', 'suic_strong']]
distal = w7.join(suic_df, how='outer')

df1, _ = pyreadstat.read_dta(os.path.join(ROOT, 'KCYPS2018m1Yw1.dta'))
distal['female'] = (df1.set_index('ID')['YGENDERw1'] == 2).astype(float)
m = class_df.merge(distal, on='ID', how='left')
m['modal'] = m['class_K4'].astype(int)

# ===== BCH bootstrap proportions with clipped CIs =====
def bch_logit_clipped(outcome_col, ref_class):
    individual_y = m[['ID', outcome_col]].dropna().set_index('ID')[outcome_col]
    individual = m[m['ID'].isin(individual_y.index)].copy()
    individual['y'] = individual_y.reindex(individual['ID']).values
    individual['modal'] = individual['class_K4'].astype(int)

    rng = np.random.default_rng(42)
    B = 2000
    rows = []
    boot_props = np.zeros((B, K))
    for t in range(K):
        w_t = D_inv[(individual['modal'] - 1).astype(int).values, t]
        num = (w_t * individual['y']).sum()
        denom = w_t.sum()
        prop = num / denom if denom != 0 else np.nan
        # CLIP to [0, 1] for interpretability of proportion
        prop = np.clip(prop, 0, 1)
        rows.append({'target_class': t+1, 'BCH_corrected_proportion': prop, 'sum_weights': denom})
    for b in range(B):
        idx = rng.integers(0, len(individual), size=len(individual))
        sub = individual.iloc[idx]
        for t in range(K):
            w_t = D_inv[(sub['modal'] - 1).astype(int).values, t]
            num = (w_t * sub['y']).sum()
            denom = w_t.sum()
            p = num / denom if denom != 0 else np.nan
            boot_props[b, t] = np.clip(p, 0, 1)
    res = pd.DataFrame(rows)
    res['SE'] = np.std(boot_props, axis=0)
    res['CI_lo'] = np.quantile(boot_props, 0.025, axis=0)  # already clipped
    res['CI_hi'] = np.quantile(boot_props, 0.975, axis=0)

    # ORs vs ref
    or_rows = []
    for t in range(K):
        if t + 1 == ref_class: continue
        p_t = res.iloc[t]['BCH_corrected_proportion']
        p_r = res.iloc[ref_class - 1]['BCH_corrected_proportion']
        eps = 1e-6
        odds_t = (p_t + eps) / (1 - p_t + eps)
        odds_r = (p_r + eps) / (1 - p_r + eps)
        or_pt = odds_t / odds_r
        boot_or = []
        for b in range(B):
            pt = boot_props[b, t]; pr = boot_props[b, ref_class - 1]
            if 0 < pt < 1 and 0 < pr < 1:
                boot_or.append((pt/(1-pt)) / (pr/(1-pr)))
        if boot_or:
            ci = np.quantile(boot_or, [0.025, 0.975])
        else:
            ci = (np.nan, np.nan)
        or_rows.append({'class': t+1, 'OR_vs_ref': or_pt, 'CI_lo': ci[0], 'CI_hi': ci[1],
                        'n_valid_boot': len(boot_or)})
    return res, pd.DataFrame(or_rows)

def bch_mean(outcome_col, ref_class):
    individual_y = m[['ID', outcome_col]].dropna().set_index('ID')[outcome_col]
    individual = m[m['ID'].isin(individual_y.index)].copy()
    individual['y'] = individual_y.reindex(individual['ID']).values
    individual['modal'] = individual['class_K4'].astype(int)
    rng = np.random.default_rng(42); B = 2000
    rows = []
    boot_means = np.zeros((B, K))
    for t in range(K):
        w_t = D_inv[(individual['modal'] - 1).astype(int).values, t]
        rows.append({'target_class': t+1, 'BCH_mean': (w_t * individual['y']).sum() / w_t.sum()})
    for b in range(B):
        idx = rng.integers(0, len(individual), size=len(individual))
        sub = individual.iloc[idx]
        for t in range(K):
            w_t = D_inv[(sub['modal'] - 1).astype(int).values, t]
            boot_means[b, t] = (w_t * sub['y']).sum() / w_t.sum()
    res = pd.DataFrame(rows)
    res['SE'] = np.std(boot_means, axis=0)
    res['CI_lo'] = np.quantile(boot_means, 0.025, axis=0)
    res['CI_hi'] = np.quantile(boot_means, 0.975, axis=0)
    return res

# Determine the C2 (Stable-Low) reference class in 9-item
# From the 9-item profiles: class 1 has intercept 1.33 (lowest) -> Stable-Low
# So C2 in old labeling = C1 in 9-item
class_profiles = pd.read_csv(os.path.join(P9, 'class_profiles_K4_9item.csv'))
print('\n9-item class profiles (for ref selection):')
print(class_profiles[['class', 'proportion', 'beta_intercept', 'mean_w1', 'mean_w7']].round(3).to_string(index=False))
# class 1 = Stable-Low (intercept 1.33)
# class 2 = Persistent-High (intercept 2.61, peak 2.62)
# class 3 = Decreasing (intercept 1.93, declining)
# class 4 = Late-Increasing (intercept 1.86, rising 1.49 -> 1.93)
ref_class_9item = 1  # Stable-Low
LABEL_MAP_9 = {1: 'Stable-Low', 2: 'Persistent-High', 3: 'Decreasing', 4: 'Late-Increasing'}
print(f'\nReference class (Stable-Low): {ref_class_9item}')

print('\n=== BCH-corrected suicidal ideation (any, ≥2) ===')
res_any, or_any = bch_logit_clipped('suic_any', ref_class_9item)
res_any['class_label'] = res_any['target_class'].map(LABEL_MAP_9)
print(res_any.round(3).to_string(index=False))
print(or_any.round(3).to_string(index=False))
res_any.to_csv(os.path.join(BCH9, 'bch_suic_any_proportions.csv'), index=False, encoding='utf-8-sig')
or_any.to_csv(os.path.join(BCH9, 'bch_suic_any_OR.csv'), index=False, encoding='utf-8-sig')

print('\n=== BCH-corrected suicidal ideation (strong, ≥3) ===')
res_strong, or_strong = bch_logit_clipped('suic_strong', ref_class_9item)
res_strong['class_label'] = res_strong['target_class'].map(LABEL_MAP_9)
print(res_strong.round(3).to_string(index=False))
print(or_strong.round(3).to_string(index=False))
res_strong.to_csv(os.path.join(BCH9, 'bch_suic_strong_proportions.csv'), index=False, encoding='utf-8-sig')
or_strong.to_csv(os.path.join(BCH9, 'bch_suic_strong_OR.csv'), index=False, encoding='utf-8-sig')

print('\n=== BCH-corrected mean depression W7 ===')
res_dep = bch_mean('depression_w7', ref_class_9item)
res_dep['class_label'] = res_dep['target_class'].map(LABEL_MAP_9)
print(res_dep.round(3).to_string(index=False))
res_dep.to_csv(os.path.join(BCH9, 'bch_depression_w7_means.csv'), index=False, encoding='utf-8-sig')

print('\n=== BCH-corrected mean life sat W7 ===')
res_lsat = bch_mean('life_sat_w7', ref_class_9item)
res_lsat['class_label'] = res_lsat['target_class'].map(LABEL_MAP_9)
print(res_lsat.round(3).to_string(index=False))
res_lsat.to_csv(os.path.join(BCH9, 'bch_life_sat_w7_means.csv'), index=False, encoding='utf-8-sig')

# ===== Longitudinal SI (W1-W7) by class =====
suic_long_rows = []
for w in range(1, 8):
    df, _ = pyreadstat.read_dta(os.path.join(ROOT, f'KCYPS2018m1Yw{w}.dta'))
    s = df[f'YPSY4E04w{w}'].where((df[f'YPSY4E04w{w}']>=1) & (df[f'YPSY4E04w{w}']<=4))
    suic_long_rows.append(pd.DataFrame({
        'ID': df['ID'].values, 'wave': w,
        'suic_raw': s.values,
        'suic_any': (s >= 2).astype(float),
        'suic_strong': (s >= 3).astype(float),
    }))
suic_long = pd.concat(suic_long_rows)
suic_long = suic_long.merge(class_df[['ID', 'class_K4']], on='ID', how='left')
suic_long.to_csv(os.path.join(LD9, 'suic_long.csv'), index=False, encoding='utf-8-sig')

rates_any = suic_long.groupby(['class_K4', 'wave'])['suic_any'].mean().unstack().round(3)
rates_strong = suic_long.groupby(['class_K4', 'wave'])['suic_strong'].mean().unstack().round(3)
rates_any.to_csv(os.path.join(LD9, 'suic_any_by_class_wave.csv'), encoding='utf-8-sig')
rates_strong.to_csv(os.path.join(LD9, 'suic_strong_by_class_wave.csv'), encoding='utf-8-sig')
print('\n=== Suicidal ideation any (≥2) rates by class × wave ===')
print(rates_any.to_string())
print('\n=== Suicidal ideation strong (≥3) rates by class × wave ===')
print(rates_strong.to_string())

# ===== Cox + Schoenfeld =====
surv_rows = []
for pid, sub in suic_long.dropna(subset=['suic_strong']).groupby('ID'):
    sub = sub.sort_values('wave')
    cls = sub['class_K4'].iloc[0]
    if pd.isna(cls): continue
    events = sub[sub['suic_strong'] == 1]
    if len(events):
        first = events.iloc[0]
        surv_rows.append({'ID': pid, 'class_K4': int(cls),
                          'time': first['wave'], 'event': 1})
    else:
        last = sub.iloc[-1]
        surv_rows.append({'ID': pid, 'class_K4': int(cls),
                          'time': last['wave'], 'event': 0})
surv_df = pd.DataFrame(surv_rows)
surv_df.to_csv(os.path.join(LD9, 'survival_data.csv'), index=False, encoding='utf-8-sig')

# Cox with class indicators (ref = class 1 Stable-Low)
df_cox = surv_df.copy()
df_cox['c2'] = (df_cox['class_K4'] == 2).astype(int)  # Persistent-High
df_cox['c3'] = (df_cox['class_K4'] == 3).astype(int)  # Decreasing
df_cox['c4'] = (df_cox['class_K4'] == 4).astype(int)  # Late-Increasing
cph = CoxPHFitter()
cph.fit(df_cox[['time', 'event', 'c2', 'c3', 'c4']], duration_col='time', event_col='event')
cph_summary = cph.summary[['exp(coef)', 'exp(coef) lower 95%', 'exp(coef) upper 95%', 'p']]
cph_summary.columns = ['HR', 'HR_lower', 'HR_upper', 'p']
cph_summary = cph_summary.round(3)
cph_summary.index = ['c2 (Persistent-High vs Stable-Low)',
                     'c3 (Decreasing vs Stable-Low)',
                     'c4 (Late-Increasing vs Stable-Low)']
print('\n=== Cox HR (vs Stable-Low ref) ===')
print(cph_summary.to_string())
cph_summary.to_csv(os.path.join(LD9, 'cox_HR_class_9item.csv'), encoding='utf-8-sig')

# Schoenfeld residual test for PH assumption
print('\n=== Schoenfeld residual PH test ===')
ph_test = proportional_hazard_test(cph, df_cox[['time', 'event', 'c2', 'c3', 'c4']],
                                    time_transform='rank')
ph_summary = ph_test.summary
print(ph_summary.round(3).to_string())
ph_summary.to_csv(os.path.join(LD9, 'cox_PH_schoenfeld_test.csv'), encoding='utf-8-sig')
with open(os.path.join(LD9, 'cox_PH_test_summary.txt'), 'w', encoding='utf-8') as f:
    f.write('Schoenfeld residual proportional-hazards test\n')
    f.write(f'Cox model: time-to-first-strong-suicidal-ideation by 9-item LCGA class\n')
    f.write(f'Reference class: Stable-Low (class 1)\n\n')
    f.write(ph_summary.round(4).to_string())
    f.write('\n\nInterpretation: PH assumption is supported if all p > .05.\n')

# Class-specific event counts
print('\n=== Event rates by 9-item class ===')
for k in [1, 2, 3, 4]:
    sub = surv_df[surv_df['class_K4'] == k]
    print(f'  C{k} ({LABEL_MAP_9[k]}): N={len(sub)}, events={sub.event.sum()} ({sub.event.mean()*100:.1f}%)')

# Plot KM curves
labels = {1: 'Stable-Low (25.4%)', 2: 'Persistent-High (10.5%)',
          3: 'Decreasing (29.3%)', 4: 'Late-Increasing (34.8%)'}
colors_9 = {1: '#2ca02c', 2: '#d62728', 3: '#1f77b4', 4: '#ff7f0e'}
fig, ax = plt.subplots(figsize=(10, 6))
for k in [1, 2, 3, 4]:
    sub = surv_df[surv_df['class_K4'] == k]
    kmf = KaplanMeierFitter()
    kmf.fit(sub['time'], sub['event'], label=labels[k])
    kmf.plot_survival_function(ax=ax, color=colors_9[k], linewidth=2.5, ci_alpha=0.1)
ax.set_xlabel('Wave (1-7)')
ax.set_ylabel('Probability free of strong suicidal ideation')
ax.set_title('KM curves (9-item primary): time to first strong suicidal ideation')
ax.set_xticks(range(1, 8))
ax.legend()
ax.grid(alpha=0.3)
plt.tight_layout()
plt.savefig(os.path.join(LD9, 'KM_curves_9item.png'), dpi=300, bbox_inches='tight')
plt.close()

print('\nAll 9-item BCH+Cox outputs saved to', BCH9, 'and', LD9)
