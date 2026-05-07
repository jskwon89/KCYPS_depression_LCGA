# -*- coding: utf-8 -*-
"""
Longitudinal distal outcomes: suicidal ideation (E04) across W4-W7 by class.
Plus: time-to-first-strong-ideation survival analysis.
"""
import pandas as pd
import numpy as np
import os
import pyreadstat
from lifelines import KaplanMeierFitter, CoxPHFitter
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

ROOT = r'D:/2026/SCI/청소년패널조사'
OUT  = os.path.join(ROOT, 'output', 'paperA_lcga')
LDIST = os.path.join(OUT, 'long_distal')
os.makedirs(LDIST, exist_ok=True)

class_df = pd.read_csv(os.path.join(OUT, 'class_assignments_K4.csv'), encoding='utf-8-sig')

# Pull suicidal ideation E04 from each wave 1..7
rows = []
for w in range(1, 8):
    df, _ = pyreadstat.read_dta(os.path.join(ROOT, f'KCYPS2018m1Yw{w}.dta'))
    col = f'YPSY4E04w{w}'
    if col not in df.columns: continue
    s = df[col].where((df[col]>=1) & (df[col]<=4))
    rows.append(pd.DataFrame({
        'ID': df['ID'].values,
        'wave': w,
        'suic_raw': s.values,
        'suic_any': (s >= 2).astype(float),
        'suic_strong': (s >= 3).astype(float),
    }))
suic_long = pd.concat(rows)
suic_long = suic_long.merge(class_df[['ID', 'class_K4']], on='ID', how='left')
suic_long.to_csv(os.path.join(LDIST, 'suic_long.csv'), index=False, encoding='utf-8-sig')

# Wave-specific endorsement rates by class
print('=== Suicidal ideation (any, ≥2) by class × wave ===')
rates_any = suic_long.groupby(['class_K4', 'wave'])['suic_any'].mean().unstack().round(3)
rates_any.to_csv(os.path.join(LDIST, 'suic_any_by_class_wave.csv'), encoding='utf-8-sig')
print(rates_any.to_string())

print('\n=== Suicidal ideation (strong, ≥3) by class × wave ===')
rates_strong = suic_long.groupby(['class_K4', 'wave'])['suic_strong'].mean().unstack().round(3)
rates_strong.to_csv(os.path.join(LDIST, 'suic_strong_by_class_wave.csv'), encoding='utf-8-sig')
print(rates_strong.to_string())

# Plot: suicidal ideation rate by class across waves
fig, axes = plt.subplots(1, 2, figsize=(14, 5))
colors = {1: '#ff7f0e', 2: '#2ca02c', 3: '#1f77b4', 4: '#d62728'}
labels = {1: 'C1 Late-Increasing (36.3%)', 2: 'C2 Stable-Low (25.6%)',
          3: 'C3 Decreasing (27.2%)', 4: 'C4 Persistent-High (10.9%)'}
for k in [1, 2, 3, 4]:
    axes[0].plot(rates_any.columns, rates_any.loc[k].values, '-o', color=colors[k],
                 linewidth=2.5, markersize=9, label=labels[k])
    axes[1].plot(rates_strong.columns, rates_strong.loc[k].values, '-o', color=colors[k],
                 linewidth=2.5, markersize=9, label=labels[k])
axes[0].set_title('Any suicidal ideation (≥2)', fontsize=12)
axes[1].set_title('Strong suicidal ideation (≥3, "often")', fontsize=12)
for ax in axes:
    ax.set_xlabel('Wave (W1=2018 Grade 7 → W7=2024)', fontsize=11)
    ax.set_ylabel('Endorsement rate', fontsize=11)
    ax.legend(loc='upper right', fontsize=9)
    ax.grid(alpha=0.3)
plt.tight_layout()
plt.savefig(os.path.join(LDIST, 'suic_by_class_across_waves.png'), dpi=300, bbox_inches='tight')
plt.close()

# Time-to-first-strong-suicidal-ideation survival
# Each individual: event = first wave where suic_strong=1; if never, censored at last observed wave
print('\n=== Time-to-first STRONG suicidal ideation ===')
surv_rows = []
for pid, sub in suic_long.dropna(subset=['suic_strong']).groupby('ID'):
    sub = sub.sort_values('wave')
    cls = sub['class_K4'].iloc[0] if pd.notna(sub['class_K4'].iloc[0]) else None
    if cls is None: continue
    events = sub[sub['suic_strong'] == 1]
    if len(events):
        first_event = events.iloc[0]
        surv_rows.append({'ID': pid, 'class_K4': int(cls),
                          'time': first_event['wave'], 'event': 1})
    else:
        last_obs = sub.iloc[-1]
        surv_rows.append({'ID': pid, 'class_K4': int(cls),
                          'time': last_obs['wave'], 'event': 0})
surv_df = pd.DataFrame(surv_rows)
surv_df.to_csv(os.path.join(LDIST, 'survival_data.csv'), index=False, encoding='utf-8-sig')

print(f'\nSurvival sample: N={len(surv_df)}, events={surv_df.event.sum()} '
      f'({surv_df.event.mean()*100:.1f}% experienced ≥1 strong suicidal ideation)')

# KM curves by class
fig, ax = plt.subplots(figsize=(9, 6))
kmfs = {}
for k in [1, 2, 3, 4]:
    sub = surv_df[surv_df['class_K4'] == k]
    kmf = KaplanMeierFitter()
    kmf.fit(sub['time'], sub['event'], label=labels[k])
    kmfs[k] = kmf
    kmf.plot_survival_function(ax=ax, color=colors[k], linewidth=2.5, ci_show=True, ci_alpha=0.1)
ax.set_xlabel('Wave (1-7)', fontsize=12)
ax.set_ylabel('Survival prob (no strong suicidal ideation)', fontsize=12)
ax.set_title('Time to first strong suicidal ideation by depression trajectory class', fontsize=13)
ax.set_xticks(range(1, 8))
ax.grid(alpha=0.3)
plt.tight_layout()
plt.savefig(os.path.join(LDIST, 'KM_curves_by_class.png'), dpi=300, bbox_inches='tight')
plt.close()

# Cox proportional hazards: class as predictor (vs ref C2)
print('\n=== Cox proportional hazards (class vs C2 ref) ===')
df_cox = surv_df.copy()
df_cox['c1'] = (df_cox['class_K4'] == 1).astype(int)
df_cox['c3'] = (df_cox['class_K4'] == 3).astype(int)
df_cox['c4'] = (df_cox['class_K4'] == 4).astype(int)
cph = CoxPHFitter()
cph.fit(df_cox[['time', 'event', 'c1', 'c3', 'c4']], duration_col='time', event_col='event')
cph_summary = cph.summary[['exp(coef)', 'exp(coef) lower 95%', 'exp(coef) upper 95%', 'p']]
cph_summary.columns = ['HR', 'HR_lower', 'HR_upper', 'p']
cph_summary = cph_summary.round(3)
print(cph_summary.to_string())
cph_summary.to_csv(os.path.join(LDIST, 'cox_HR_class.csv'), encoding='utf-8-sig')

# Median wave to event by class
print('\n=== Median time to event (KM) ===')
for k in [1, 2, 3, 4]:
    sub = surv_df[surv_df['class_K4'] == k]
    n_event = sub.event.sum()
    p_event = sub.event.mean()
    print(f'  C{k}: N={len(sub)}, events={n_event} ({p_event*100:.1f}%)')

print('\nDone. Longitudinal distal results in', LDIST)
