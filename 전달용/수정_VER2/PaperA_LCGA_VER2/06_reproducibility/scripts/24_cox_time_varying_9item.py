# -*- coding: utf-8 -*-
"""
Address PH violation for C3 (Decreasing) and C4 (Late-Increasing) in 9-item Cox.
Approach 1: Cox with time-varying class effects (time-by-class interaction)
Approach 2: Stratified Cox (treat class as strata) - not interpretable
Approach 3: Piecewise-constant HRs by wave window

We use approach 1 (time-varying coefficient) and approach 3 (piecewise HR).
"""
import pandas as pd
import numpy as np
import os
from lifelines import CoxPHFitter, KaplanMeierFitter
from lifelines.statistics import logrank_test
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

ROOT = r'D:/2026/SCI/청소년패널조사'
P9   = os.path.join(ROOT, 'output', 'paperA_lcga', 'primary_9item')
LD9  = os.path.join(P9, 'long_distal')

surv = pd.read_csv(os.path.join(LD9, 'survival_data.csv'), encoding='utf-8-sig')
LABEL = {1: 'Stable-Low', 2: 'Persistent-High', 3: 'Decreasing', 4: 'Late-Increasing'}

# Approach 3 — piecewise-constant HR by wave window
# Split: W1-W3 (early adolescence, age 13-15) vs W4-W7 (mid-late adolescence)
# Use lifelines.utils.to_long_format with cut points
from lifelines.utils import to_episodic_format

# Build long-form
surv_long_rows = []
for _, r in surv.iterrows():
    pid = r['ID']; cls = int(r['class_K4']); time = int(r['time']); event = int(r['event'])
    cuts = [3]
    intervals = []
    prev = 0
    for c in cuts:
        if time >= c:
            intervals.append((prev, c, 0))
            prev = c
    intervals.append((prev, time, event))
    for start, stop, ev in intervals:
        surv_long_rows.append({'ID': pid, 'class_K4': cls,
                               'start': start, 'stop': stop, 'event': ev,
                               'period': 'early' if stop <= 3 else 'late'})
surv_long = pd.DataFrame(surv_long_rows)
print(f'Long-form survival: {len(surv_long)} rows from {surv_long.ID.nunique()} ids')

# Fit period-stratified Cox: HR(class | period)
print('\n=== Piecewise Cox: early (W1-W3) vs late (W4-W7) ===')
results = []
for period in ['early', 'late']:
    sub = surv_long[surv_long['period'] == period].copy()
    sub['c2'] = (sub['class_K4'] == 2).astype(int)
    sub['c3'] = (sub['class_K4'] == 3).astype(int)
    sub['c4'] = (sub['class_K4'] == 4).astype(int)
    cph = CoxPHFitter()
    try:
        cph.fit(sub[['start', 'stop', 'event', 'c2', 'c3', 'c4']],
                duration_col='stop', event_col='event', entry_col='start')
        s = cph.summary[['exp(coef)', 'exp(coef) lower 95%', 'exp(coef) upper 95%', 'p']]
        s.columns = ['HR', 'HR_lo', 'HR_hi', 'p']
        s['period'] = period
        s['class_label'] = [LABEL[2], LABEL[3], LABEL[4]]
        s['class'] = ['c2', 'c3', 'c4']
        results.append(s.reset_index(drop=True))
        print(f'\n  Period: {period}')
        print(s.round(3).to_string())
    except Exception as e:
        print(f'  fit failed for period {period}: {e}')

if results:
    pw = pd.concat(results, ignore_index=True)
    pw.to_csv(os.path.join(LD9, 'cox_HR_piecewise_9item.csv'), index=False, encoding='utf-8-sig')

# Approach 1 — Cox with time-varying coefficient via tt() (lifelines doesn't have native tt() but we can fake via long-format with time-by-class interaction)
print('\n=== Cox with time-varying class effect (time × class interaction) ===')
# Build extended long-form: each row is (start, stop, event), with time-by-class interaction
# We'll create one row per individual per wave they were at risk
surv_extended = []
for _, r in surv.iterrows():
    pid = r['ID']; cls = int(r['class_K4']); time = int(r['time']); event = int(r['event'])
    for w in range(1, time + 1):
        ev = 1 if (w == time and event == 1) else 0
        surv_extended.append({'ID': pid, 'class_K4': cls,
                              'start': w - 1, 'stop': w, 'event': ev,
                              'time_var': w})
sx = pd.DataFrame(surv_extended)
sx['c2'] = (sx['class_K4'] == 2).astype(int)
sx['c3'] = (sx['class_K4'] == 3).astype(int)
sx['c4'] = (sx['class_K4'] == 4).astype(int)
sx['c2_x_t'] = sx['c2'] * sx['time_var']
sx['c3_x_t'] = sx['c3'] * sx['time_var']
sx['c4_x_t'] = sx['c4'] * sx['time_var']

cph_tv = CoxPHFitter()
try:
    cph_tv.fit(sx[['start', 'stop', 'event', 'c2', 'c3', 'c4',
                    'c2_x_t', 'c3_x_t', 'c4_x_t']],
                duration_col='stop', event_col='event', entry_col='start')
    s_tv = cph_tv.summary[['exp(coef)', 'exp(coef) lower 95%', 'exp(coef) upper 95%', 'coef', 'se(coef)', 'p']]
    s_tv.columns = ['HR', 'HR_lo', 'HR_hi', 'beta', 'se', 'p']
    s_tv = s_tv.round(3)
    s_tv.to_csv(os.path.join(LD9, 'cox_time_varying_9item.csv'), encoding='utf-8-sig')
    print(s_tv.to_string())

    # Predicted HR per wave for each class
    print('\nImplied HRs by wave (from time-varying model):')
    rows = []
    for w in range(1, 8):
        for cls_name, cls_label in [('c2', 'Persistent-High'), ('c3', 'Decreasing'), ('c4', 'Late-Increasing')]:
            beta_main = cph_tv.params_[cls_name]
            beta_int = cph_tv.params_[f'{cls_name}_x_t']
            beta_tot = beta_main + beta_int * w
            hr = np.exp(beta_tot)
            rows.append({'class': cls_label, 'wave': w, 'HR_at_wave': round(hr, 2)})
    pred_df = pd.DataFrame(rows)
    pred_df_wide = pred_df.pivot(index='class', columns='wave', values='HR_at_wave')
    pred_df_wide.to_csv(os.path.join(LD9, 'cox_HR_per_wave_9item.csv'), encoding='utf-8-sig')
    print(pred_df_wide.to_string())

    # Plot: implied HR by wave
    fig, ax = plt.subplots(figsize=(9, 6))
    colors = {'Persistent-High': '#d62728', 'Decreasing': '#1f77b4', 'Late-Increasing': '#ff7f0e'}
    for cls_label in ['Persistent-High', 'Decreasing', 'Late-Increasing']:
        hrs = [pred_df_wide.loc[cls_label, w] for w in range(1, 8)]
        ax.plot(range(1, 8), hrs, '-o', color=colors[cls_label], linewidth=2.5, markersize=10,
                label=f'{cls_label} vs Stable-Low')
    ax.axhline(1.0, color='gray', linestyle='--', alpha=0.5)
    ax.set_yscale('log')
    ax.set_xlabel('Wave')
    ax.set_ylabel('Hazard ratio for first strong suicidal ideation (log)')
    ax.set_title('Time-varying class hazard ratios: 9-item LCGA primary\nC4 (Late-Increasing) crosses C2 (Persistent-High) ~W4-5')
    ax.set_xticks(range(1, 8))
    ax.legend(loc='upper right')
    ax.grid(alpha=0.3)
    plt.tight_layout()
    plt.savefig(os.path.join(LD9, 'cox_HR_per_wave_9item.png'), dpi=300, bbox_inches='tight')
    plt.close()

except Exception as e:
    print(f'TV Cox failed: {e}')

print('\nDone — time-varying Cox results in', LD9)
