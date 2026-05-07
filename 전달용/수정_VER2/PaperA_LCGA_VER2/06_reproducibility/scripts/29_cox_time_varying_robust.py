# -*- coding: utf-8 -*-
"""
Time-varying Cox with robust (sandwich/cluster) standard errors and ID clustering.

Reviewer P2: previously the class × wave interaction Cox split each individual into
multiple intervals, but the model used independent-row likelihood. Lifelines'
CoxPHFitter supports cluster_col + robust=True for sandwich SEs that account for
within-individual correlation, which is the appropriate adjustment.

Also use CoxTimeVaryingFitter for cross-check.
"""
import pandas as pd
import numpy as np
import os
from lifelines import CoxPHFitter, CoxTimeVaryingFitter, KaplanMeierFitter
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

ROOT = r'D:/2026/SCI/청소년패널조사'
P9   = os.path.join(ROOT, 'output', 'paperA_lcga', 'primary_9item')
LD9  = os.path.join(P9, 'long_distal')

surv = pd.read_csv(os.path.join(LD9, 'survival_data.csv'), encoding='utf-8-sig')
LABEL = {1: 'Stable-Low', 2: 'Persistent-High', 3: 'Decreasing', 4: 'Late-Increasing'}
COLORS = {1: '#2ca02c', 2: '#d62728', 3: '#1f77b4', 4: '#ff7f0e'}

# Build extended long-form: each row is one wave-interval at risk
sx_rows = []
for _, r in surv.iterrows():
    pid = int(r['ID']); cls = int(r['class_K4'])
    time = int(r['time']); event = int(r['event'])
    for w in range(1, time + 1):
        ev = 1 if (w == time and event == 1) else 0
        sx_rows.append({'ID': pid, 'class_K4': cls,
                        'start': w - 1, 'stop': w, 'event': ev,
                        'time_var': w})
sx = pd.DataFrame(sx_rows)
sx['c2'] = (sx['class_K4'] == 2).astype(int)
sx['c3'] = (sx['class_K4'] == 3).astype(int)
sx['c4'] = (sx['class_K4'] == 4).astype(int)
sx['c2_x_t'] = sx['c2'] * sx['time_var']
sx['c3_x_t'] = sx['c3'] * sx['time_var']
sx['c4_x_t'] = sx['c4'] * sx['time_var']
print(f'Long-form Cox data: {len(sx)} intervals from {sx.ID.nunique()} unique IDs')

# Approach A: CoxPHFitter with cluster_col + robust=True
print('\n=== Approach A: CoxPHFitter with cluster_col=ID, robust=True ===')
cph_robust = CoxPHFitter()
cph_robust.fit(
    sx[['ID', 'start', 'stop', 'event', 'c2', 'c3', 'c4',
        'c2_x_t', 'c3_x_t', 'c4_x_t']],
    duration_col='stop', event_col='event', entry_col='start',
    cluster_col='ID', robust=True,
)
s_robust = cph_robust.summary[['exp(coef)', 'exp(coef) lower 95%', 'exp(coef) upper 95%',
                                'coef', 'se(coef)', 'p']]
s_robust.columns = ['HR', 'HR_lo', 'HR_hi', 'beta', 'se_robust', 'p']
print(s_robust.round(3).to_string())
s_robust.to_csv(os.path.join(LD9, 'cox_time_varying_9item_ROBUST.csv'), encoding='utf-8-sig')

# Approach B: CoxTimeVaryingFitter
print('\n=== Approach B: CoxTimeVaryingFitter (id_col=ID) ===')
ctv = CoxTimeVaryingFitter()
ctv.fit(
    sx[['ID', 'start', 'stop', 'event', 'c2', 'c3', 'c4',
        'c2_x_t', 'c3_x_t', 'c4_x_t']],
    id_col='ID', event_col='event', start_col='start', stop_col='stop',
)
s_ctv = ctv.summary[['exp(coef)', 'exp(coef) lower 95%', 'exp(coef) upper 95%',
                     'coef', 'se(coef)', 'p']]
s_ctv.columns = ['HR', 'HR_lo', 'HR_hi', 'beta', 'se', 'p']
print(s_ctv.round(3).to_string())
s_ctv.to_csv(os.path.join(LD9, 'cox_time_varying_9item_CTV.csv'), encoding='utf-8-sig')

# Compute implied per-wave HRs using ROBUST coefficients and CIs via delta method
print('\n=== Implied HRs by wave (ROBUST model) ===')
beta_main = {'c2': cph_robust.params_['c2'], 'c3': cph_robust.params_['c3'], 'c4': cph_robust.params_['c4']}
beta_int  = {'c2': cph_robust.params_['c2_x_t'], 'c3': cph_robust.params_['c3_x_t'], 'c4': cph_robust.params_['c4_x_t']}
# Robust covariance matrix
cov_mat = cph_robust.variance_matrix_
# Indices for the 6 params
params_order = list(cph_robust.params_.index)
def idx(name): return params_order.index(name)

rows = []
for w_idx in range(1, 8):
    for cls_name, label in [('c2', 'Persistent-High'),
                             ('c3', 'Decreasing'),
                             ('c4', 'Late-Increasing')]:
        b_total = beta_main[cls_name] + beta_int[cls_name] * w_idx
        # Var(beta_main + w*beta_int) = var(b_main) + w² var(b_int) + 2w cov(b_main, b_int)
        i_m = idx(cls_name); i_i = idx(f'{cls_name}_x_t')
        var_total = cov_mat.iloc[i_m, i_m] + (w_idx**2) * cov_mat.iloc[i_i, i_i] + 2*w_idx*cov_mat.iloc[i_m, i_i]
        se_total = np.sqrt(max(var_total, 0))
        hr = np.exp(b_total)
        hr_lo = np.exp(b_total - 1.96*se_total)
        hr_hi = np.exp(b_total + 1.96*se_total)
        z = b_total / se_total if se_total > 0 else np.nan
        from scipy import stats as sps
        p_w = 2 * (1 - sps.norm.cdf(abs(z)))
        rows.append({'class': label, 'wave': w_idx,
                     'HR': round(hr, 2), 'HR_lo': round(hr_lo, 2), 'HR_hi': round(hr_hi, 2),
                     'p': round(p_w, 4)})
pred_df = pd.DataFrame(rows)
pred_wide = pred_df.pivot(index='class', columns='wave', values='HR')
pred_wide.to_csv(os.path.join(LD9, 'cox_HR_per_wave_9item_ROBUST.csv'), encoding='utf-8-sig')
print(pred_wide.to_string())

# Save full table with CI per wave
pred_df.to_csv(os.path.join(LD9, 'cox_HR_per_wave_with_CI_9item_ROBUST.csv'),
               index=False, encoding='utf-8-sig')
print('\n=== Detailed per-wave HRs with robust 95% CIs ===')
print(pred_df.to_string(index=False))

# Plot Figure 6 again with robust HRs and CI envelopes
fig, ax = plt.subplots(figsize=(9, 6))
for cls_label, color_key in [('Persistent-High', 2), ('Decreasing', 3), ('Late-Increasing', 4)]:
    sub = pred_df[pred_df['class'] == cls_label].sort_values('wave')
    ax.plot(sub['wave'], sub['HR'], '-o',
            color=COLORS[color_key], linewidth=2.5, markersize=10,
            label=f'{cls_label} vs Stable-Low')
    ax.fill_between(sub['wave'], sub['HR_lo'], sub['HR_hi'],
                    color=COLORS[color_key], alpha=0.15)
ax.axhline(1.0, color='gray', linestyle='--', alpha=0.5)
ax.set_yscale('log')
ax.set_xlabel('Wave (1-7)'); ax.set_ylabel('Hazard ratio (log scale, 95% robust CI)')
ax.set_title('Figure 6 (REVISED). Time-varying class HR for first strong SI\n'
             'Cox class × wave with cluster_col=ID, robust SE')
ax.set_xticks(range(1, 8)); ax.legend(loc='upper right'); ax.grid(alpha=0.3, which='both')
plt.tight_layout()
plt.savefig(os.path.join(P9, 'figures', 'Figure6_time_varying_HR_9item_ROBUST.png'),
            dpi=300, bbox_inches='tight')
plt.savefig(os.path.join(P9, 'figures', 'Figure6_time_varying_HR_9item_ROBUST.pdf'),
            bbox_inches='tight')
plt.close()
print('\nFigure 6 (robust) saved')

# Compare: naive vs robust SE
print('\n=== Comparison: naive vs robust SE for time-varying coefficients ===')
naive = pd.read_csv(os.path.join(LD9, 'cox_time_varying_9item.csv'))
naive_idx = naive[naive.columns[0]].tolist()
print(f'Naive (no cluster):')
print(naive.to_string(index=False))
print(f'\nRobust (cluster_col=ID):')
print(s_robust.round(3).to_string())
