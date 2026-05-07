# -*- coding: utf-8 -*-
"""
Cluster bootstrap for time-varying Cox per-wave HR with PROPER robust CIs.

Reviewer P1: lifelines `variance_matrix_` is the model-based (Hessian) covariance,
not the cluster-robust sandwich. So delta-method per-wave CIs computed from
variance_matrix_ are NOT cluster-robust regardless of `robust=True` flag at fit.

Solution: ID-level cluster bootstrap (B=200). Resample individuals with replacement,
rebuild long-form intervals, refit Cox with class × time interaction, compute
per-wave HR. Take percentile 2.5/97.5 for CI; bootstrap p = fraction of bootstrap
HRs ≤ 1 (vs. > 1) for one-sided test, two-sided p = 2·min(p_left, p_right).
"""
import pandas as pd
import numpy as np
import os, time
from lifelines import CoxPHFitter
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import warnings
warnings.filterwarnings('ignore')

ROOT = r'D:/2026/SCI/청소년패널조사'
P9   = os.path.join(ROOT, 'output', 'paperA_lcga', 'primary_9item')
LD9  = os.path.join(P9, 'long_distal')

surv = pd.read_csv(os.path.join(LD9, 'survival_data.csv'), encoding='utf-8-sig')
print(f'Survival data: {len(surv)} individuals, {surv.event.sum()} events')

# Build long-form once for the original data
def build_long(s):
    rows = []
    for _, r in s.iterrows():
        pid = int(r['ID']); cls = int(r['class_K4'])
        t = int(r['time']); ev = int(r['event'])
        for w in range(1, t + 1):
            e = 1 if (w == t and ev == 1) else 0
            rows.append({'ID': pid, 'class_K4': cls,
                         'start': w - 1, 'stop': w, 'event': e,
                         'time_var': w})
    df = pd.DataFrame(rows)
    df['c2'] = (df['class_K4'] == 2).astype(int)
    df['c3'] = (df['class_K4'] == 3).astype(int)
    df['c4'] = (df['class_K4'] == 4).astype(int)
    df['c2_x_t'] = df['c2'] * df['time_var']
    df['c3_x_t'] = df['c3'] * df['time_var']
    df['c4_x_t'] = df['c4'] * df['time_var']
    return df

# Original fit (point estimates)
sx_orig = build_long(surv)
print(f'Long-form: {len(sx_orig)} intervals')
cph0 = CoxPHFitter(penalizer=1e-6)
cph0.fit(sx_orig[['start', 'stop', 'event', 'c2', 'c3', 'c4',
                   'c2_x_t', 'c3_x_t', 'c4_x_t']],
         duration_col='stop', event_col='event', entry_col='start')
beta0 = cph0.params_.to_dict()
print('Point estimates (beta):')
for k, v in beta0.items(): print(f'  {k}: {v:.4f}')

def per_wave_hr(beta_dict):
    """Compute per-wave HR for c2, c3, c4 vs ref."""
    out = {}
    for cls in ['c2', 'c3', 'c4']:
        out[cls] = {}
        for w in range(1, 8):
            b = beta_dict[cls] + beta_dict[f'{cls}_x_t'] * w
            out[cls][w] = np.exp(b)
    return out

hr_point = per_wave_hr(beta0)
print('\nPoint per-wave HRs:')
for cls in ['c2', 'c3', 'c4']:
    print(f'  {cls}: ' + ' '.join(f'{hr_point[cls][w]:.2f}' for w in range(1, 8)))

# Cluster bootstrap
B = 300
print(f'\n=== Cluster bootstrap B={B} ===')
unique_ids = surv['ID'].unique()
N_id = len(unique_ids)
print(f'  N unique IDs: {N_id}')
rng = np.random.default_rng(42)
boot_betas = []
fail = 0
t0 = time.time()
for b in range(B):
    boot_ids = rng.choice(unique_ids, size=N_id, replace=True)
    boot_surv = surv.set_index('ID').loc[boot_ids].reset_index()
    # Build long-form for bootstrap sample
    sx_b = build_long(boot_surv)
    try:
        cph_b = CoxPHFitter(penalizer=1e-6)
        cph_b.fit(sx_b[['start', 'stop', 'event', 'c2', 'c3', 'c4',
                         'c2_x_t', 'c3_x_t', 'c4_x_t']],
                   duration_col='stop', event_col='event', entry_col='start')
        boot_betas.append(cph_b.params_.to_dict())
    except Exception as e:
        fail += 1
    if (b + 1) % 25 == 0:
        elapsed = time.time() - t0
        eta = elapsed * B / (b + 1) - elapsed
        print(f'  b={b+1}/{B}, fail={fail}, elapsed={elapsed:.0f}s, ETA={eta:.0f}s')

print(f'\nBootstrap done: {len(boot_betas)} successful, {fail} failed')

# Compute per-wave HR distribution
boot_hr = {cls: {w: [] for w in range(1, 8)} for cls in ['c2', 'c3', 'c4']}
for b in boot_betas:
    hb = per_wave_hr(b)
    for cls in ['c2', 'c3', 'c4']:
        for w in range(1, 8):
            boot_hr[cls][w].append(hb[cls][w])

# Summarize
LABEL = {'c2': 'Persistent-High', 'c3': 'Decreasing', 'c4': 'Late-Increasing'}
rows = []
for cls in ['c2', 'c3', 'c4']:
    for w in range(1, 8):
        hr_vals = np.array(boot_hr[cls][w])
        hr_pt = hr_point[cls][w]
        ci_lo = np.quantile(hr_vals, 0.025)
        ci_hi = np.quantile(hr_vals, 0.975)
        # Two-sided percentile bootstrap p for HR=1
        p_left = (hr_vals <= 1).mean()
        p_right = (hr_vals >= 1).mean()
        p_two = 2 * min(p_left, p_right)
        p_two = max(p_two, 1 / (len(hr_vals) + 1))  # add-one floor
        rows.append({'class': LABEL[cls], 'wave': w, 'class_short': cls,
                     'HR': round(hr_pt, 2),
                     'HR_lo': round(ci_lo, 2), 'HR_hi': round(ci_hi, 2),
                     'boot_p': round(p_two, 4)})

pred_df = pd.DataFrame(rows)
pred_df.to_csv(os.path.join(LD9, 'cox_HR_per_wave_with_CI_9item_BOOT.csv'),
               index=False, encoding='utf-8-sig')
print('\n=== Per-wave HRs with cluster-bootstrap 95% percentile CIs ===')
print(pred_df.to_string(index=False))

# Wide
pred_wide = pred_df.pivot(index='class', columns='wave', values='HR')
pred_wide.to_csv(os.path.join(LD9, 'cox_HR_per_wave_9item_BOOT.csv'), encoding='utf-8-sig')
print('\nHR by wave (point):')
print(pred_wide.to_string())

# Bootstrap interaction-coefficient p
print('\n=== Bootstrap interaction coefficient distribution ===')
boot_int = pd.DataFrame([
    {'c2_x_t': b['c2_x_t'], 'c3_x_t': b['c3_x_t'], 'c4_x_t': b['c4_x_t']}
    for b in boot_betas
])
for cls in ['c2_x_t', 'c3_x_t', 'c4_x_t']:
    pt = beta0[cls]
    vals = boot_int[cls].values
    ci = np.quantile(vals, [0.025, 0.975])
    p_left = (vals <= 0).mean()
    p_right = (vals >= 0).mean()
    p_two = 2 * min(p_left, p_right)
    p_two = max(p_two, 1 / (len(vals) + 1))
    print(f'  {cls}: beta = {pt:.4f}, CI = [{ci[0]:.4f}, {ci[1]:.4f}], boot_p = {p_two:.4f}')
boot_int.to_csv(os.path.join(LD9, 'cox_interaction_bootstrap_dist_9item.csv'),
                index=False, encoding='utf-8-sig')

# Plot Figure 6 BOOT
COLORS = {'Persistent-High': '#d62728', 'Decreasing': '#1f77b4', 'Late-Increasing': '#ff7f0e'}
fig, ax = plt.subplots(figsize=(9, 6))
for cls_label in ['Persistent-High', 'Decreasing', 'Late-Increasing']:
    sub = pred_df[pred_df['class'] == cls_label].sort_values('wave')
    ax.plot(sub['wave'], sub['HR'], '-o',
            color=COLORS[cls_label], linewidth=2.5, markersize=10,
            label=f'{cls_label} vs Stable-Low')
    ax.fill_between(sub['wave'], sub['HR_lo'], sub['HR_hi'],
                    color=COLORS[cls_label], alpha=0.18)
ax.axhline(1.0, color='gray', linestyle='--', alpha=0.5)
ax.set_yscale('log')
ax.set_xlabel('Wave (1-7)'); ax.set_ylabel('Hazard ratio (log scale, 95% cluster-bootstrap CI)')
ax.set_title('Figure 6 (BOOTSTRAP). Time-varying class HR for first strong SI\n'
             'Cox class × wave; ID-level cluster bootstrap (B=300) percentile 95% CI')
ax.set_xticks(range(1, 8)); ax.legend(loc='upper right')
ax.grid(alpha=0.3, which='both')
plt.tight_layout()
plt.savefig(os.path.join(P9, 'figures', 'Figure6_time_varying_HR_9item_BOOT.png'),
            dpi=600, bbox_inches='tight')
plt.close()
print(f'\nSaved Figure 6 BOOT to figures/')
print('Cluster bootstrap done.')
