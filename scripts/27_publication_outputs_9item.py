# -*- coding: utf-8 -*-
"""
Publication-quality tables and figures using 9-item PRIMARY analysis.
Fixes:
  - Table 1 uses analytic-sample N=2,452 consistently
  - BCH CIs already clipped to [0,1] in 23_bch_cox_9item.py
  - Figures use 9-item primary classes
"""
import pandas as pd
import numpy as np
import os
import pyreadstat
from scipy import stats
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

ROOT = r'D:/2026/SCI/청소년패널조사'
P9   = os.path.join(ROOT, 'output', 'paperA_lcga', 'primary_9item')
TBLS = os.path.join(P9, 'tables')
FIGS = os.path.join(P9, 'figures')
os.makedirs(TBLS, exist_ok=True); os.makedirs(FIGS, exist_ok=True)

plt.rcParams.update({'font.family': 'DejaVu Sans', 'font.size': 11,
                     'axes.titlesize': 13, 'axes.labelsize': 12,
                     'xtick.labelsize': 10, 'ytick.labelsize': 10,
                     'legend.fontsize': 10, 'savefig.dpi': 600})

# 9-item primary class labels (from class_profiles_K4_9item.csv)
# class 1 = Stable-Low, class 2 = Persistent-High, class 3 = Decreasing, class 4 = Late-Increasing
CLASS_LABEL = {1: 'C1 Stable-Low', 2: 'C2 Persistent-High',
               3: 'C3 Decreasing', 4: 'C4 Late-Increasing'}
COLORS = {1: '#2ca02c', 2: '#d62728', 3: '#1f77b4', 4: '#ff7f0e'}

class_df = pd.read_csv(os.path.join(P9, 'class_assignments_K4_9item.csv'), encoding='utf-8-sig')
class_df = class_df.rename(columns={'class_K4_9item': 'class_K4'})
long_df = pd.read_parquet(os.path.join(ROOT, 'output', 'desc', 'long_panel.parquet'))
fit_df = pd.read_csv(os.path.join(P9, 'fit_summary_9item.csv'))
profiles = pd.read_csv(os.path.join(P9, 'class_profiles_K4_9item.csv'))

w1 = long_df[long_df.wave == 1].copy()
df1, _ = pyreadstat.read_dta(os.path.join(ROOT, 'KCYPS2018m1Yw1.dta'))
df1 = df1[['ID', 'YGENDERw1']].rename(columns={'YGENDERw1': 'gender'})
w1m = w1.merge(df1, on='ID', how='left').merge(class_df[['ID', 'class_K4']], on='ID', how='left')

# =============== Table 1: Sample characteristics — analytic sample N=2,452 ===============
analytic = w1m[w1m['class_K4'].notna()].copy()
N_full = 2590  # original m1 cohort
N_analytic = len(analytic)
n_class = analytic.groupby('class_K4').size()

t1_rows = []
t1_rows.append({'Variable': 'N (Original m1 cohort) - context only', 'Overall': f'{N_full}',
                **{CLASS_LABEL[k]: '' for k in [1,2,3,4]}})
t1_rows.append({'Variable': 'N (analytic sample, ≥3 valid waves)', 'Overall': f'{N_analytic}',
                **{CLASS_LABEL[k]: f'{n_class.get(k, 0)} ({n_class.get(k, 0)/N_analytic*100:.1f}%)' for k in [1,2,3,4]}})

# Gender (within analytic)
gen_overall = (analytic.gender == 2).mean()
gen_class = analytic.groupby('class_K4').apply(lambda d: (d.gender == 2).mean())
t1_rows.append({'Variable': 'Female, %', 'Overall': f'{gen_overall*100:.1f}',
                **{CLASS_LABEL[k]: f'{gen_class.get(k, np.nan)*100:.1f}' for k in [1,2,3,4]}})

for var, label in [
    ('depression9', 'Depression W1 (9-item primary), M (SD)'),
    ('depression', 'Depression W1 (10-item, sensitivity), M (SD)'),
    ('aggression', 'Aggression W1, M (SD)'),
    ('withdrawal', 'Social withdrawal W1, M (SD)'),
    ('somatic', 'Somatic symptoms W1, M (SD)'),
    ('attention', 'Attention problems W1, M (SD)'),
    ('self_esteem', 'Self-esteem W1, M (SD)'),
    ('life_sat', 'Life satisfaction W1, M (SD)'),
    ('happy', 'Subjective happiness W1, M (SD)'),
    ('grit', 'Grit W1, M (SD)'),
    ('parent_warmth', 'Parental warmth W1, M (SD)'),
    ('parent_reject', 'Parental rejection W1, M (SD)'),
    ('parent_coerce', 'Parental coercion W1, M (SD)'),
    ('parent_chaos', 'Parental chaos W1, M (SD)'),
    ('parent_autonomy', 'Parental autonomy W1, M (SD)'),
    ('parent_structure', 'Parental structure W1, M (SD)'),
    ('peer_attach', 'Peer attachment W1, M (SD)'),
    ('teacher_rel', 'Teacher relationship W1, M (SD)'),
    ('sm_addiction', 'Smartphone dependence W1, M (SD)'),
]:
    overall = f'{analytic[var].mean():.2f} ({analytic[var].std():.2f})'
    by_class = {}
    for k in [1, 2, 3, 4]:
        sub = analytic[analytic.class_K4 == k]
        by_class[k] = f'{sub[var].mean():.2f} ({sub[var].std():.2f})'
    groups = [analytic[analytic.class_K4 == k][var].dropna() for k in [1,2,3,4]]
    if all(len(g) > 5 for g in groups):
        F, p = stats.f_oneway(*groups)
        all_vals = pd.concat(groups)
        ss_b = sum(len(g) * (g.mean() - all_vals.mean())**2 for g in groups)
        ss_t = ((all_vals - all_vals.mean())**2).sum()
        eta2 = ss_b / ss_t
        sig = '***' if p < 0.001 else ('**' if p < 0.01 else ('*' if p < 0.05 else 'ns'))
        overall_with_eta = f'{analytic[var].mean():.2f} ({analytic[var].std():.2f}) [F={F:.1f}, η²={eta2:.3f}{sig}]'
    else:
        overall_with_eta = overall
    t1_rows.append({'Variable': label, 'Overall': overall_with_eta,
                    **{CLASS_LABEL[k]: by_class[k] for k in [1,2,3,4]}})
t1 = pd.DataFrame(t1_rows)
t1.to_csv(os.path.join(TBLS, 'Table1_sample_characteristics_9item.csv'), index=False, encoding='utf-8-sig')
print('Table 1 saved (analytic-sample N=2,452 consistent)')

# =============== Table 2: Model fit indices ===============
t2 = fit_df.copy()
t2.columns = ['K', 'logL', 'BIC', 'aBIC', 'AIC', 'Entropy', 'Min %', 'Max %']
t2['Min %'] = (t2['Min %'] * 100).round(1)
t2['Max %'] = (t2['Max %'] * 100).round(1)
for c in ['logL', 'BIC', 'aBIC', 'AIC']:
    t2[c] = t2[c].round(2)
t2['Entropy'] = t2['Entropy'].round(3)
# BLRT note
t2['BLRT'] = ['—' if k == 1 else ('K=2 vs K=1: LRT_obs=2,174.7, B=30, p=.032 (lowest attainable for B=30; null max=11.5 vs obs 2,174.7)' if k == 2 else 'BIC-based; not bootstrapped') for k in t2['K']]
t2.to_csv(os.path.join(TBLS, 'Table2_model_fit_9item.csv'), index=False, encoding='utf-8-sig')
print('Table 2 saved')

# =============== Table 3: Class profiles ===============
t3 = profiles.copy()
t3['proportion'] = (t3['proportion'] * 100).round(1)
t3.columns = ['Class', 'Proportion (%)', 'Intercept (β0)', 'Linear (β1)', 'Quadratic (β2)',
              'M W1', 'M W2', 'M W3', 'M W4', 'M W5', 'M W6', 'M W7']
t3['Class Label'] = t3['Class'].map(CLASS_LABEL)
cols = ['Class', 'Class Label', 'Proportion (%)', 'Intercept (β0)', 'Linear (β1)', 'Quadratic (β2)',
        'M W1', 'M W2', 'M W3', 'M W4', 'M W5', 'M W6', 'M W7']
t3 = t3[cols]
for c in cols[3:]:
    t3[c] = t3[c].round(3)
t3.to_csv(os.path.join(TBLS, 'Table3_class_profiles_9item.csv'), index=False, encoding='utf-8-sig')
print('Table 3 saved')

# =============== Table 4: Baseline by class with eta² ===============
t4_rows = []
for var, label in [
    ('depression9', 'Depression (9-item primary)'),
    ('self_esteem', 'Self-esteem'),
    ('aggression', 'Aggression'),
    ('withdrawal', 'Social withdrawal'),
    ('somatic', 'Somatic symptoms'),
    ('attention', 'Attention problems'),
    ('life_sat', 'Life satisfaction'),
    ('happy', 'Subjective happiness'),
    ('grit', 'Grit'),
    ('peer_attach', 'Peer attachment'),
    ('teacher_rel', 'Teacher relationship'),
    ('sm_addiction', 'Smartphone dependence'),
    ('parent_warmth', 'Parental warmth'),
    ('parent_reject', 'Parental rejection'),
    ('parent_coerce', 'Parental coercion'),
    ('parent_chaos', 'Parental chaos'),
    ('parent_autonomy', 'Parental autonomy'),
    ('parent_structure', 'Parental structure'),
]:
    row = {'Construct (W1)': label}
    groups = [analytic[analytic.class_K4 == k][var].dropna() for k in [1,2,3,4]]
    for k, g in zip([1,2,3,4], groups):
        row[CLASS_LABEL[k]] = f'{g.mean():.2f} ({g.std():.2f})'
    F, p = stats.f_oneway(*groups)
    all_vals = pd.concat(groups)
    ss_b = sum(len(g) * (g.mean() - all_vals.mean())**2 for g in groups)
    ss_t = ((all_vals - all_vals.mean())**2).sum()
    eta2 = ss_b / ss_t
    row['F'] = f'{F:.1f}'
    row['p'] = '<.001' if p < 0.001 else f'{p:.3f}'
    row['η²'] = f'{eta2:.3f}'
    t4_rows.append(row)
t4 = pd.DataFrame(t4_rows)
t4.to_csv(os.path.join(TBLS, 'Table4_baseline_by_class_9item.csv'), index=False, encoding='utf-8-sig')
print('Table 4 saved')

# =============== Table 5: BCH-corrected distal outcomes ===============
bch_any  = pd.read_csv(os.path.join(P9, 'bch', 'bch_suic_any_proportions.csv'))
or_any  = pd.read_csv(os.path.join(P9, 'bch', 'bch_suic_any_OR.csv'))
bch_strong = pd.read_csv(os.path.join(P9, 'bch', 'bch_suic_strong_proportions.csv'))
or_strong = pd.read_csv(os.path.join(P9, 'bch', 'bch_suic_strong_OR.csv'))
bch_dep = pd.read_csv(os.path.join(P9, 'bch', 'bch_depression_w7_means.csv'))
bch_lsat = pd.read_csv(os.path.join(P9, 'bch', 'bch_life_sat_w7_means.csv'))

REF = 1  # Stable-Low

t5_rows = []

def add_bch_block(prop_df, or_df, header, ref=REF, fmt_pct=True):
    t5_rows.append({'Outcome (W7)': header, 'Statistic': 'BCH-corrected proportion [95% CI]; OR vs C1 Stable-Low [95% CI]'})
    for k in [1, 2, 3, 4]:
        row = prop_df.iloc[k-1]
        if k == ref:
            or_str = '1.00 (ref)'
        else:
            o = or_df[or_df['class'] == k]
            if len(o):
                ov = o.iloc[0]
                or_str = f'OR = {ov["OR_vs_ref"]:.2f} [{ov["CI_lo"]:.2f}, {ov["CI_hi"]:.2f}]'
            else:
                or_str = '-'
        if fmt_pct:
            stat = f'{row["BCH_corrected_proportion"]*100:.1f}% [{row["CI_lo"]*100:.1f}, {row["CI_hi"]*100:.1f}]    {or_str}'
        else:
            stat = f'{row["BCH_mean"]:.2f} [{row["CI_lo"]:.2f}, {row["CI_hi"]:.2f}]'
        t5_rows.append({'Outcome (W7)': f'  {CLASS_LABEL[k]}', 'Statistic': stat})

add_bch_block(bch_any, or_any, 'Suicidal ideation (any, ≥2)', fmt_pct=True)
add_bch_block(bch_strong, or_strong, 'Suicidal ideation (strong, ≥3)', fmt_pct=True)

t5_rows.append({'Outcome (W7)': 'Mean depression W7 (1-4)', 'Statistic': 'BCH-corrected mean [95% CI]'})
for k in [1, 2, 3, 4]:
    row = bch_dep.iloc[k-1]
    t5_rows.append({'Outcome (W7)': f'  {CLASS_LABEL[k]}',
                    'Statistic': f'{row["BCH_mean"]:.2f} [{row["CI_lo"]:.2f}, {row["CI_hi"]:.2f}]'})

t5_rows.append({'Outcome (W7)': 'Mean life satisfaction W7 (1-4)', 'Statistic': 'BCH-corrected mean [95% CI]'})
for k in [1, 2, 3, 4]:
    row = bch_lsat.iloc[k-1]
    t5_rows.append({'Outcome (W7)': f'  {CLASS_LABEL[k]}',
                    'Statistic': f'{row["BCH_mean"]:.2f} [{row["CI_lo"]:.2f}, {row["CI_hi"]:.2f}]'})

t5 = pd.DataFrame(t5_rows)
t5.to_csv(os.path.join(TBLS, 'Table5_BCH_distal_outcomes_9item.csv'), index=False, encoding='utf-8-sig')
print('Table 5 saved')

# =============== Table 6: Cox HR with Schoenfeld test + time-varying HR ===============
cox = pd.read_csv(os.path.join(P9, 'long_distal', 'cox_HR_class_9item.csv'))
ph = pd.read_csv(os.path.join(P9, 'long_distal', 'cox_PH_schoenfeld_test.csv'))
ph = ph.rename(columns={ph.columns[0]: 'covariate'})
tv = pd.read_csv(os.path.join(P9, 'long_distal', 'cox_HR_per_wave_9item.csv'))

# Make table 6
t6_rows = []
t6_rows.append({'Class': 'Stable-Low (ref)', 'HR_avg [95% CI]': '1.00 (ref)',
                'p': '—', 'Schoenfeld p': '—', 'HR W1 → W7': '1.00 → 1.00'})
cox_idx = ['c2 (Persistent-High vs Stable-Low)', 'c3 (Decreasing vs Stable-Low)', 'c4 (Late-Increasing vs Stable-Low)']
ph_idx = ['c2', 'c3', 'c4']
class_names = ['Persistent-High', 'Decreasing', 'Late-Increasing']
for cox_i, ph_i, name in zip(cox_idx, ph_idx, class_names):
    row = cox.set_index('Unnamed: 0').loc[cox_i] if 'Unnamed: 0' in cox.columns else cox.set_index(cox.columns[0]).loc[cox_i]
    pp = ph[ph['covariate'] == ph_i].iloc[0]
    tv_row = tv[tv['class'] == name].iloc[0]
    t6_rows.append({
        'Class': f'{name} vs Stable-Low',
        'HR_avg [95% CI]': f'{row["HR"]:.2f} [{row["HR_lower"]:.2f}, {row["HR_upper"]:.2f}]',
        'p': '<.001' if row['p'] == 0.0 else f'{row["p"]:.3f}',
        'Schoenfeld p': '<.001' if pp['p'] < 0.001 else f'{pp["p"]:.3f}',
        'HR W1 → W7': f'{tv_row["1"]:.1f} → {tv_row["7"]:.1f}',
    })
t6 = pd.DataFrame(t6_rows)
t6.to_csv(os.path.join(TBLS, 'Table6_cox_HR_9item.csv'), index=False, encoding='utf-8-sig')
print('Table 6 saved')

# =============== FIGURES ===============

# Figure 1: 9-item trajectories
y9 = long_df.pivot_table(index='ID', columns='wave', values='depression9')[list(range(1,8))]
y_with_class = y9.merge(class_df.set_index('ID')[['class_K4']], left_index=True, right_index=True)

print('Bootstrap CI for trajectories ...')
B = 500
rng = np.random.default_rng(42)
boot_means = {k: np.zeros((B, 7)) for k in [1,2,3,4]}
for b in range(B):
    if b % 100 == 0: print(f'  b={b}')
    for k in [1, 2, 3, 4]:
        sub = y_with_class[y_with_class.class_K4 == k]
        idx = rng.integers(0, len(sub), size=len(sub))
        sample = sub.iloc[idx, :7]
        boot_means[k][b] = sample.mean(axis=0).values

fig, ax = plt.subplots(figsize=(10, 7))
times = np.arange(1, 8)
prop_pct = {k: profiles.iloc[k-1]['proportion']*100 for k in [1,2,3,4]}
n_pct = {k: int((class_df['class_K4'] == k).sum()) for k in [1,2,3,4]}
for k in [1, 2, 3, 4]:
    sub = y_with_class[y_with_class.class_K4 == k]
    means = sub.iloc[:, :7].mean(axis=0).values
    lo = np.quantile(boot_means[k], 0.025, axis=0)
    hi = np.quantile(boot_means[k], 0.975, axis=0)
    ax.fill_between(times, lo, hi, color=COLORS[k], alpha=0.18)
    label = f'{CLASS_LABEL[k]} ({prop_pct[k]:.1f}%, n={n_pct[k]})'
    ax.plot(times, means, '-o', color=COLORS[k], linewidth=3, markersize=10,
            label=label, markeredgecolor='white', markeredgewidth=1.2)

ax.axhspan(2.0, 3.0, alpha=0.05, color='red', zorder=0)
ax.text(7.1, 2.5, 'Clinical\nrange', fontsize=9, color='darkred', alpha=0.7)
ax.set_xlabel('Wave (W1=2018 Grade 7 → W7=2024 ≈ 1 yr post Grade 12)', fontsize=12)
ax.set_ylabel('9-item depression score (1-4; YPSY4E excl. E04 SI)', fontsize=12)
ax.set_title('Figure 1. PRIMARY: LCGA of Korean adolescent depression trajectories\n(9-item, criterion-clean; KCYPS m1 cohort, N=2,452)', fontsize=13)
ax.set_xticks(times)
ax.set_xticklabels(['W1\n(13y)', 'W2\n(14y)', 'W3\n(15y)', 'W4\n(16y)', 'W5\n(17y)', 'W6\n(18y)', 'W7\n(19y)'])
ax.grid(alpha=0.3)
ax.legend(loc='upper right', fontsize=10, framealpha=0.95)
ax.set_ylim(1.0, 3.0)
ax.set_xlim(0.7, 7.6)
plt.tight_layout()
plt.savefig(os.path.join(FIGS, 'Figure1_LCGA_trajectories_9item.png'), dpi=600, bbox_inches='tight')
plt.close()
print('Figure 1 saved')

# Figure 2: SI by class and wave + KM
import sys
suic = pd.read_csv(os.path.join(P9, 'long_distal', 'suic_long.csv'), encoding='utf-8-sig')
rates_any = suic.groupby(['class_K4', 'wave'])['suic_any'].mean().unstack()
rates_strong = suic.groupby(['class_K4', 'wave'])['suic_strong'].mean().unstack()

from lifelines import KaplanMeierFitter
surv = pd.read_csv(os.path.join(P9, 'long_distal', 'survival_data.csv'), encoding='utf-8-sig')

fig, axes = plt.subplots(1, 3, figsize=(18, 6))
for k in [1, 2, 3, 4]:
    axes[0].plot(rates_any.columns, rates_any.loc[float(k)].values, '-o', color=COLORS[k],
                 linewidth=2.5, markersize=8, label=CLASS_LABEL[k])
    axes[1].plot(rates_strong.columns, rates_strong.loc[float(k)].values, '-o', color=COLORS[k],
                 linewidth=2.5, markersize=8, label=CLASS_LABEL[k])
axes[0].set_title('A. Any suicidal ideation (≥2)')
axes[1].set_title('B. Strong suicidal ideation (≥3)')
for ax in axes[:2]:
    ax.set_xlabel('Wave'); ax.set_ylabel('Endorsement rate')
    ax.set_xticks(range(1, 8)); ax.legend(fontsize=9, loc='upper right'); ax.grid(alpha=0.3)
axes[0].set_ylim(0, 0.85); axes[1].set_ylim(0, 0.40)

ax = axes[2]
for k in [1, 2, 3, 4]:
    sub = surv[surv['class_K4'] == k]
    kmf = KaplanMeierFitter()
    kmf.fit(sub['time'], sub['event'], label=CLASS_LABEL[k])
    kmf.plot_survival_function(ax=ax, color=COLORS[k], linewidth=2.5, ci_alpha=0.1)
ax.set_xlabel('Wave'); ax.set_ylabel('P(free of strong SI)')
ax.set_title('C. Time to first strong SI (KM)')
ax.set_xticks(range(1, 8)); ax.legend(fontsize=9, loc='lower left'); ax.grid(alpha=0.3); ax.set_ylim(0, 1.05)
plt.suptitle('Figure 2. Suicidal ideation by 9-item depression trajectory class (PRIMARY)', fontsize=14, weight='bold')
plt.tight_layout()
plt.savefig(os.path.join(FIGS, 'Figure2_suicidal_ideation_9item.png'), dpi=600, bbox_inches='tight')
plt.close()
print('Figure 2 saved')

# Figure 3: covariate heatmap
covars = ['parent_warmth', 'parent_reject', 'parent_coerce', 'parent_chaos',
          'parent_autonomy', 'parent_structure', 'peer_attach', 'teacher_rel',
          'self_esteem', 'life_sat', 'happy', 'grit', 'depression9', 'somatic',
          'aggression', 'withdrawal', 'attention', 'sm_addiction']
covar_labels = ['Parent: Warmth', 'Parent: Rejection', 'Parent: Coercion',
                'Parent: Chaos', 'Parent: Autonomy', 'Parent: Structure',
                'Peer attachment', 'Teacher relationship',
                'Self-esteem', 'Life satisfaction', 'Subjective happiness', 'Grit',
                'Depression 9-item W1', 'Somatic symptoms', 'Aggression', 'Social withdrawal',
                'Attention problems', 'Smartphone dependence']

zmat = np.zeros((len(covars), 4))
for i, v in enumerate(covars):
    vz = (analytic[v] - analytic[v].mean()) / analytic[v].std()
    for j, k in enumerate([1, 2, 3, 4]):
        zmat[i, j] = vz[analytic['class_K4'] == k].mean()

fig, ax = plt.subplots(figsize=(8, 10))
im = ax.imshow(zmat, cmap='RdBu_r', vmin=-1, vmax=1, aspect='auto')
ax.set_xticks(range(4))
ax.set_xticklabels([CLASS_LABEL[k].split(' ', 1)[1].replace(' ', '\n') for k in [1,2,3,4]], fontsize=10)
ax.set_yticks(range(len(covars)))
ax.set_yticklabels(covar_labels)
for i in range(len(covars)):
    for j in range(4):
        v = zmat[i, j]
        color = 'white' if abs(v) > 0.5 else 'black'
        ax.text(j, i, f'{v:+.2f}', ha='center', va='center', color=color, fontsize=9)
ax.set_title('Figure 3. Baseline (W1) z-score profile by 9-item trajectory class', fontsize=12)
cbar = plt.colorbar(im, ax=ax, fraction=0.04)
cbar.set_label('z-score (deviation from grand mean)')
plt.tight_layout()
plt.savefig(os.path.join(FIGS, 'Figure3_baseline_heatmap_9item.png'), dpi=600, bbox_inches='tight')
plt.close()
print('Figure 3 saved')

# Figure 4: forest plot — BCH ORs + Cox time-varying HRs
fig, ax = plt.subplots(figsize=(10, 7))
forest_rows = []
for _, r in or_any.iterrows():
    forest_rows.append({'group': 'Suic ideation any (BCH OR W7)',
                       'subgroup': f'{CLASS_LABEL[int(r["class"])]} vs C1',
                       'est': r['OR_vs_ref'], 'lo': r['CI_lo'], 'hi': r['CI_hi'],
                       'cls': int(r['class'])})
for _, r in or_strong.iterrows():
    forest_rows.append({'group': 'Suic ideation strong (BCH OR W7)',
                       'subgroup': f'{CLASS_LABEL[int(r["class"])]} vs C1',
                       'est': r['OR_vs_ref'], 'lo': r['CI_lo'], 'hi': r['CI_hi'],
                       'cls': int(r['class'])})
for cls_i, cls_name in [(2, 'Persistent-High'), (3, 'Decreasing'), (4, 'Late-Increasing')]:
    cox_row = cox.set_index(cox.columns[0]).iloc[cls_i - 2]
    forest_rows.append({'group': 'Time to strong SI (Cox HR avg)',
                       'subgroup': f'{CLASS_LABEL[cls_i]} vs C1',
                       'est': cox_row['HR'], 'lo': cox_row['HR_lower'], 'hi': cox_row['HR_upper'],
                       'cls': cls_i})

fdf = pd.DataFrame(forest_rows)
y_pos = np.arange(len(fdf))[::-1]
for i, (_, r) in enumerate(fdf.iterrows()):
    color = COLORS[r['cls']]
    ax.errorbar(r['est'], y_pos[i], xerr=[[r['est'] - r['lo']], [r['hi'] - r['est']]],
                fmt='o', color=color, markersize=8, capsize=4, linewidth=2)
ax.set_yticks(y_pos)
ax.set_yticklabels([f"{r['group']}\n  {r['subgroup']}" for _, r in fdf.iterrows()], fontsize=8)
ax.axvline(1.0, color='gray', linestyle='--', alpha=0.6)
ax.set_xscale('log')
ax.set_xlabel('Effect (OR or HR), log scale; ref = C1 Stable-Low', fontsize=11)
ax.set_title('Figure 4. BCH-corrected ORs and Cox HRs (9-item primary)', fontsize=12)
ax.grid(alpha=0.3, which='both')
plt.tight_layout()
plt.savefig(os.path.join(FIGS, 'Figure4_distal_forest_9item.png'), dpi=600, bbox_inches='tight')
plt.close()
print('Figure 4 saved')

# Figure 5: model fit indices
fig, axes = plt.subplots(1, 2, figsize=(13, 5))
ax = axes[0]
ax.plot(fit_df['K'], fit_df['BIC'], '-o', label='BIC', linewidth=2, markersize=10)
ax.plot(fit_df['K'], fit_df['aBIC'], '-s', label='aBIC', linewidth=2, markersize=8)
ax.plot(fit_df['K'], fit_df['AIC'], '-^', label='AIC', linewidth=2, markersize=8)
ax.axvline(4, color='red', linestyle='--', alpha=0.5)
ax.text(4.05, ax.get_ylim()[1]*0.95, 'K=4 selected (5% rule)', color='red', fontsize=9)
ax.set_xlabel('K'); ax.set_ylabel('Information criterion')
ax.set_title('A. Information criteria (9-item primary)')
ax.set_xticks([1,2,3,4,5]); ax.legend(); ax.grid(alpha=0.3)

ax = axes[1]
ax.plot(fit_df['K'], fit_df['entropy'], '-o', color='#d62728', linewidth=2, markersize=10)
ax.axvline(4, color='red', linestyle='--', alpha=0.5)
ax.set_xlabel('K'); ax.set_ylabel('Entropy')
ax.set_title('B. Entropy (≥0.7 = good)')
ax.set_xticks([1,2,3,4,5]); ax.set_ylim(0.6, 1.05)
ax.axhline(0.7, color='gray', linestyle=':', alpha=0.5)
ax.grid(alpha=0.3)
plt.suptitle('Figure 5. LCGA model fit indices (9-item primary)', fontsize=13, weight='bold')
plt.tight_layout()
plt.savefig(os.path.join(FIGS, 'Figure5_fit_indices_9item.png'), dpi=600, bbox_inches='tight')
plt.close()
print('Figure 5 saved')

# Figure 6: time-varying Cox HRs by wave
fig, ax = plt.subplots(figsize=(9, 6))
for cls_label, color_key in [('Persistent-High', 2), ('Decreasing', 3), ('Late-Increasing', 4)]:
    hrs = [tv[tv['class'] == cls_label].iloc[0][str(w)] for w in range(1, 8)]
    ax.plot(range(1, 8), hrs, '-o', color=COLORS[color_key], linewidth=2.5, markersize=10,
            label=f'{cls_label} vs Stable-Low')
ax.axhline(1.0, color='gray', linestyle='--', alpha=0.5)
ax.set_yscale('log')
ax.set_xlabel('Wave (1-7)'); ax.set_ylabel('Hazard ratio (log scale)')
ax.set_title('Figure 6. Time-varying class hazard ratios for first strong SI\n(Cox with class × time interaction; addresses PH violation)')
ax.set_xticks(range(1, 8)); ax.legend(loc='upper right'); ax.grid(alpha=0.3, which='both')
plt.tight_layout()
plt.savefig(os.path.join(FIGS, 'Figure6_time_varying_HR_9item.png'), dpi=600, bbox_inches='tight')
plt.close()
print('Figure 6 saved')

print('\n=== All 9-item primary tables and figures saved ===')
print('Tables:', TBLS)
print('Figures:', FIGS)
