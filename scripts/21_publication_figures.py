# -*- coding: utf-8 -*-
"""
Publication-quality figures (300+ dpi) for Paper A.
  Figure 1 — 4-class depression trajectories with 95% bootstrap CI envelopes
  Figure 2 — Suicidal ideation trajectory and Cox-KM curve by class
  Figure 3 — Baseline covariate profile heatmap (z-score by class)
  Figure 4 — BCH-corrected distal outcomes forest plot
"""
import pandas as pd
import numpy as np
import os
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.patches import Patch
import warnings
warnings.filterwarnings('ignore')

ROOT = r'D:/2026/SCI/청소년패널조사'
OUT = os.path.join(ROOT, 'output', 'paperA_lcga')
FIGS = os.path.join(OUT, 'figures')
os.makedirs(FIGS, exist_ok=True)

# Set publication style
plt.rcParams.update({
    'font.family': 'DejaVu Sans',
    'font.size': 11,
    'axes.titlesize': 13,
    'axes.labelsize': 12,
    'xtick.labelsize': 10,
    'ytick.labelsize': 10,
    'legend.fontsize': 10,
    'savefig.dpi': 300,
    'savefig.bbox': 'tight',
})

CLASS_LABELS = {
    1: 'C1 Late-Increasing\n(36.3%, n=911)',
    2: 'C2 Stable-Low\n(25.6%, n=617)',
    3: 'C3 Decreasing\n(27.2%, n=668)',
    4: 'C4 Persistent-High\n(10.9%, n=256)',
}
COLORS = {1: '#ff7f0e', 2: '#2ca02c', 3: '#1f77b4', 4: '#d62728'}

class_df = pd.read_csv(os.path.join(OUT, 'class_assignments_K4.csv'), encoding='utf-8-sig')
long_df = pd.read_parquet(os.path.join(ROOT, 'output', 'desc', 'long_panel.parquet'))
profiles = pd.read_csv(os.path.join(OUT, 'class_profiles_K4.csv'))

# ===== Figure 1: depression trajectories with bootstrap CIs =====
y = long_df.pivot_table(index='ID', columns='wave', values='depression')[list(range(1,8))]
y_with_class = y.merge(class_df.set_index('ID')[['class_K4']], left_index=True, right_index=True)

# Bootstrap CI for class trajectories
print('Computing bootstrap CIs ...')
B = 500
rng = np.random.default_rng(42)
boot_means = {k: np.zeros((B, 7)) for k in [1,2,3,4]}
for b in range(B):
    if b % 100 == 0: print(f'  b={b}/{B}')
    for k in [1, 2, 3, 4]:
        sub = y_with_class[y_with_class.class_K4 == k]
        idx = rng.integers(0, len(sub), size=len(sub))
        sample = sub.iloc[idx, :7]
        boot_means[k][b] = sample.mean(axis=0).values

fig, ax = plt.subplots(figsize=(10, 7))
times = np.arange(1, 8)
for k in [1, 2, 3, 4]:
    sub = y_with_class[y_with_class.class_K4 == k]
    means = sub.iloc[:, :7].mean(axis=0).values
    lo = np.quantile(boot_means[k], 0.025, axis=0)
    hi = np.quantile(boot_means[k], 0.975, axis=0)
    ax.fill_between(times, lo, hi, color=COLORS[k], alpha=0.18)
    ax.plot(times, means, '-o', color=COLORS[k], linewidth=3, markersize=10,
            label=CLASS_LABELS[k], markeredgecolor='white', markeredgewidth=1.2)

# Annotate suicide-risk bands
ax.axhspan(2.0, 3.0, alpha=0.05, color='red', zorder=0)
ax.text(7.1, 2.5, 'Clinical\nrange', fontsize=9, color='darkred', ha='left', va='center', alpha=0.7)

ax.set_xlabel('Wave (W1=2018 Grade 7 → W7=2024 ≈ Grade 12+1)', fontsize=12)
ax.set_ylabel('Depression score (mean of 10 items, 1-4)', fontsize=12)
ax.set_title('Figure 1. Latent Class Growth Analysis of Korean adolescent depression\ntrajectories from age 13 to 19 (KCYPS-2018 m1 cohort, N=2,452)', fontsize=13)
ax.set_xticks(times)
ax.set_xticklabels(['W1\n(13y)', 'W2\n(14y)', 'W3\n(15y)', 'W4\n(16y)', 'W5\n(17y)', 'W6\n(18y)', 'W7\n(19y)'])
ax.grid(alpha=0.3)
ax.legend(loc='upper right', fontsize=10, framealpha=0.95)
ax.set_ylim(1.0, 3.0)
ax.set_xlim(0.7, 7.6)
plt.tight_layout()
plt.savefig(os.path.join(FIGS, 'Figure1_LCGA_trajectories.png'), dpi=300, bbox_inches='tight')
plt.savefig(os.path.join(FIGS, 'Figure1_LCGA_trajectories.pdf'), bbox_inches='tight')
plt.close()
print('Figure 1 saved')

# ===== Figure 2: suicidal ideation trajectory + KM curve =====
suic = pd.read_csv(os.path.join(OUT, 'long_distal', 'suic_long.csv'), encoding='utf-8-sig')
rates_any = suic.groupby(['class_K4', 'wave'])['suic_any'].mean().unstack()
rates_strong = suic.groupby(['class_K4', 'wave'])['suic_strong'].mean().unstack()

from lifelines import KaplanMeierFitter
surv = pd.read_csv(os.path.join(OUT, 'long_distal', 'survival_data.csv'), encoding='utf-8-sig')

fig, axes = plt.subplots(1, 3, figsize=(18, 6))

ax = axes[0]
for k in [1, 2, 3, 4]:
    ax.plot(rates_any.columns, rates_any.loc[float(k)].values, '-o', color=COLORS[k],
            linewidth=2.5, markersize=8, label=CLASS_LABELS[k].replace('\n', ' '))
ax.set_xlabel('Wave', fontsize=12)
ax.set_ylabel('Endorsement rate', fontsize=12)
ax.set_title('A. Any suicidal ideation (≥2 = "sometimes")', fontsize=12)
ax.set_xticks(range(1, 8))
ax.legend(fontsize=8, loc='upper right')
ax.grid(alpha=0.3)
ax.set_ylim(0, 0.85)

ax = axes[1]
for k in [1, 2, 3, 4]:
    ax.plot(rates_strong.columns, rates_strong.loc[float(k)].values, '-o', color=COLORS[k],
            linewidth=2.5, markersize=8, label=CLASS_LABELS[k].replace('\n', ' '))
ax.set_xlabel('Wave', fontsize=12)
ax.set_ylabel('Endorsement rate', fontsize=12)
ax.set_title('B. Strong suicidal ideation (≥3 = "often")', fontsize=12)
ax.set_xticks(range(1, 8))
ax.legend(fontsize=8, loc='upper right')
ax.grid(alpha=0.3)
ax.set_ylim(0, 0.40)

ax = axes[2]
for k in [1, 2, 3, 4]:
    sub = surv[surv['class_K4'] == k]
    kmf = KaplanMeierFitter()
    kmf.fit(sub['time'], sub['event'], label=CLASS_LABELS[k].replace('\n', ' '))
    kmf.plot_survival_function(ax=ax, color=COLORS[k], linewidth=2.5, ci_alpha=0.1)
ax.set_xlabel('Wave', fontsize=12)
ax.set_ylabel('Probability free of strong suicidal ideation', fontsize=12)
ax.set_title('C. Time to first strong suicidal ideation\n(Kaplan-Meier estimator)', fontsize=12)
ax.set_xticks(range(1, 8))
ax.legend(fontsize=8, loc='lower left')
ax.grid(alpha=0.3)
ax.set_ylim(0, 1.05)

plt.suptitle('Figure 2. Suicidal ideation by depression trajectory class', fontsize=14, weight='bold')
plt.tight_layout()
plt.savefig(os.path.join(FIGS, 'Figure2_suicidal_ideation_by_class.png'), dpi=300, bbox_inches='tight')
plt.savefig(os.path.join(FIGS, 'Figure2_suicidal_ideation_by_class.pdf'), bbox_inches='tight')
plt.close()
print('Figure 2 saved')

# ===== Figure 3: baseline covariate heatmap (z-scores) =====
w1 = long_df[long_df.wave == 1].set_index('ID')
covars = ['parent_warmth', 'parent_reject', 'parent_coerce', 'parent_chaos',
          'parent_autonomy', 'parent_structure', 'peer_attach', 'teacher_rel',
          'self_esteem', 'life_sat', 'happy', 'grit', 'depression', 'somatic',
          'aggression', 'withdrawal', 'attention', 'sm_addiction']
covar_labels = ['Parent: Warmth', 'Parent: Rejection', 'Parent: Coercion',
                'Parent: Chaos', 'Parent: Autonomy', 'Parent: Structure',
                'Peer attachment', 'Teacher relationship',
                'Self-esteem', 'Life satisfaction', 'Subjective happiness', 'Grit',
                'Depression W1', 'Somatic symptoms', 'Aggression', 'Social withdrawal',
                'Attention problems', 'Smartphone dependence']

w1m = w1.merge(class_df.set_index('ID')[['class_K4']], left_index=True, right_index=True)

# Standardize each covariate, then mean by class
zmat = np.zeros((len(covars), 4))
for i, v in enumerate(covars):
    vz = (w1m[v] - w1m[v].mean()) / w1m[v].std()
    for j, k in enumerate([1, 2, 3, 4]):
        zmat[i, j] = vz[w1m['class_K4'] == k].mean()

fig, ax = plt.subplots(figsize=(8, 10))
im = ax.imshow(zmat, cmap='RdBu_r', vmin=-1, vmax=1, aspect='auto')
ax.set_xticks(range(4))
ax.set_xticklabels(['C1\nLate-Inc', 'C2\nStable-Low', 'C3\nDecreasing', 'C4\nPersistent-High'])
ax.set_yticks(range(len(covars)))
ax.set_yticklabels(covar_labels)
for i in range(len(covars)):
    for j in range(4):
        v = zmat[i, j]
        color = 'white' if abs(v) > 0.5 else 'black'
        ax.text(j, i, f'{v:+.2f}', ha='center', va='center', color=color, fontsize=9)
ax.set_title('Figure 3. Baseline (W1) covariate profile by depression trajectory\n(values are z-scores: cell mean − overall mean / SD)', fontsize=12)
cbar = plt.colorbar(im, ax=ax, fraction=0.04)
cbar.set_label('z-score (deviation from grand mean)')
plt.tight_layout()
plt.savefig(os.path.join(FIGS, 'Figure3_baseline_covariate_heatmap.png'), dpi=300, bbox_inches='tight')
plt.savefig(os.path.join(FIGS, 'Figure3_baseline_covariate_heatmap.pdf'), bbox_inches='tight')
plt.close()
print('Figure 3 saved')

# ===== Figure 4: BCH-corrected distal forest plot =====
or_any = pd.read_csv(os.path.join(OUT, 'bch', 'bch_suic_any_OR.csv'))
or_strong = pd.read_csv(os.path.join(OUT, 'bch', 'bch_suic_strong_OR.csv'))
cox = pd.read_csv(os.path.join(OUT, 'long_distal', 'cox_HR_class.csv'))

fig, ax = plt.subplots(figsize=(10, 7))

# Build forest plot data
forest_rows = []
# Suicidal ideation any (BCH-corrected ORs)
for _, r in or_any.iterrows():
    forest_rows.append({
        'group': f'Suicidal ideation any (BCH OR)',
        'subgroup': f'C{int(r["class"])} vs C2',
        'est': r['OR_vs_ref'], 'lo': r['CI_lo'], 'hi': r['CI_hi']
    })
# Suicidal ideation strong
for _, r in or_strong.iterrows():
    forest_rows.append({
        'group': f'Suicidal ideation strong (BCH OR)',
        'subgroup': f'C{int(r["class"])} vs C2',
        'est': r['OR_vs_ref'], 'lo': r['CI_lo'], 'hi': r['CI_hi']
    })
# Cox HR
cox_clean = cox.rename(columns={'covariate': 'cov'})
for _, r in cox_clean.iterrows():
    cls = {'c1': 'C1', 'c3': 'C3', 'c4': 'C4'}[r['cov']]
    forest_rows.append({
        'group': 'Time to strong SI (Cox HR)',
        'subgroup': f'{cls} vs C2',
        'est': r['exp(coef)'] if 'exp(coef)' in r else r.get('HR'),
        'lo': r['exp(coef) lower 95%'] if 'exp(coef) lower 95%' in r else r.get('HR_lower'),
        'hi': r['exp(coef) upper 95%'] if 'exp(coef) upper 95%' in r else r.get('HR_upper'),
    })

fdf = pd.DataFrame(forest_rows)
y_pos = np.arange(len(fdf))[::-1]

# Color by class
sub_colors = []
for s in fdf['subgroup']:
    if 'C1' in s: sub_colors.append(COLORS[1])
    elif 'C3' in s: sub_colors.append(COLORS[3])
    elif 'C4' in s: sub_colors.append(COLORS[4])
    else: sub_colors.append('gray')

for i, (_, r) in enumerate(fdf.iterrows()):
    ypos = y_pos[i]
    ax.errorbar(r['est'], ypos, xerr=[[r['est'] - r['lo']], [r['hi'] - r['est']]],
                fmt='o', color=sub_colors[i], markersize=8, capsize=4, linewidth=2)

ax.set_yticks(y_pos)
ax.set_yticklabels([f"{r['group']}\n  {r['subgroup']}" for _, r in fdf.iterrows()],
                   fontsize=9)
ax.axvline(1.0, color='gray', linestyle='--', alpha=0.6)
ax.set_xscale('log')
ax.set_xlabel('Effect (OR or HR), log scale; reference = C2 Stable-Low', fontsize=11)
ax.set_title('Figure 4. BCH-corrected distal outcome ORs and Cox HRs by trajectory class', fontsize=12)
ax.grid(alpha=0.3, which='both')
plt.tight_layout()
plt.savefig(os.path.join(FIGS, 'Figure4_distal_forest.png'), dpi=300, bbox_inches='tight')
plt.savefig(os.path.join(FIGS, 'Figure4_distal_forest.pdf'), bbox_inches='tight')
plt.close()
print('Figure 4 saved')

# ===== Figure 5: Model fit indices visualization =====
fit = pd.read_csv(os.path.join(OUT, 'fit_summary.csv'))
fig, axes = plt.subplots(1, 2, figsize=(13, 5))

ax = axes[0]
ax.plot(fit['K'], fit['BIC'], '-o', color='#1f77b4', linewidth=2, markersize=10, label='BIC')
ax.plot(fit['K'], fit['aBIC'], '-s', color='#ff7f0e', linewidth=2, markersize=8, label='aBIC')
ax.plot(fit['K'], fit['AIC'], '-^', color='#2ca02c', linewidth=2, markersize=8, label='AIC')
# Highlight K=4
ax.axvline(4, color='red', linestyle='--', alpha=0.5)
ax.text(4.1, ax.get_ylim()[1]*0.95, 'K=4 selected', color='red', fontsize=10, va='top')
ax.set_xlabel('Number of latent classes (K)')
ax.set_ylabel('Information criterion')
ax.set_title('A. BIC, aBIC, AIC across K')
ax.set_xticks([1, 2, 3, 4, 5])
ax.legend()
ax.grid(alpha=0.3)

ax = axes[1]
ax.plot(fit['K'], fit['entropy'], '-o', color='#d62728', linewidth=2, markersize=10)
ax.axvline(4, color='red', linestyle='--', alpha=0.5)
ax.set_xlabel('Number of latent classes (K)')
ax.set_ylabel('Entropy')
ax.set_title('B. Classification entropy')
ax.set_xticks([1, 2, 3, 4, 5])
ax.set_ylim(0.6, 1.05)
ax.axhline(0.7, color='gray', linestyle=':', alpha=0.5)
ax.text(1.1, 0.71, 'Entropy ≥ 0.7 (good)', fontsize=9, color='gray')
ax.grid(alpha=0.3)

plt.suptitle('Figure 5. Model fit indices for K=1 to 5', fontsize=13, weight='bold')
plt.tight_layout()
plt.savefig(os.path.join(FIGS, 'Figure5_fit_indices.png'), dpi=300, bbox_inches='tight')
plt.savefig(os.path.join(FIGS, 'Figure5_fit_indices.pdf'), bbox_inches='tight')
plt.close()
print('Figure 5 saved')

print('\n=== All publication figures saved to', FIGS, '===')
