# -*- coding: utf-8 -*-
"""
Build publication-quality tables for Paper A:
  Table 1 - Sample characteristics (overall and by class)
  Table 2 - LCGA model fit indices (K=1..5)
  Table 3 - Class profiles (intercept, slope, quadratic + estimated means by wave)
  Table 4 - Baseline (W1) covariates by class with effect sizes
  Table 5 - BCH-corrected distal outcomes at W7 with ORs (95% CI)
  Table 6 - Cox HR for time-to-strong-suicidal-ideation
"""
import pandas as pd
import numpy as np
import os
import pyreadstat
from scipy import stats

ROOT = r'D:/2026/SCI/청소년패널조사'
OUT  = os.path.join(ROOT, 'output', 'paperA_lcga')
TBLS = os.path.join(OUT, 'tables')
os.makedirs(TBLS, exist_ok=True)

class_df = pd.read_csv(os.path.join(OUT, 'class_assignments_K4.csv'), encoding='utf-8-sig')
long_df = pd.read_parquet(os.path.join(ROOT, 'output', 'desc', 'long_panel.parquet'))
fit_summary = pd.read_csv(os.path.join(OUT, 'fit_summary.csv'))
class_profiles = pd.read_csv(os.path.join(OUT, 'class_profiles_K4.csv'))

CLASS_LABELS = {
    1: 'C1 Late-Increasing',
    2: 'C2 Stable-Low',
    3: 'C3 Decreasing',
    4: 'C4 Persistent-High',
}

# ===== Table 1: Sample characteristics =====
w1 = long_df[long_df.wave==1].copy()
df1, _ = pyreadstat.read_dta(os.path.join(ROOT, 'KCYPS2018m1Yw1.dta'))
df1 = df1[['ID', 'YGENDERw1', 'ARA1Aw1', 'ARA2Aw1']].rename(columns={
    'YGENDERw1': 'gender', 'ARA1Aw1': 'region', 'ARA2Aw1': 'urban_size'
})
w1m = w1.merge(df1, on='ID', how='left')
w1m = w1m.merge(class_df[['ID', 'class_K4']], on='ID', how='left')

def make_descriptive_row(name, vals_overall, vals_by_class):
    """Format row for Table 1. Returns dict for overall and per-class."""
    row = {'Variable': name, 'Overall': vals_overall}
    for k in [1, 2, 3, 4]:
        row[CLASS_LABELS[k]] = vals_by_class.get(k, '')
    return row

t1_rows = []
n_total = len(w1m)
n_class = w1m.groupby('class_K4').size()
t1_rows.append({'Variable': 'N',
                'Overall': f'{n_total}',
                **{CLASS_LABELS[k]: f'{n_class.get(k, 0)} ({n_class.get(k, 0)/n_total*100:.1f}%)' for k in [1,2,3,4]}})

# Gender
gen_overall = (w1m.gender == 2).mean()
gen_class = w1m.groupby('class_K4').apply(lambda d: (d.gender == 2).mean())
t1_rows.append({'Variable': 'Female, %',
                'Overall': f'{gen_overall*100:.1f}',
                **{CLASS_LABELS[k]: f'{gen_class.get(k, np.nan)*100:.1f}' for k in [1,2,3,4]}})

# Continuous baseline scales — Mean (SD)
for var, label in [
    ('depression', 'Depression W1, M (SD)'),
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
    ('parent_autonomy', 'Parental autonomy support W1, M (SD)'),
    ('parent_structure', 'Parental structure W1, M (SD)'),
    ('peer_attach', 'Peer attachment W1, M (SD)'),
    ('teacher_rel', 'Teacher relationship W1, M (SD)'),
    ('sm_addiction', 'Smartphone dependence W1, M (SD)'),
]:
    overall = f'{w1m[var].mean():.2f} ({w1m[var].std():.2f})'
    by_class = {}
    for k in [1, 2, 3, 4]:
        sub = w1m[w1m.class_K4 == k]
        by_class[k] = f'{sub[var].mean():.2f} ({sub[var].std():.2f})'
    # Eta squared (effect size for ANOVA across classes)
    groups = [w1m[w1m.class_K4 == k][var].dropna() for k in [1,2,3,4]]
    if all(len(g) > 5 for g in groups):
        F, p = stats.f_oneway(*groups)
        # Eta-squared
        all_vals = pd.concat(groups)
        ss_between = sum(len(g) * (g.mean() - all_vals.mean())**2 for g in groups)
        ss_total = ((all_vals - all_vals.mean())**2).sum()
        eta2 = ss_between / ss_total
        sig = '***' if p < 0.001 else ('**' if p < 0.01 else ('*' if p < 0.05 else 'ns'))
        overall_with_eta = f'{w1m[var].mean():.2f} ({w1m[var].std():.2f}) [F={F:.1f}, η²={eta2:.3f}{sig}]'
    else:
        overall_with_eta = overall
    t1_rows.append({'Variable': label,
                    'Overall': overall_with_eta,
                    **{CLASS_LABELS[k]: by_class[k] for k in [1,2,3,4]}})

t1 = pd.DataFrame(t1_rows)
t1.to_csv(os.path.join(TBLS, 'Table1_sample_characteristics.csv'), index=False, encoding='utf-8-sig')
print('Table 1 created')

# ===== Table 2: Model fit indices =====
# Add BLRT p-values if available
blrt_path = os.path.join(OUT, 'blrt', 'blrt_summary.csv')
blrt = None
if os.path.exists(blrt_path):
    blrt = pd.read_csv(blrt_path)
    blrt = blrt.set_index('K_alt')

t2 = fit_summary.copy()
t2.columns = ['K', 'logL', 'BIC', 'aBIC', 'AIC', 'Entropy', 'Min %', 'Max %']
t2['Min %'] = (t2['Min %'] * 100).round(1)
t2['Max %'] = (t2['Max %'] * 100).round(1)
t2['logL'] = t2['logL'].round(2)
t2['BIC'] = t2['BIC'].round(2)
t2['aBIC'] = t2['aBIC'].round(2)
t2['AIC'] = t2['AIC'].round(2)
t2['Entropy'] = t2['Entropy'].round(3)
if blrt is not None:
    t2['BLRT p'] = t2['K'].map(lambda k: blrt['p_BLRT'].get(k, np.nan))
    t2['BLRT p'] = t2['BLRT p'].apply(lambda p: '<.001' if pd.notna(p) and p < 0.001 else (f'{p:.3f}' if pd.notna(p) else '-'))
t2.to_csv(os.path.join(TBLS, 'Table2_model_fit.csv'), index=False, encoding='utf-8-sig')
print('Table 2 created')
print(t2.to_string(index=False))

# ===== Table 3: Class profiles =====
t3 = class_profiles.copy()
t3['proportion'] = (t3['proportion'] * 100).round(1)
t3.columns = ['Class', 'Proportion (%)', 'Intercept (β0)', 'Linear slope (β1)', 'Quadratic (β2)',
              'Mean W1', 'Mean W2', 'Mean W3', 'Mean W4', 'Mean W5', 'Mean W6', 'Mean W7']
t3['Class Label'] = t3['Class'].map(CLASS_LABELS)
cols = ['Class', 'Class Label', 'Proportion (%)', 'Intercept (β0)', 'Linear slope (β1)', 'Quadratic (β2)',
        'Mean W1', 'Mean W2', 'Mean W3', 'Mean W4', 'Mean W5', 'Mean W6', 'Mean W7']
t3 = t3[cols]
for c in cols[3:]:
    t3[c] = t3[c].round(3)
t3.to_csv(os.path.join(TBLS, 'Table3_class_profiles.csv'), index=False, encoding='utf-8-sig')
print('Table 3 created')

# ===== Table 4: Baseline covariates by class with eta² =====
# This is essentially a transposed/cleaner version of Table 1
t4_rows = []
for var, label in [
    ('depression', 'Depression'),
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
    ('parent_autonomy', 'Parental autonomy support'),
    ('parent_structure', 'Parental structure'),
]:
    row = {'Construct (W1)': label}
    groups = [w1m[w1m.class_K4 == k][var].dropna() for k in [1,2,3,4]]
    for k, g in zip([1,2,3,4], groups):
        row[f'{CLASS_LABELS[k]}'] = f'{g.mean():.2f} ({g.std():.2f})'
    # F and eta²
    F, p = stats.f_oneway(*groups)
    all_vals = pd.concat(groups)
    ss_b = sum(len(g) * (g.mean() - all_vals.mean())**2 for g in groups)
    ss_t = ((all_vals - all_vals.mean())**2).sum()
    eta2 = ss_b / ss_t
    p_str = '<.001' if p < 0.001 else f'{p:.3f}'
    row['F'] = f'{F:.1f}'
    row['p'] = p_str
    row['η²'] = f'{eta2:.3f}'
    t4_rows.append(row)
t4 = pd.DataFrame(t4_rows)
t4.to_csv(os.path.join(TBLS, 'Table4_baseline_by_class.csv'), index=False, encoding='utf-8-sig')
print('Table 4 created')

# ===== Table 5: BCH-corrected distal outcomes =====
bch_props_any = pd.read_csv(os.path.join(OUT, 'bch', 'bch_suic_any_proportions.csv'))
bch_or_any = pd.read_csv(os.path.join(OUT, 'bch', 'bch_suic_any_OR.csv'))
bch_props_strong = pd.read_csv(os.path.join(OUT, 'bch', 'bch_suic_strong_proportions.csv'))
bch_or_strong = pd.read_csv(os.path.join(OUT, 'bch', 'bch_suic_strong_OR.csv'))
bch_dep_means = pd.read_csv(os.path.join(OUT, 'bch', 'bch_depression_w7_means.csv'))
bch_lsat_means = pd.read_csv(os.path.join(OUT, 'bch', 'bch_life_sat_w7_means.csv'))

t5_rows = []
# Suicidal ideation (any) by class
t5_rows.append({'Outcome (W7)': 'Suicidal ideation (any, ≥2)',
                'Statistic': 'BCH-corrected proportion [95% CI]'})
for k in [1, 2, 3, 4]:
    row = bch_props_any.iloc[k-1]
    or_row = bch_or_any[bch_or_any['class'] == k]
    or_str = ''
    if len(or_row):
        orv = or_row.iloc[0]
        or_str = f'OR = {orv["OR_vs_ref"]:.2f} [{orv["CI_lo"]:.2f}, {orv["CI_hi"]:.2f}]'
    elif k == 2:
        or_str = '1.00 (ref)'
    p = row['BCH_corrected_proportion']
    ci = (row['CI_lo'], row['CI_hi'])
    t5_rows.append({'Outcome (W7)': f'  {CLASS_LABELS[k]}',
                    'Statistic': f'{p*100:.1f}% [{ci[0]*100:.1f}, {ci[1]*100:.1f}]    {or_str}'})

t5_rows.append({'Outcome (W7)': 'Suicidal ideation (strong, ≥3)',
                'Statistic': 'BCH-corrected proportion [95% CI]'})
for k in [1, 2, 3, 4]:
    row = bch_props_strong.iloc[k-1]
    or_row = bch_or_strong[bch_or_strong['class'] == k]
    or_str = ''
    if len(or_row):
        orv = or_row.iloc[0]
        or_str = f'OR = {orv["OR_vs_ref"]:.2f} [{orv["CI_lo"]:.2f}, {orv["CI_hi"]:.2f}]'
    elif k == 2:
        or_str = '1.00 (ref)'
    p = row['BCH_corrected_proportion']
    ci = (row['CI_lo'], row['CI_hi'])
    t5_rows.append({'Outcome (W7)': f'  {CLASS_LABELS[k]}',
                    'Statistic': f'{p*100:.1f}% [{ci[0]*100:.1f}, {ci[1]*100:.1f}]    {or_str}'})

t5_rows.append({'Outcome (W7)': 'Depression score (1-4)',
                'Statistic': 'BCH-corrected mean [95% CI]'})
for k in [1, 2, 3, 4]:
    row = bch_dep_means.iloc[k-1]
    t5_rows.append({'Outcome (W7)': f'  {CLASS_LABELS[k]}',
                    'Statistic': f'{row["BCH_mean"]:.2f} [{row["CI_lo"]:.2f}, {row["CI_hi"]:.2f}]'})

t5_rows.append({'Outcome (W7)': 'Life satisfaction (1-4)',
                'Statistic': 'BCH-corrected mean [95% CI]'})
for k in [1, 2, 3, 4]:
    row = bch_lsat_means.iloc[k-1]
    t5_rows.append({'Outcome (W7)': f'  {CLASS_LABELS[k]}',
                    'Statistic': f'{row["BCH_mean"]:.2f} [{row["CI_lo"]:.2f}, {row["CI_hi"]:.2f}]'})

t5 = pd.DataFrame(t5_rows)
t5.to_csv(os.path.join(TBLS, 'Table5_BCH_distal_outcomes.csv'), index=False, encoding='utf-8-sig')
print('Table 5 created')

# ===== Table 6: Cox HR for time-to-strong-SI =====
cox = pd.read_csv(os.path.join(OUT, 'long_distal', 'cox_HR_class.csv'))
cox.columns = ['Class', 'HR', 'HR_lower', 'HR_upper', 'p']
cox['Class'] = cox['Class'].map({'c1': 'C1 Late-Increasing',
                                  'c3': 'C3 Decreasing',
                                  'c4': 'C4 Persistent-High'})
# Add reference
ref = pd.DataFrame([{'Class': 'C2 Stable-Low (ref)', 'HR': 1.00, 'HR_lower': '-', 'HR_upper': '-', 'p': '-'}])
cox = pd.concat([ref, cox], ignore_index=True)
cox['HR (95% CI)'] = cox.apply(
    lambda r: f'{r["HR"]:.2f} [{r["HR_lower"]:.2f}, {r["HR_upper"]:.2f}]' if r['HR_lower'] != '-' else '1.00 (ref)',
    axis=1
)
t6 = cox[['Class', 'HR (95% CI)', 'p']].copy()
t6['p'] = t6['p'].apply(lambda p: '<.001' if (p == 0.0 or p == '<.001') else p)
t6.to_csv(os.path.join(TBLS, 'Table6_cox_HR.csv'), index=False, encoding='utf-8-sig')
print('Table 6 created')

# Print summary
print('\n=== Tables 1-6 generated in', TBLS, '===')
print('\nTable 2 — Model fit:')
print(t2.to_string(index=False))
print('\nTable 5 — BCH-corrected distal outcomes:')
print(t5.to_string(index=False))
print('\nTable 6 — Cox HR for time to first strong suicidal ideation:')
print(t6.to_string(index=False))
