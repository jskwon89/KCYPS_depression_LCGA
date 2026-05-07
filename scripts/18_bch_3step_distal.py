# -*- coding: utf-8 -*-
"""
BCH 3-step (Bolck, Croon, Hagenaars 2004; Vermunt 2010) corrected distal outcome analysis.

Procedure:
  Step 1. Fit LCGA on observed Y (depression) -> posterior P(C|Y)
  Step 2. Build classification error matrix D (K x K):
            D[s, t] = E[ class_pred = t | true_class = s ]
          Estimated as D[s, t] = sum_i p_post[i, s] * I(modal[i] == t) / sum_i p_post[i, s]
  Step 3. Build BCH weights for each individual:
            w_BCH[i, t] = sum_s D^{-1}[s, t] * p_post[i, s]
            (one weight column per "treating each individual as if in class t")
  Step 4. For distal Y_dist, fit weighted regression separately for each class t,
            using w_BCH[:, t] as the analytic weight, observations get treated as
            class-t with classification-error correction.

For binary outcome (suicidal ideation):
   Method: re-shape data: each individual contributes K rows, one per class,
            with weight w_BCH[i, t] and class indicator t.
   Then run logistic with class as predictor and BCH weights.

Reference:
  Bakk & Vermunt (2016). Robustness of Stepwise LCA Approaches with Continuous
    Distal Outcomes. SEM. 23(1), 20-31.
"""
import pandas as pd
import numpy as np
import os
from numpy.linalg import inv
import statsmodels.api as sm
import statsmodels.formula.api as smf
import warnings
warnings.filterwarnings('ignore')

ROOT = r'D:/2026/SCI/청소년패널조사'
OUT  = os.path.join(ROOT, 'output', 'paperA_lcga')
BCH  = os.path.join(OUT, 'bch')
os.makedirs(BCH, exist_ok=True)

# Load posterior probabilities and class assignments
class_df = pd.read_csv(os.path.join(OUT, 'class_assignments_K4.csv'), encoding='utf-8-sig')
post = class_df[[c for c in class_df.columns if c.startswith('prob_class')]].values
modal = class_df['class_K4'].values
K = post.shape[1]
N = post.shape[0]
print(f'N = {N}, K = {K}')

# --- Build classification-error matrix D ---
# D[s, t] = P(modal=t | true class=s)
# Estimate via empirical: weight individuals by post[i, s] (= P(true class=s | data))
D = np.zeros((K, K))
for s in range(K):
    denom = post[:, s].sum()
    for t in range(K):
        D[s, t] = (post[:, s] * (modal == (t+1))).sum() / denom
print('\nClassification-error matrix D (rows = true class, cols = modal):')
print(pd.DataFrame(D, index=[f's={k+1}' for k in range(K)],
                   columns=[f't={k+1}' for k in range(K)]).round(3).to_string())

# Diagonal entries should be highest if entropy is good
print('Diagonal (correct classification):', np.round(np.diag(D), 3))

# --- BCH weights (Vermunt's formula) ---
# Method: use modal class assignment as the "predicted", and back out via D^-1
# w_BCH[i, t] = sum_s D^{-1}[modal_i=t, s] * p_post[i, s] is an alternate.
# Standard BCH weights are (D^{-T})_{s, t} for individuals modal-class t
# and target latent class s.
# Concretely: each individual i has weight matrix W[i, t] given their MODAL class.
# Build W[modal_class, target_class] = D^{-1} entries.
D_inv = inv(D)
print('\nD^-1:')
print(np.round(D_inv, 3))

# For BCH 3-step: each individual contributes K rows (one per target class t).
# Weight = D_inv[modal_class - 1, t]
# Then logistic regression: y ~ class_indicator with these weights captures
# the relation between (corrected) latent class and outcome.
print('\n=== Building BCH-weighted long-form data ===')

# Load distal outcome (suicidal ideation W7)
import pyreadstat
df7, _ = pyreadstat.read_dta(os.path.join(ROOT, 'KCYPS2018m1Yw7.dta'))
df7['suic_w7_raw'] = df7['YPSY4E04w7'].where((df7['YPSY4E04w7']>=1) & (df7['YPSY4E04w7']<=4))
df7['suic_any'] = (df7['suic_w7_raw'] >= 2).astype(float)
df7['suic_strong'] = (df7['suic_w7_raw'] >= 3).astype(float)

# Get also depression score W7, life sat W7, etc.
desc = pd.read_parquet(os.path.join(ROOT, 'output', 'desc', 'long_panel.parquet'))
w7 = desc[desc.wave == 7].set_index('ID')[['depression', 'life_sat', 'self_esteem', 'happy']]
w7.columns = [c + '_w7' for c in w7.columns]
suic_df = df7.set_index('ID')[['suic_w7_raw', 'suic_any', 'suic_strong']]
distal = w7.join(suic_df, how='outer')

# Add YGENDER from W1
df1, _ = pyreadstat.read_dta(os.path.join(ROOT, 'KCYPS2018m1Yw1.dta'))
distal['female'] = (df1.set_index('ID')['YGENDERw1'] == 2).astype(float)

# Merge with class info
m = class_df.merge(distal, on='ID', how='left')
m['modal'] = m['class_K4'].astype(int)

# Build long form
long_rows = []
for _, row in m.iterrows():
    modal_i = int(row['modal']) - 1
    for t in range(K):
        long_rows.append({
            'ID': row['ID'],
            'target_class': t + 1,
            'BCH_weight': D_inv[modal_i, t],
            'suic_any': row.get('suic_any', np.nan),
            'suic_strong': row.get('suic_strong', np.nan),
            'depression_w7': row.get('depression_w7', np.nan),
            'life_sat_w7': row.get('life_sat_w7', np.nan),
            'female': row.get('female', np.nan),
        })
long_bch = pd.DataFrame(long_rows)
long_bch.to_csv(os.path.join(BCH, 'bch_long_form.csv'), index=False, encoding='utf-8-sig')
print(f'BCH long form: {long_bch.shape[0]} rows ({len(m)} individuals × {K} classes)')

# --- Distal outcome analyses ---
def bch_logit(outcome_col, ref_class=2):
    """Logistic regression with BCH weights, target_class as predictor."""
    df = long_bch.dropna(subset=[outcome_col, 'BCH_weight']).copy()
    df['y'] = df[outcome_col].astype(int)
    df['target_class'] = df['target_class'].astype(int)
    # Frequency-weighted GLM with logit link
    # Note: BCH weights can be NEGATIVE (matrix inverse) — that's intentional
    # but standard GLM software requires non-negative. Use survey-style or
    # restrict to absolute weights (interpretable as "expected counts").
    # Workaround: aggregate weighted counts and run logistic on grouped data.
    # The standard BCH approach in Mplus uses the weighted likelihood directly
    # which CAN handle negative weights via pseudo-likelihood.
    # We use the alternative: logistic with class indicators and BCH weights as
    # frequency weights, allowing negative weights. statsmodels GLM will fail
    # on neg weights, so we use survey design (via linearization) or pseudo-MLE.

    # Simple alternative: aggregate to predicted P(outcome=1 | class=t)
    # using classification-corrected estimator:
    #   P(Y=1 | class=t) = [sum_i D_inv[modal_i, t] * y_i] / [sum_i D_inv[modal_i, t]]
    # This is the BCH MEAN/PROPORTION estimator (Vermunt 2010, Eq. 7).
    print(f'\n--- BCH-corrected proportions: {outcome_col} ---')
    rows = []
    individual = df.drop_duplicates('ID').copy()
    individual_y = m[['ID', outcome_col]].dropna().set_index('ID')[outcome_col]
    individual = individual[individual['ID'].isin(individual_y.index)].copy()
    individual['y'] = individual_y.reindex(individual['ID']).values
    individual['modal'] = individual['ID'].map(m.set_index('ID')['modal'])
    # Compute corrected proportion per target class
    for t in range(K):
        # weights for individuals
        w_t = D_inv[(individual['modal'] - 1).astype(int).values, t]
        num = (w_t * individual['y']).sum()
        denom = w_t.sum()
        prop = num / denom if denom != 0 else np.nan
        rows.append({'target_class': t+1, 'BCH_corrected_proportion': prop, 'sum_weights': denom})
    res = pd.DataFrame(rows)

    # Bootstrap SE for the proportion
    rng = np.random.default_rng(42)
    B = 2000
    boot_props = np.zeros((B, K))
    for b in range(B):
        idx = rng.integers(0, len(individual), size=len(individual))
        sub = individual.iloc[idx]
        for t in range(K):
            w_t = D_inv[(sub['modal'] - 1).astype(int).values, t]
            num = (w_t * sub['y']).sum()
            denom = w_t.sum()
            boot_props[b, t] = num / denom if denom != 0 else np.nan
    res['SE'] = np.std(boot_props, axis=0)
    res['CI_lo'] = np.quantile(boot_props, 0.025, axis=0)
    res['CI_hi'] = np.quantile(boot_props, 0.975, axis=0)
    print(res.round(3).to_string(index=False))

    # Pairwise OR vs reference
    print(f'\nORs vs class {ref_class} (reference):')
    or_rows = []
    p_ref = res.iloc[ref_class - 1]['BCH_corrected_proportion']
    odds_ref = p_ref / (1 - p_ref)
    for t in range(K):
        if t + 1 == ref_class: continue
        p_t = res.iloc[t]['BCH_corrected_proportion']
        odds_t = p_t / (1 - p_t)
        or_t = odds_t / odds_ref
        # Bootstrap CI for OR
        boot_or = []
        for b in range(B):
            pt = boot_props[b, t]; pr = boot_props[b, ref_class - 1]
            if 0 < pt < 1 and 0 < pr < 1:
                boot_or.append((pt/(1-pt)) / (pr/(1-pr)))
        if boot_or:
            ci = np.quantile(boot_or, [0.025, 0.975])
        else:
            ci = (np.nan, np.nan)
        or_rows.append({'class': t+1, 'OR_vs_ref': or_t, 'CI_lo': ci[0], 'CI_hi': ci[1]})
    or_df = pd.DataFrame(or_rows)
    print(or_df.round(3).to_string(index=False))
    return res, or_df

# Run for each outcome
print('=' * 80)
res_suic, or_suic = bch_logit('suic_any', ref_class=2)
res_suic.to_csv(os.path.join(BCH, 'bch_suic_any_proportions.csv'), index=False, encoding='utf-8-sig')
or_suic.to_csv(os.path.join(BCH, 'bch_suic_any_OR.csv'), index=False, encoding='utf-8-sig')

print('=' * 80)
res_sst, or_sst = bch_logit('suic_strong', ref_class=2)
res_sst.to_csv(os.path.join(BCH, 'bch_suic_strong_proportions.csv'), index=False, encoding='utf-8-sig')
or_sst.to_csv(os.path.join(BCH, 'bch_suic_strong_OR.csv'), index=False, encoding='utf-8-sig')

# Also continuous outcomes via BCH-corrected mean
def bch_mean(outcome_col, ref_class=2):
    print(f'\n--- BCH-corrected mean: {outcome_col} ---')
    individual_y = m[['ID', outcome_col]].dropna().set_index('ID')[outcome_col]
    individual = m[m['ID'].isin(individual_y.index)].copy()
    individual['y'] = individual_y.reindex(individual['ID']).values
    individual['modal'] = individual['class_K4'].astype(int)
    rows = []
    rng = np.random.default_rng(42)
    B = 2000
    means_boot = np.zeros((B, K))
    for t in range(K):
        w_t = D_inv[(individual['modal'] - 1).astype(int).values, t]
        m_t = (w_t * individual['y']).sum() / w_t.sum()
        rows.append({'target_class': t+1, 'BCH_mean': m_t})
    for b in range(B):
        idx = rng.integers(0, len(individual), size=len(individual))
        sub = individual.iloc[idx]
        for t in range(K):
            w_t = D_inv[(sub['modal'] - 1).astype(int).values, t]
            means_boot[b, t] = (w_t * sub['y']).sum() / w_t.sum()
    res = pd.DataFrame(rows)
    res['SE'] = np.std(means_boot, axis=0)
    res['CI_lo'] = np.quantile(means_boot, 0.025, axis=0)
    res['CI_hi'] = np.quantile(means_boot, 0.975, axis=0)
    print(res.round(3).to_string(index=False))
    # Pairwise diff vs ref
    diffs = []
    for t in range(K):
        if t + 1 == ref_class: continue
        diff = res.iloc[t]['BCH_mean'] - res.iloc[ref_class - 1]['BCH_mean']
        boot_diff = means_boot[:, t] - means_boot[:, ref_class - 1]
        ci = np.quantile(boot_diff, [0.025, 0.975])
        diffs.append({'class': t+1, 'diff_vs_ref': diff, 'CI_lo': ci[0], 'CI_hi': ci[1]})
    print('\nDiff vs ref:')
    print(pd.DataFrame(diffs).round(3).to_string(index=False))
    return res, diffs

print('=' * 80)
res_d, _ = bch_mean('depression_w7', ref_class=2)
res_d.to_csv(os.path.join(BCH, 'bch_depression_w7_means.csv'), index=False, encoding='utf-8-sig')
print('=' * 80)
res_l, _ = bch_mean('life_sat_w7', ref_class=2)
res_l.to_csv(os.path.join(BCH, 'bch_life_sat_w7_means.csv'), index=False, encoding='utf-8-sig')

# Save D and D_inv
pd.DataFrame(D).to_csv(os.path.join(BCH, 'classification_error_matrix_D.csv'), index=False, encoding='utf-8-sig')
pd.DataFrame(D_inv).to_csv(os.path.join(BCH, 'classification_error_matrix_Dinv.csv'), index=False, encoding='utf-8-sig')

print('\n=== BCH 3-step done ===')
print('Outputs at:', BCH)
