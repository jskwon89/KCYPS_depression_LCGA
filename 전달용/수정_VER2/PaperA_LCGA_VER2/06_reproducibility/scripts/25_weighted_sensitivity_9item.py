# -*- coding: utf-8 -*-
"""
Weighted LCGA sensitivity for 9-item depression.
Fixes: Hungarian label alignment + correct Cohen kappa save.
"""
import pandas as pd
import numpy as np
import os
import pyreadstat
from scipy.special import logsumexp
from scipy.optimize import linear_sum_assignment
from sklearn.metrics import cohen_kappa_score

ROOT = r'D:/2026/SCI/청소년패널조사'
DESC = os.path.join(ROOT, 'output', 'desc')
P9   = os.path.join(ROOT, 'output', 'paperA_lcga', 'primary_9item')
SENS = os.path.join(P9, 'sensitivity_weighted')
os.makedirs(SENS, exist_ok=True)

long_df = pd.read_parquet(os.path.join(DESC, 'long_panel.parquet'))
y = long_df.pivot_table(index='ID', columns='wave', values='depression9')[list(range(1,8))]
ids_full = y.index.values
Y_full = y.values

df7, _ = pyreadstat.read_dta(os.path.join(ROOT, 'KCYPS2018m1Yw7.dta'))
w_col = next((c for c in df7.columns if c.startswith('WEIGHTB2w7')), None)
weights_w7 = df7.set_index('ID')[w_col].reindex(ids_full)
mean_w = np.nanmean(weights_w7)
weights_norm = weights_w7.fillna(mean_w).values / mean_w

valid_count = np.sum(~np.isnan(Y_full), axis=1)
keep = valid_count >= 3
Y = Y_full[keep]; ids = ids_full[keep]; w = weights_norm[keep]
N = Y.shape[0]
print(f'N (>=3 waves): {N}')
print(f'Weight: median={np.median(w):.3f}, range [{w.min():.3f}, {w.max():.3f}]')

times = np.arange(1, 8, dtype=float) - 4
X = np.column_stack([np.ones(7), times, times**2])

def fit_lcga_weighted(Y, K, X, weights, max_iter=300, n_starts=8, seed=42):
    N, T = Y.shape; P = X.shape[1]
    rng = np.random.default_rng(seed)
    w_norm = np.asarray(weights, dtype=float)
    w_norm = w_norm * (N / w_norm.sum())
    best = None
    for start in range(n_starts):
        if start == 0:
            mean_y = np.nanmean(Y, axis=1)
            cuts = np.nanquantile(mean_y, np.linspace(0, 1, K+1))
            init = np.zeros(N, dtype=int)
            for k in range(K):
                if k == K-1: init[mean_y >= cuts[k]] = k
                else:        init[(mean_y >= cuts[k]) & (mean_y < cuts[k+1])] = k
        else:
            init = rng.integers(0, K, size=N)
        beta = np.zeros((K, P))
        for k in range(K):
            mask = init == k
            if mask.sum() < 5:
                beta[k] = [np.nanmean(Y), 0, 0]; continue
            xs, ys = [], []
            for i in np.where(mask)[0]:
                obs = ~np.isnan(Y[i])
                xs.append(X[obs]); ys.append(Y[i][obs])
            beta[k] = np.linalg.lstsq(np.vstack(xs), np.concatenate(ys), rcond=None)[0]
        sigma2 = np.nanvar(Y); pi = np.full(K, 1/K)
        prev_ll = -np.inf
        for it in range(max_iter):
            mu = np.einsum('kp,tp->kt', beta, X)
            log_lik_ik = np.full((N, K), -np.inf)
            for i in range(N):
                obs = ~np.isnan(Y[i])
                if not obs.any(): continue
                resid2 = (Y[i][obs][None, :] - mu[:, obs])**2
                log_lik_ik[i] = -0.5*obs.sum()*np.log(2*np.pi*sigma2) - 0.5*resid2.sum(axis=1)/sigma2
            log_post_unnorm = log_lik_ik + np.log(pi + 1e-300)
            log_z = logsumexp(log_post_unnorm, axis=1)
            post = np.exp(log_post_unnorm - log_z[:, None])
            ll_w = (w_norm * log_z).sum()
            if abs(ll_w - prev_ll) < 1e-5: break
            prev_ll = ll_w
            wpost = post * w_norm[:, None]
            pi = wpost.sum(axis=0) / wpost.sum()
            pi = np.clip(pi, 1e-6, None); pi /= pi.sum()
            for k in range(K):
                XtWX = np.zeros((P, P)); XtWy = np.zeros(P)
                for i in range(N):
                    obs = ~np.isnan(Y[i])
                    if not obs.any(): continue
                    Xi = X[obs]; ww = wpost[i, k]
                    XtWX += ww * Xi.T @ Xi; XtWy += ww * Xi.T @ Y[i][obs]
                if np.linalg.matrix_rank(XtWX) < P: continue
                beta[k] = np.linalg.solve(XtWX, XtWy)
            ss_total = w_total = 0.0
            mu = np.einsum('kp,tp->kt', beta, X)
            for i in range(N):
                obs = ~np.isnan(Y[i])
                if not obs.any(): continue
                resid2 = (Y[i][obs][None, :] - mu[:, obs])**2
                ss_total += (wpost[i] * resid2.sum(axis=1)).sum()
                w_total  += (wpost[i] * obs.sum()).sum()
            sigma2 = ss_total / max(w_total, 1)
        if best is None or ll_w > best['ll']:
            best = dict(K=K, ll=ll_w, pi=pi.copy(), beta=beta.copy(), sigma2=sigma2,
                        post=post.copy(), assign=np.argmax(post, axis=1))
    return best

print('\n=== Weighted LCGA K=4 (9-item) ===')
m4 = fit_lcga_weighted(Y, 4, X, w, n_starts=8, seed=42)
print(f'  proportions: {m4["pi"].round(3)}')

# Save profiles
profiles = pd.DataFrame({
    'class': range(1, 5),
    'proportion_weighted': m4['pi'],
    'beta_intercept': m4['beta'][:, 0],
    'beta_linear': m4['beta'][:, 1],
    'beta_quadratic': m4['beta'][:, 2],
})
mu = np.einsum('kp,tp->kt', m4['beta'], X)
for w_idx in range(7):
    profiles[f'mean_w{w_idx+1}'] = mu[:, w_idx]
profiles.to_csv(os.path.join(SENS, 'weighted_class_profiles_K4_9item.csv'), index=False, encoding='utf-8-sig')
print('\nWeighted class profiles:')
print(profiles.round(3).to_string(index=False))

# Compare to unweighted 9-item
unw = pd.read_csv(os.path.join(P9, 'class_assignments_K4_9item.csv'))
unw_aligned = unw.set_index('ID')['class_K4_9item'].reindex(ids).values

# Crosstab BEFORE alignment
ct_before = pd.crosstab(pd.Series(unw_aligned, name='unweighted'),
                         pd.Series(m4['assign'] + 1, name='weighted'),
                         dropna=False, margins=True)
ct_before.to_csv(os.path.join(SENS, 'crosstab_unweighted_vs_weighted_BEFORE_align.csv'), encoding='utf-8-sig')
kappa_raw = cohen_kappa_score(unw_aligned, m4['assign'] + 1)

# Hungarian alignment: maximize agreement
ct_arr = ct_before.iloc[:-1, :-1].values
row_ind, col_ind = linear_sum_assignment(-ct_arr)
mapping = dict(zip(col_ind + 1, row_ind + 1))  # weighted label -> unweighted label
weighted_aligned = np.array([mapping[c] for c in (m4['assign'] + 1)])

ct_after = pd.crosstab(pd.Series(unw_aligned, name='unweighted'),
                        pd.Series(weighted_aligned, name='weighted_aligned'),
                        margins=True)
ct_after.to_csv(os.path.join(SENS, 'crosstab_unweighted_vs_weighted_AFTER_align.csv'), encoding='utf-8-sig')

agreement = (unw_aligned == weighted_aligned).mean()
kappa_aligned = cohen_kappa_score(unw_aligned, weighted_aligned)

print(f'\n=== Crosstab BEFORE Hungarian alignment ===')
print(ct_before.to_string())
print(f'\nRaw (label-naive) Cohen kappa: {kappa_raw:.4f}')
print(f'\nHungarian mapping (weighted label -> unweighted label): {mapping}')
print(f'\n=== Crosstab AFTER alignment ===')
print(ct_after.to_string())
print(f'\nPost-alignment percent agreement: {agreement:.4f}')
print(f'Post-alignment Cohen kappa: {kappa_aligned:.4f}')

with open(os.path.join(SENS, 'kappa_aligned.txt'), 'w', encoding='utf-8') as f:
    f.write(f'Weighted vs Unweighted LCGA (9-item) — Cohen kappa\n')
    f.write(f'='*60 + '\n')
    f.write(f'Weight column: WEIGHTB2w7 (KCYPS standardized longitudinal weight)\n')
    f.write(f'N: {N}\n')
    f.write(f'Raw (label-naive) kappa: {kappa_raw:.4f}  [NOT RECOMMENDED — different label permutations]\n')
    f.write(f'Hungarian label alignment (weighted -> unweighted): {mapping}\n')
    f.write(f'Post-alignment percent agreement: {agreement:.4f}\n')
    f.write(f'Post-alignment Cohen kappa: {kappa_aligned:.4f}\n')
    f.write('\nINTERPRETATION: post-alignment kappa reflects whether the 4-class structure '
            'is robust to differential attrition. Agreement >= .90 indicates strong robustness.\n')

# Save aligned assignments
out = pd.DataFrame({
    'ID': ids,
    'class_K4_unweighted': unw_aligned,
    'class_K4_weighted_aligned': weighted_aligned,
})
out.to_csv(os.path.join(SENS, 'aligned_assignments.csv'), index=False, encoding='utf-8-sig')

print('\nDone. Aligned weighted sensitivity saved to', SENS)
