# -*- coding: utf-8 -*-
"""
Sensitivity analysis: re-fit LCGA with KCYPS longitudinal weights.

KCYPS provides per-wave longitudinal weights:
  WEIGHTB1w{w}: 종단면 원 가중치
  WEIGHTB2w{w}: 종단면 표준화 가중치 (표준화: mean=1)

Approach:
  - Pull the standardized longitudinal weight WEIGHTB2 from each wave (or W7 as primary)
  - Apply weighted likelihood: weighted EM where each individual's likelihood
    is multiplied by their weight in the M-step.
  - Compare class structure, proportions, and class-defining intercepts/slopes.
"""
import pandas as pd
import numpy as np
import os
import pyreadstat
from scipy.special import logsumexp

ROOT = r'D:/2026/SCI/청소년패널조사'
DESC = os.path.join(ROOT, 'output', 'desc')
OUT = os.path.join(ROOT, 'output', 'paperA_lcga')
SENS_OUT = os.path.join(OUT, 'sensitivity_weighted')
os.makedirs(SENS_OUT, exist_ok=True)

long_df = pd.read_parquet(os.path.join(DESC, 'long_panel.parquet'))
y = long_df.pivot_table(index='ID', columns='wave', values='depression')[list(range(1,8))]
ids_full = y.index.values
Y_full = y.values

# Get W7 longitudinal weight (most encompassing — survey through to W7)
df7, _ = pyreadstat.read_dta(os.path.join(ROOT, 'KCYPS2018m1Yw7.dta'))
w_col = next((c for c in df7.columns if c.startswith('WEIGHTB2w7')), None)
if w_col is None:
    # Fallback to WEIGHTB1
    w_col = next((c for c in df7.columns if c.startswith('WEIGHTB1w7')), None)
print('Using weight column:', w_col)
weights_w7 = df7.set_index('ID')[w_col].reindex(ids_full)
# Normalize: mean=1 across non-NaN, NaN -> use 1.0
mean_w = np.nanmean(weights_w7)
weights_norm = weights_w7.fillna(mean_w).values / mean_w  # so unweighted recovers when all=1

print(f'Weight summary: mean={np.mean(weights_norm):.3f}, median={np.median(weights_norm):.3f}, '
      f'min={np.min(weights_norm):.3f}, max={np.max(weights_norm):.3f}, '
      f'fraction NaN replaced: {(weights_w7.isna().sum() / len(weights_w7)):.2%}')

valid_count = np.sum(~np.isnan(Y_full), axis=1)
keep = valid_count >= 3
Y = Y_full[keep]; ids = ids_full[keep]; w = weights_norm[keep]
N = Y.shape[0]
print(f'After filter: N = {N}')

times = np.arange(1, 8, dtype=float) - 4
X = np.column_stack([np.ones(7), times, times**2])

def fit_lcga_weighted(Y, K, X, weights, max_iter=300, tol=1e-5, n_starts=8, seed=42):
    N, T = Y.shape; P = X.shape[1]
    rng = np.random.default_rng(seed)
    # Normalize weights to sum to N (so weighted N == N)
    w = np.asarray(weights, dtype=float)
    w = w * (N / w.sum())
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
            # Weighted log-likelihood
            ll_w = (w * log_z).sum()
            if abs(ll_w - prev_ll) < tol: break
            prev_ll = ll_w
            # Weighted M-step
            wpost = post * w[:, None]
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

print('\n=== LCGA with longitudinal weight (WEIGHTB2w7) ===')
results_w = {}
for K in [2, 3, 4, 5]:
    m = fit_lcga_weighted(Y, K, X, w, max_iter=300, n_starts=8, seed=42)
    n_par = (K - 1) + K * X.shape[1] + 1
    bic = -2 * m['ll'] + n_par * np.log(N)
    abic = -2 * m['ll'] + n_par * np.log((N + 2) / 24)
    if K > 1:
        ent = -np.sum(m['post'] * np.log(m['post'] + 1e-300)) / (N * np.log(K))
        entropy = 1 - ent
    else:
        entropy = 1.0
    print(f'  K={K}: ll(weighted)={m["ll"]:.2f}, BIC={bic:.2f}, aBIC={abic:.2f}, '
          f'entropy={entropy:.3f}, props={m["pi"].round(3)}')
    results_w[K] = dict(model=m, bic=bic, abic=abic, entropy=entropy)

# Build comparison table (unweighted vs weighted) for K=4
fit_orig = pd.read_csv(os.path.join(OUT, 'fit_summary.csv'))
print('\nOriginal (unweighted) K=4:')
print(fit_orig[fit_orig.K == 4].to_string(index=False))

# Save weighted K=4
m4 = results_w[4]['model']
prof = pd.DataFrame({
    'class': range(1, 5),
    'proportion_weighted': m4['pi'],
    'beta_intercept': m4['beta'][:, 0],
    'beta_linear': m4['beta'][:, 1],
    'beta_quadratic': m4['beta'][:, 2],
})
mu4 = np.einsum('kp,tp->kt', m4['beta'], X)
for w_idx in range(7):
    prof[f'mean_w{w_idx+1}'] = mu4[:, w_idx]
prof.to_csv(os.path.join(SENS_OUT, 'weighted_class_profiles_K4.csv'), index=False, encoding='utf-8-sig')

# Compare: align classes between unweighted and weighted by mean trajectory ordering
print('\n=== Weighted K=4 class profiles ===')
print(prof.round(3).to_string(index=False))

# Save assignments
df_assign = pd.DataFrame({
    'ID': ids,
    'class_K4_weighted': m4['assign'] + 1,
    'maxprob_weighted': m4['post'].max(axis=1),
    'weight_used': w,
})
for k in range(4):
    df_assign[f'prob_class{k+1}_weighted'] = m4['post'][:, k]
df_assign.to_csv(os.path.join(SENS_OUT, 'weighted_assignments_K4.csv'), index=False, encoding='utf-8-sig')

# Compare to unweighted via crosstab
unw = pd.read_csv(os.path.join(OUT, 'class_assignments_K4.csv'))
unw_aligned = unw.set_index('ID')['class_K4'].reindex(ids)
crosstab = pd.crosstab(unw_aligned, m4['assign'] + 1, dropna=False, margins=True)
crosstab.to_csv(os.path.join(SENS_OUT, 'crosstab_unweighted_vs_weighted.csv'), encoding='utf-8-sig')
print('\n=== Crosstab: unweighted (rows) vs weighted (cols) class assignment K=4 ===')
print(crosstab.to_string())

# Agreement: Kappa
from sklearn.metrics import cohen_kappa_score
mask_both = unw_aligned.notna() & ~np.isnan(m4['assign'])
kappa = cohen_kappa_score(unw_aligned[mask_both].astype(int), (m4['assign'] + 1)[mask_both.values])
print(f'\nCohen kappa (unweighted vs weighted K=4 assignment): {kappa:.3f}')

with open(os.path.join(SENS_OUT, 'kappa.txt'), 'w', encoding='utf-8') as f:
    f.write(f'Cohen kappa unweighted vs weighted K=4: {kappa:.4f}\n')
    f.write(f'Weight column used: {w_col}\n')
    f.write(f'N: {N}\n')

print('\nDone. Sensitivity weighted analysis saved to', SENS_OUT)
