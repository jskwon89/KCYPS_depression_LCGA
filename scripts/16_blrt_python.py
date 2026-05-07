# -*- coding: utf-8 -*-
"""
Parametric Bootstrap Likelihood-Ratio Test (BLRT) for LCGA class enumeration.
Tests K-1 vs K under the null that K-1 is the true model.

Procedure (McLachlan 1987; Nylund-Gibson, Asparouhov & Muthen 2007):
  1. Fit K-1 model on observed data, record LL_(K-1).
  2. Fit K model on observed data, record LL_K. Compute LRT_obs = 2*(LL_K - LL_{K-1}).
  3. For b = 1..B:
     a. Simulate Y* from K-1 model parameters.
     b. Refit K-1 and K on Y*. Compute LRT_b = 2*(LL_K* - LL_{K-1}*).
  4. p_BLRT = (1 + sum(LRT_b >= LRT_obs)) / (B + 1)

We use B = 100 for time. Refits use multiple random starts.
"""
import pandas as pd
import numpy as np
import os, json, time
from scipy.special import logsumexp

ROOT = r'D:/2026/SCI/청소년패널조사'
DESC = os.path.join(ROOT, 'output', 'desc')
OUT  = os.path.join(ROOT, 'output', 'paperA_lcga')
BLRT_OUT = os.path.join(OUT, 'blrt')
os.makedirs(BLRT_OUT, exist_ok=True)

long_df = pd.read_parquet(os.path.join(DESC, 'long_panel.parquet'))
y = long_df.pivot_table(index='ID', columns='wave', values='depression')[list(range(1,8))]
ids = y.index.values
Y = y.values
times = np.arange(1, 8, dtype=float) - 4
X_design = np.column_stack([np.ones(7), times, times**2])

# Filter: at least 3 valid waves
valid_count = np.sum(~np.isnan(Y), axis=1)
keep = valid_count >= 3
Y = Y[keep]; ids = ids[keep]
N = Y.shape[0]
print(f'N = {N}')

# --- LCGA fitting ---
def fit_lcga(Y, K, X, max_iter=200, tol=1e-5, seed=0, n_starts=5, verbose=False):
    N, T = Y.shape; P = X.shape[1]
    rng = np.random.default_rng(seed)
    best = None
    for start in range(n_starts):
        if start == 0:
            mean_y = np.nanmean(Y, axis=1)
            cuts = np.nanquantile(mean_y, np.linspace(0, 1, K+1))
            init = np.zeros(N, dtype=int)
            for k in range(K):
                if k == K-1:
                    init[mean_y >= cuts[k]] = k
                else:
                    init[(mean_y >= cuts[k]) & (mean_y < cuts[k+1])] = k
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
            ll = log_z.sum()
            if abs(ll - prev_ll) < tol: break
            prev_ll = ll
            pi = post.sum(axis=0) / N
            pi = np.clip(pi, 1e-6, None); pi /= pi.sum()
            for k in range(K):
                XtWX = np.zeros((P, P)); XtWy = np.zeros(P)
                for i in range(N):
                    obs = ~np.isnan(Y[i])
                    if not obs.any(): continue
                    Xi = X[obs]; w = post[i, k]
                    XtWX += w * Xi.T @ Xi; XtWy += w * Xi.T @ Y[i][obs]
                if np.linalg.matrix_rank(XtWX) < P: continue
                beta[k] = np.linalg.solve(XtWX, XtWy)
            ss_total = w_total = 0.0
            mu = np.einsum('kp,tp->kt', beta, X)
            for i in range(N):
                obs = ~np.isnan(Y[i])
                if not obs.any(): continue
                resid2 = (Y[i][obs][None, :] - mu[:, obs])**2
                ss_total += (post[i] * resid2.sum(axis=1)).sum()
                w_total += (post[i] * obs.sum()).sum()
            sigma2 = ss_total / max(w_total, 1)
        if best is None or ll > best['ll']:
            best = dict(K=K, ll=ll, pi=pi.copy(), beta=beta.copy(), sigma2=sigma2,
                        post=post.copy(), assign=np.argmax(post, axis=1), iters=it+1)
    return best

# --- Simulate from K-1 model ---
def simulate_from_model(model, X, miss_pattern, rng):
    N, T = miss_pattern.shape
    K = model['K']
    pi = model['pi']
    beta = model['beta']
    sigma = np.sqrt(model['sigma2'])
    mu = np.einsum('kp,tp->kt', beta, X)
    cls = rng.choice(K, size=N, p=pi)
    Y_sim = np.full((N, T), np.nan)
    for i in range(N):
        eps = rng.normal(0, sigma, size=T)
        Y_sim[i] = mu[cls[i]] + eps
    Y_sim[~miss_pattern] = np.nan  # apply same missingness
    return Y_sim, cls

# --- Run BLRT for K-1 vs K ---
miss_pattern = ~np.isnan(Y)

# First fit observed K=1..5 once
print('\nFitting observed models K=1..5 ...')
obs_models = {}
for K in [1, 2, 3, 4, 5]:
    t0 = time.time()
    m = fit_lcga(Y, K, X_design, max_iter=300, n_starts=8, seed=42)
    obs_models[K] = m
    print(f'  K={K}: ll={m["ll"]:.2f}, sigma2={m["sigma2"]:.4f}, t={time.time()-t0:.1f}s')

# BLRT loop
B = 100  # bootstrap samples
blrt_results = []
for K_alt in [2, 3, 4, 5]:
    K_null = K_alt - 1
    obs_lrt = 2 * (obs_models[K_alt]['ll'] - obs_models[K_null]['ll'])
    print(f'\n=== BLRT K={K_null} vs K={K_alt}: observed LRT = {obs_lrt:.2f} ===')
    rng = np.random.default_rng(42 + K_alt)
    lrt_dist = []
    fail = 0
    t0 = time.time()
    for b in range(B):
        Y_sim, _ = simulate_from_model(obs_models[K_null], X_design, miss_pattern, rng)
        try:
            m_null = fit_lcga(Y_sim, K_null, X_design, max_iter=80, n_starts=2, seed=b)
            m_alt  = fit_lcga(Y_sim, K_alt,  X_design, max_iter=80, n_starts=2, seed=b)
            lrt_b = 2 * (m_alt['ll'] - m_null['ll'])
            lrt_dist.append(lrt_b)
        except Exception as e:
            fail += 1
            continue
        if (b + 1) % 20 == 0:
            elapsed = time.time() - t0
            print(f'  b={b+1}/{B}, mean LRT={np.mean(lrt_dist):.2f}, fail={fail}, elapsed={elapsed:.1f}s')
    lrt_dist = np.array(lrt_dist)
    p_blrt = (1 + np.sum(lrt_dist >= obs_lrt)) / (len(lrt_dist) + 1)
    blrt_results.append({
        'K_null': K_null, 'K_alt': K_alt,
        'obs_LRT': obs_lrt, 'B': len(lrt_dist), 'fail': fail,
        'mean_null_LRT': np.mean(lrt_dist), 'sd_null_LRT': np.std(lrt_dist),
        'q95_null_LRT': np.quantile(lrt_dist, 0.95),
        'p_BLRT': p_blrt,
    })
    print(f'  observed LRT={obs_lrt:.2f}, q95(null)={np.quantile(lrt_dist, 0.95):.2f}, p_BLRT={p_blrt:.4f}')
    np.save(os.path.join(BLRT_OUT, f'blrt_dist_K{K_null}vs{K_alt}.npy'), lrt_dist)

blrt_df = pd.DataFrame(blrt_results)
blrt_df.to_csv(os.path.join(BLRT_OUT, 'blrt_summary.csv'), index=False, encoding='utf-8-sig')
print('\n=== BLRT SUMMARY ===')
print(blrt_df.round(3).to_string(index=False))

# Decision rule: stop at smallest K where K vs K+1 BLRT is non-significant (p > .05).
# By Mplus/Nylund convention, the K-class solution is supported if K vs K-1 has p < .05
# AND K+1 vs K has p > .05 (or close).

print('\n=== K-class enumeration verdict (BLRT) ===')
for r in blrt_results:
    sig = '✓ K is supported over K-1' if r['p_BLRT'] < 0.05 else '✗ K not supported (parsimony favors K-1)'
    print(f'  K={r["K_alt"]}: p={r["p_BLRT"]:.4f} - {sig}')
print('\nDone.')
