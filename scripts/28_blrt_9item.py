# -*- coding: utf-8 -*-
"""
Bootstrap LRT (BLRT) for 9-item primary depression LCGA.
K=1 vs K=2 only. B=30 to fit within ~10-15 minutes.

Observed:
  K=1 logL = -14143.19, K=2 logL = -13055.84
  LRT_obs = 2 × (13055.84 - 14143.19) = 2174.70
"""
import pandas as pd
import numpy as np
import os, time
from scipy.special import logsumexp

ROOT = r'D:/2026/SCI/청소년패널조사'
DESC = os.path.join(ROOT, 'output', 'desc')
P9   = os.path.join(ROOT, 'output', 'paperA_lcga', 'primary_9item')
BLRT_OUT = os.path.join(P9, 'blrt')
os.makedirs(BLRT_OUT, exist_ok=True)

long_df = pd.read_parquet(os.path.join(DESC, 'long_panel.parquet'))
y = long_df.pivot_table(index='ID', columns='wave', values='depression9')[list(range(1,8))]
ids = y.index.values; Y = y.values
times = np.arange(1, 8, dtype=float) - 4
X = np.column_stack([np.ones(7), times, times**2])
keep = np.sum(~np.isnan(Y), axis=1) >= 3
Y = Y[keep]; ids = ids[keep]
N = Y.shape[0]
print(f'N = {N}')

def fit_lcga(Y, K, X, max_iter=200, tol=1e-5, n_starts=3, seed=42):
    N, T = Y.shape; P = X.shape[1]
    rng = np.random.default_rng(seed)
    best = None
    for start in range(n_starts):
        if start == 0:
            mean_y = np.nanmean(Y, axis=1)
            cuts = np.nanquantile(mean_y, np.linspace(0, 1, K+1))
            init = np.zeros(N, dtype=int)
            for k in range(K):
                init[(mean_y >= cuts[k]) & (mean_y <= cuts[k+1] + 1e-9)] = k
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
                w_total  += (post[i] * obs.sum()).sum()
            sigma2 = ss_total / max(w_total, 1)
        if best is None or ll > best['ll']:
            best = dict(K=K, ll=ll, pi=pi.copy(), beta=beta.copy(),
                        sigma2=sigma2, post=post.copy())
    return best

# Fit observed K=1, K=2
print('Fitting observed K=1, K=2 ...')
m1 = fit_lcga(Y, 1, X, n_starts=3)
m2 = fit_lcga(Y, 2, X, n_starts=5, max_iter=300)
obs_lrt = 2 * (m2['ll'] - m1['ll'])
print(f'  K=1 ll = {m1["ll"]:.2f}')
print(f'  K=2 ll = {m2["ll"]:.2f}')
print(f'  LRT_obs = {obs_lrt:.2f}')

# Bootstrap
miss = ~np.isnan(Y)
B = 30
print(f'\nBootstrap LRT under K=1 null (B={B}) ...')
rng = np.random.default_rng(42)
beta1 = m1['beta']; sigma1 = np.sqrt(m1['sigma2'])
mu1 = np.einsum('kp,tp->kt', beta1, X).flatten()  # K=1 has single trajectory

lrt_dist = []
fail = 0
t0 = time.time()
for b in range(B):
    Y_sim = np.full_like(Y, np.nan)
    eps = rng.normal(0, sigma1, size=Y.shape)
    Y_sim_full = mu1[None, :] + eps
    Y_sim[miss] = Y_sim_full[miss]
    try:
        mb1 = fit_lcga(Y_sim, 1, X, n_starts=2, max_iter=80)
        mb2 = fit_lcga(Y_sim, 2, X, n_starts=2, max_iter=80)
        lrt_b = 2 * (mb2['ll'] - mb1['ll'])
        lrt_dist.append(lrt_b)
    except Exception as e:
        fail += 1
        continue
    if (b + 1) % 5 == 0:
        elapsed = time.time() - t0
        print(f'  b={b+1}/{B}, mean LRT={np.mean(lrt_dist):.2f}, max={np.max(lrt_dist):.2f}, '
              f'fail={fail}, elapsed={elapsed:.0f}s')

lrt_dist = np.array(lrt_dist)
p_blrt = (1 + np.sum(lrt_dist >= obs_lrt)) / (len(lrt_dist) + 1)
print(f'\n=== BLRT 9-item K=1 vs K=2 ===')
print(f'  observed LRT = {obs_lrt:.2f}')
print(f'  null bootstrap (B={len(lrt_dist)}): mean={np.mean(lrt_dist):.2f}, '
      f'sd={np.std(lrt_dist):.2f}, max={np.max(lrt_dist):.2f}')
print(f'  p_BLRT = {p_blrt:.4f}')

# Save
summary = pd.DataFrame([{
    'K_null': 1, 'K_alt': 2,
    'obs_LRT': obs_lrt, 'B': len(lrt_dist), 'fail': fail,
    'mean_null_LRT': np.mean(lrt_dist), 'sd_null_LRT': np.std(lrt_dist),
    'max_null_LRT': np.max(lrt_dist),
    'q95_null_LRT': np.quantile(lrt_dist, 0.95) if len(lrt_dist) >= 20 else np.nan,
    'p_BLRT': p_blrt,
    'note': '9-item primary; B=30 parametric bootstrap; K>2 not bootstrapped (BIC-based decision)',
}])
summary.to_csv(os.path.join(BLRT_OUT, 'blrt_summary_9item.csv'), index=False, encoding='utf-8-sig')
np.save(os.path.join(BLRT_OUT, 'blrt_dist_K1vs2_9item.npy'), lrt_dist)
print(f'\nSaved summary to {BLRT_OUT}')
