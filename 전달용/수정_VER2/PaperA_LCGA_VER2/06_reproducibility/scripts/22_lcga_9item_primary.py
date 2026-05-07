# -*- coding: utf-8 -*-
"""
PRIMARY ANALYSIS (Reviewer P0 fix): LCGA on 9-item depression (excluding YPSY4E04 suicidal ideation)
to avoid criterion contamination with the W7 suicidal ideation distal outcome.

Pipeline:
  1. Load 9-item depression matrix
  2. Fit LCGA K=1..5 (Python EM, 8 random starts)
  3. Compute model fit (BIC, aBIC, AIC, entropy, smallest-class-prop)
  4. Save class assignments + posterior probs for K=4 (primary)
  5. Build comparison table to 10-item solution
"""
import pandas as pd
import numpy as np
import os, json, time
from scipy.special import logsumexp
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

ROOT = r'D:/2026/SCI/청소년패널조사'
DESC = os.path.join(ROOT, 'output', 'desc')
P9   = os.path.join(ROOT, 'output', 'paperA_lcga', 'primary_9item')
os.makedirs(P9, exist_ok=True)

long_df = pd.read_parquet(os.path.join(DESC, 'long_panel.parquet'))
y = long_df.pivot_table(index='ID', columns='wave', values='depression9')[list(range(1,8))]
ids = y.index.values
Y = y.values
print(f'9-item depression matrix: {Y.shape[0]} persons × {Y.shape[1]} waves')
print(f'Per-wave non-NaN: {np.sum(~np.isnan(Y), axis=0)}')

times = np.arange(1, 8, dtype=float) - 4
X = np.column_stack([np.ones(7), times, times**2])
valid_count = np.sum(~np.isnan(Y), axis=1)
keep = valid_count >= 3
Y = Y[keep]; ids = ids[keep]
N = Y.shape[0]
print(f'After filter (>=3 valid waves): N = {N}')

def fit_lcga(Y, K, X, max_iter=300, tol=1e-5, n_starts=8, seed=42):
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
            n_par = (K - 1) + K * P + 1
            bic = -2 * ll + n_par * np.log(N)
            abic = -2 * ll + n_par * np.log((N + 2) / 24)
            aic = -2 * ll + 2 * n_par
            if K > 1:
                ent = -np.sum(post * np.log(post + 1e-300)) / (N * np.log(K))
                entropy = 1 - ent
            else:
                entropy = 1.0
            best = dict(K=K, ll=ll, bic=bic, abic=abic, aic=aic, entropy=entropy,
                        pi=pi.copy(), beta=beta.copy(), sigma2=sigma2,
                        post=post.copy(), assign=np.argmax(post, axis=1))
    return best

results = {}
print('\nFitting K=1..5 on 9-item depression ...')
for K in [1, 2, 3, 4, 5]:
    t0 = time.time()
    m = fit_lcga(Y, K, X, n_starts=8, seed=42)
    results[K] = m
    print(f'  K={K}: ll={m["ll"]:.2f}, BIC={m["bic"]:.2f}, aBIC={m["abic"]:.2f}, entropy={m["entropy"]:.3f}, '
          f'sizes={(m["pi"]*100).round(1).tolist()}, t={time.time()-t0:.1f}s')

# Save fit summary
fit_rows = []
for K in [1, 2, 3, 4, 5]:
    m = results[K]
    fit_rows.append({
        'K': K, 'logL': m['ll'], 'BIC': m['bic'], 'aBIC': m['abic'], 'AIC': m['aic'],
        'entropy': m['entropy'], 'min_prop': m['pi'].min(), 'max_prop': m['pi'].max()
    })
fit_df = pd.DataFrame(fit_rows)
fit_df.to_csv(os.path.join(P9, 'fit_summary_9item.csv'), index=False, encoding='utf-8-sig')
print('\n=== Model fit summary (9-item) ===')
print(fit_df.round(2).to_string(index=False))

# Save K=4 (primary)
primary = results[4]
class_df = pd.DataFrame({'ID': ids, 'class_K4_9item': primary['assign'] + 1,
                          'maxprob': primary['post'].max(axis=1)})
for k in range(4):
    class_df[f'prob_class{k+1}'] = primary['post'][:, k]
class_df.to_csv(os.path.join(P9, 'class_assignments_K4_9item.csv'), index=False, encoding='utf-8-sig')

# Class profiles
profiles = pd.DataFrame({
    'class': range(1, 5),
    'proportion': primary['pi'],
    'beta_intercept': primary['beta'][:, 0],
    'beta_linear': primary['beta'][:, 1],
    'beta_quadratic': primary['beta'][:, 2],
})
mu = np.einsum('kp,tp->kt', primary['beta'], X)
for w in range(7):
    profiles[f'mean_w{w+1}'] = mu[:, w]
profiles.to_csv(os.path.join(P9, 'class_profiles_K4_9item.csv'), index=False, encoding='utf-8-sig')
print('\n=== K=4 (primary, 9-item) class profiles ===')
print(profiles.round(3).to_string(index=False))

# Compare to 10-item solution
old_class = pd.read_csv(os.path.join(ROOT, 'output', 'paperA_lcga', 'class_assignments_K4.csv'))
merged = class_df.merge(old_class[['ID', 'class_K4']].rename(columns={'class_K4': 'class_K4_10item'}),
                         on='ID', how='inner')

# Crosstab
ct = pd.crosstab(merged.class_K4_9item, merged.class_K4_10item, margins=True)
ct.to_csv(os.path.join(P9, 'crosstab_9item_vs_10item.csv'), encoding='utf-8-sig')
print('\n=== Crosstab: 9-item (rows) vs 10-item (cols) ===')
print(ct.to_string())

# Hungarian alignment for kappa
from scipy.optimize import linear_sum_assignment
ct_arr = ct.iloc[:-1, :-1].values
# Maximize match -> minimize negative
row_ind, col_ind = linear_sum_assignment(-ct_arr)
print(f'\nLabel alignment (9-item -> 10-item): {dict(zip(row_ind+1, col_ind+1))}')

# Apply alignment then compute kappa
mapping = dict(zip(row_ind + 1, col_ind + 1))  # 9-item label -> 10-item label
merged['class_9item_aligned'] = merged['class_K4_9item'].map(mapping)
from sklearn.metrics import cohen_kappa_score
kappa_aligned = cohen_kappa_score(merged['class_K4_10item'], merged['class_9item_aligned'])
agreement = (merged['class_K4_10item'] == merged['class_9item_aligned']).mean()
print(f'After Hungarian alignment: agreement = {agreement:.4f}, Cohen kappa = {kappa_aligned:.4f}')

with open(os.path.join(P9, 'comparison_to_10item.txt'), 'w', encoding='utf-8') as f:
    f.write(f'9-item vs 10-item comparison\n')
    f.write(f'Total individuals: {len(merged)}\n')
    f.write(f'Hungarian label alignment (9-item label -> 10-item label): {mapping}\n')
    f.write(f'Post-alignment agreement: {agreement:.4f}\n')
    f.write(f'Post-alignment Cohen kappa: {kappa_aligned:.4f}\n')

# Plot trajectories
fig, ax = plt.subplots(figsize=(10, 6))
labels_4 = {1: 'C1', 2: 'C2', 3: 'C3', 4: 'C4'}
colors = ['#ff7f0e', '#2ca02c', '#1f77b4', '#d62728']
times_w = np.arange(1, 8)
for k in range(4):
    ax.plot(times_w, mu[k], '-o', color=colors[k], linewidth=3, markersize=10,
            label=f'C{k+1} ({primary["pi"][k]*100:.1f}%)')
ax.set_xlabel('Wave')
ax.set_ylabel('9-item depression score (mean of items E01-E03, E05-E10; 1-4)')
ax.set_title('LCGA 9-item depression (criterion-clean)\nPRIMARY ANALYSIS, N={}'.format(N))
ax.legend(loc='upper right')
ax.grid(alpha=0.3)
ax.set_ylim(1.0, 3.0)
plt.tight_layout()
plt.savefig(os.path.join(P9, 'lcga_K4_9item_trajectories.png'), dpi=300, bbox_inches='tight')
plt.close()

print('\nDone. Primary 9-item LCGA results in', P9)
