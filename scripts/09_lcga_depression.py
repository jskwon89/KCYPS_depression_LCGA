# -*- coding: utf-8 -*-
"""
Latent Class Growth Analysis (LCGA) of depression trajectories (KCYPS m1 cohort, 7 waves).

LCGA / Group-Based Trajectory Model (Nagin):
  Y_it = β0_g + β1_g * t + β2_g * t^2 + ε_it
  where g = latent class, ε_it ~ N(0, σ²) i.i.d.
  Within-class variance of growth parameters = 0.

EM algorithm:
  E-step: P(class = g | Y_i) ∝ π_g * Π_t f(Y_it; μ_gt, σ²)
  M-step:
    π_g = Σ_i posterior_ig / N
    β_g = WLS of {(Y_it, t, t²)} weighted by posterior_ig
    σ² = Σ_i Σ_g posterior_ig * Σ_t (Y_it - μ_gt)²  /  Σ_i posterior_ig * n_obs_i

We use FIML (use available cases per person — no listwise deletion).
"""
import pandas as pd
import numpy as np
import os, json, sys
from scipy.special import logsumexp
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

ROOT = r'D:/2026/SCI/청소년패널조사'
OUT = os.path.join(ROOT, 'output', 'desc')
LCGA_OUT = os.path.join(ROOT, 'output', 'paperA_lcga')
os.makedirs(LCGA_OUT, exist_ok=True)

long_df = pd.read_parquet(os.path.join(OUT, 'long_panel.parquet'))

# Prepare wide depression matrix
y = long_df.pivot_table(index='ID', columns='wave', values='depression')[list(range(1,8))]
ids = y.index.values
Y = y.values  # (N, T) with NaN for missing
N, T = Y.shape
print(f'Depression matrix: {N} persons × {T} waves')
print(f'Per-wave N (non-NaN): {np.sum(~np.isnan(Y), axis=0)}')
print(f'Persons with >=1 wave: {(np.sum(~np.isnan(Y), axis=1) >= 1).sum()}')
print(f'Persons with all 7: {(np.sum(~np.isnan(Y), axis=1) == 7).sum()}')

# Time codes: center on wave 4 (mid-point) for stability
times = np.arange(1, T+1, dtype=float)
t_centered = times - 4  # so wave 1 = -3, wave 7 = +3
# Polynomial design matrix (intercept, linear, quadratic)
# Standardize quadratic
X_design_full = np.column_stack([np.ones(T), t_centered, t_centered**2])  # T x 3

# Filter: persons with at least 3 waves of depression
valid_count = np.sum(~np.isnan(Y), axis=1)
keep = valid_count >= 3
Y = Y[keep]
ids = ids[keep]
N = Y.shape[0]
print(f'After filter (>=3 waves): N = {N}')

def fit_lcga(Y, K, X, max_iter=300, tol=1e-5, seed=0, n_starts=10):
    """Fit LCGA with K classes by EM. Returns best result over n_starts random initializations."""
    N, T = Y.shape
    P = X.shape[1]
    rng = np.random.default_rng(seed)
    best = None
    for start in range(n_starts):
        # Initialize: random class assignment, fit OLS per class
        if start == 0:
            # init: cluster on observed mean depression
            mean_dep = np.nanmean(Y, axis=1)
            quantiles = np.linspace(0, 1, K+1)
            cuts = np.nanquantile(mean_dep, quantiles)
            init_class = np.zeros(N, dtype=int)
            for k in range(K):
                if k == K - 1:
                    init_class[mean_dep >= cuts[k]] = k
                else:
                    init_class[(mean_dep >= cuts[k]) & (mean_dep < cuts[k+1])] = k
        else:
            init_class = rng.integers(0, K, size=N)
        # Initial class params via OLS per class
        beta = np.zeros((K, P))
        for k in range(K):
            mask = init_class == k
            if mask.sum() < 5:
                beta[k] = np.array([np.nanmean(Y), 0, 0])
                continue
            # Stack obs in this class
            x_list, y_list = [], []
            for i in np.where(mask)[0]:
                yi = Y[i]
                obs = ~np.isnan(yi)
                x_list.append(X[obs])
                y_list.append(yi[obs])
            Xstack = np.vstack(x_list); ystack = np.concatenate(y_list)
            beta[k] = np.linalg.lstsq(Xstack, ystack, rcond=None)[0]
        sigma2 = np.nanvar(Y)  # initial
        pi = np.full(K, 1.0/K)
        prev_ll = -np.inf
        log_post = None
        for it in range(max_iter):
            # E-step
            mu = np.einsum('kp,tp->kt', beta, X)  # K x T
            # log-lik per (i, k): sum_{t obs} log N(Y_it | mu_kt, sigma2)
            log_lik_ik = np.full((N, K), -np.inf)
            for i in range(N):
                yi = Y[i]
                obs = ~np.isnan(yi)
                if obs.sum() == 0: continue
                resid2 = (yi[obs][None, :] - mu[:, obs])**2  # K x T_obs
                log_lik_ik[i] = -0.5*obs.sum()*np.log(2*np.pi*sigma2) - 0.5*resid2.sum(axis=1)/sigma2
            log_post_unnorm = log_lik_ik + np.log(pi + 1e-300)
            log_z = logsumexp(log_post_unnorm, axis=1)
            log_post = log_post_unnorm - log_z[:, None]
            post = np.exp(log_post)  # N x K
            ll = log_z.sum()
            if abs(ll - prev_ll) < tol:
                break
            prev_ll = ll
            # M-step
            pi = post.sum(axis=0) / N
            pi = np.clip(pi, 1e-6, None); pi /= pi.sum()
            # Update beta_k via weighted regression
            for k in range(K):
                # build weighted OLS
                XtWX = np.zeros((P, P)); XtWy = np.zeros(P)
                for i in range(N):
                    yi = Y[i]
                    obs = ~np.isnan(yi)
                    if obs.sum() == 0: continue
                    Xi = X[obs]; yi_obs = yi[obs]
                    w = post[i, k]
                    XtWX += w * Xi.T @ Xi
                    XtWy += w * Xi.T @ yi_obs
                if np.linalg.matrix_rank(XtWX) < P:
                    continue
                beta[k] = np.linalg.solve(XtWX, XtWy)
            # Update sigma²
            ss_total = 0.0; w_total = 0.0
            mu = np.einsum('kp,tp->kt', beta, X)
            for i in range(N):
                yi = Y[i]
                obs = ~np.isnan(yi)
                if obs.sum() == 0: continue
                resid2 = (yi[obs][None, :] - mu[:, obs])**2
                ss_total += (post[i] * resid2.sum(axis=1)).sum()
                w_total += (post[i] * obs.sum()).sum()
            sigma2 = ss_total / max(w_total, 1)
        # done EM
        # BIC: -2*LL + k_params * log(N)
        # parameters: K-1 (mixing) + K*P (betas) + 1 (sigma²)
        k_params = (K - 1) + K * P + 1
        bic = -2 * ll + k_params * np.log(N)
        aic = -2 * ll + 2 * k_params
        # adjusted BIC: log((N+2)/24)
        abic = -2 * ll + k_params * np.log((N + 2) / 24)
        # Entropy
        if K > 1:
            ent = -np.sum(post * np.log(post + 1e-300)) / (N * np.log(K))
            entropy = 1 - ent
        else:
            entropy = 1.0
        result = dict(K=K, ll=ll, bic=bic, aic=aic, abic=abic, entropy=entropy,
                      pi=pi.copy(), beta=beta.copy(), sigma2=sigma2, post=post.copy(),
                      class_assignment=np.argmax(post, axis=1), iters=it+1, start=start)
        if best is None or ll > best['ll']:
            best = result
    return best

# Fit K=1..5
results = {}
for K in [1, 2, 3, 4, 5]:
    print(f'\n--- Fitting LCGA K={K} ---')
    res = fit_lcga(Y, K, X_design_full, max_iter=300, tol=1e-5, seed=42, n_starts=8)
    print(f'  best LL = {res["ll"]:.2f}, BIC = {res["bic"]:.2f}, aBIC = {res["abic"]:.2f}, entropy = {res["entropy"]:.3f}')
    print(f'  class proportions: {res["pi"].round(3)}')
    print(f'  mean trajectories (wave 1 -> 7):')
    mu = np.einsum('kp,tp->kt', res['beta'], X_design_full)
    for k in range(K):
        print(f'    class {k+1} (n={res["pi"][k]*100:.1f}%): {mu[k].round(3)}')
    results[K] = res

# Save model fit summary
fit_summary = pd.DataFrame([
    {'K': k, 'logL': r['ll'], 'BIC': r['bic'], 'aBIC': r['abic'], 'AIC': r['aic'],
     'entropy': r['entropy'],
     'min_prop': r['pi'].min(), 'max_prop': r['pi'].max()}
    for k, r in results.items()
])
fit_summary.to_csv(os.path.join(LCGA_OUT, 'fit_summary.csv'), index=False, encoding='utf-8-sig')
print('\n=== Model fit summary ===')
print(fit_summary.round(2).to_string(index=False))

# Choose best K based on BIC + interpretability
best_K = int(fit_summary.sort_values('BIC').iloc[0]['K'])
print(f'\nBest K by BIC: {best_K}')
# Often elbow is at smaller K — let's also save K=4 for convention
chosen = results[best_K]

# Plot mean trajectories for K=2..5
fig, axes = plt.subplots(1, 4, figsize=(18, 4), sharey=True)
for idx, K in enumerate([2, 3, 4, 5]):
    ax = axes[idx]
    res = results[K]
    mu = np.einsum('kp,tp->kt', res['beta'], X_design_full)
    for k in range(K):
        ax.plot(times, mu[k], '-o', label=f'C{k+1} ({res["pi"][k]*100:.1f}%)', linewidth=2)
    ax.set_title(f'K={K}, BIC={res["bic"]:.0f}, ent={res["entropy"]:.2f}')
    ax.set_xlabel('Wave (1=중1 2018, 7=2024)')
    ax.set_ylabel('우울 (1~4)')
    ax.legend(fontsize=8)
    ax.grid(alpha=0.3)
plt.suptitle('LCGA: Mean Depression Trajectories by Latent Class (K=2..5)')
plt.tight_layout()
plt.savefig(os.path.join(LCGA_OUT, 'lcga_trajectories_K2_5.png'), dpi=150, bbox_inches='tight')
plt.close()

# Also save individual class assignments for chosen K (for distal analyses)
# Use K=4 as primary model unless BIC strongly favors otherwise
PRIMARY_K = 4 if best_K >= 3 else 3
print(f'\nUsing primary K = {PRIMARY_K} for distal outcome analyses')
chosen = results[PRIMARY_K]
class_df = pd.DataFrame({
    'ID': ids,
    f'class_K{PRIMARY_K}': chosen['class_assignment'] + 1,
    f'maxprob_K{PRIMARY_K}': chosen['post'].max(axis=1),
})
for k in range(PRIMARY_K):
    class_df[f'prob_class{k+1}'] = chosen['post'][:, k]
class_df.to_csv(os.path.join(LCGA_OUT, f'class_assignments_K{PRIMARY_K}.csv'), index=False, encoding='utf-8-sig')

# Class profiles
profiles = pd.DataFrame({
    'class': range(1, PRIMARY_K+1),
    'proportion': chosen['pi'],
    'beta_intercept': chosen['beta'][:, 0],
    'beta_linear': chosen['beta'][:, 1],
    'beta_quadratic': chosen['beta'][:, 2],
})
mu_traj = np.einsum('kp,tp->kt', chosen['beta'], X_design_full)
for w in range(T):
    profiles[f'mean_w{w+1}'] = mu_traj[:, w]
profiles.to_csv(os.path.join(LCGA_OUT, f'class_profiles_K{PRIMARY_K}.csv'), index=False, encoding='utf-8-sig')
print(f'\nClass profiles (K={PRIMARY_K}):')
print(profiles.round(3).to_string(index=False))

# Plot the chosen K with raw-data overlay
fig, ax = plt.subplots(figsize=(10, 6))
mu_traj = np.einsum('kp,tp->kt', chosen['beta'], X_design_full)
colors = ['#2ca02c', '#1f77b4', '#ff7f0e', '#d62728', '#9467bd']
for k in range(PRIMARY_K):
    # Spaghetti subset
    members = np.where(chosen['class_assignment'] == k)[0]
    if len(members) > 80:
        members = np.random.choice(members, size=80, replace=False)
    for i in members:
        ax.plot(times, Y[i], color=colors[k], alpha=0.05, linewidth=0.6)
    ax.plot(times, mu_traj[k], color=colors[k], linewidth=3,
            label=f'C{k+1} ({chosen["pi"][k]*100:.1f}%)', marker='o', markersize=8)
ax.set_xlabel('Wave (1=2018 중1 ~ 7=2024)', fontsize=12)
ax.set_ylabel('우울 (Mean of 10 items, 1~4 scale)', fontsize=12)
ax.set_title(f'Latent Class Growth Analysis: Depression Trajectories (K={PRIMARY_K})\nN={N}, FIML', fontsize=13)
ax.legend(loc='upper right')
ax.grid(alpha=0.3)
plt.tight_layout()
plt.savefig(os.path.join(LCGA_OUT, f'lcga_trajectories_K{PRIMARY_K}_with_spaghetti.png'), dpi=150, bbox_inches='tight')
plt.close()

print(f'\nLCGA results saved to {LCGA_OUT}')
