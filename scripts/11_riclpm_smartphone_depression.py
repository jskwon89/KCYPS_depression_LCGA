# -*- coding: utf-8 -*-
"""
RI-CLPM: Smartphone addiction <-> Depression (continuous Likert, well-suited for RI-CLPM).
Compare: constrained vs free auto/cross-lag paths to evaluate stationarity.
"""
import pandas as pd
import numpy as np
import os, json
import semopy

ROOT = r'D:/2026/SCI/청소년패널조사'
DESC = os.path.join(ROOT, 'output', 'desc')
OUT = os.path.join(ROOT, 'output', 'paperB_riclpm')

long_df = pd.read_parquet(os.path.join(DESC, 'long_panel.parquet'))

X = long_df.pivot_table(index='ID', columns='wave', values='sm_addiction')
X.columns = [f'X{i}' for i in X.columns]
Y = long_df.pivot_table(index='ID', columns='wave', values='depression')
Y.columns = [f'Y{i}' for i in Y.columns]
dat = X.join(Y, how='inner').dropna(thresh=10)
print(f'Smartphone-Depression panel: N = {len(dat)}')
print(dat.describe().round(3).iloc[:, :8].to_string())

# Standardize for stable estimation
for c in dat.columns:
    dat[c] = (dat[c] - dat[c].mean()) / dat[c].std()
dat.to_csv(os.path.join(OUT, 'riclpm_sm_dep_data.csv'), encoding='utf-8-sig')

def make_riclpm_spec(constrain=True):
    """Build RI-CLPM lavaan spec for 7 waves. constrain=True applies equality constraints across waves."""
    body = ['# Random Intercepts',
            'RIx =~ 1*X1 + 1*X2 + 1*X3 + 1*X4 + 1*X5 + 1*X6 + 1*X7',
            'RIy =~ 1*Y1 + 1*Y2 + 1*Y3 + 1*Y4 + 1*Y5 + 1*Y6 + 1*Y7',
            '# Within-person factors']
    for t in range(1, 8):
        body.append(f'wx{t} =~ 1*X{t}')
        body.append(f'wy{t} =~ 1*Y{t}')
    body.append('# Constrain residual variance of observed to 0')
    for t in range(1, 8):
        body.append(f'X{t} ~~ 0*X{t}')
        body.append(f'Y{t} ~~ 0*Y{t}')
    body.append('# Cross-lagged paths')
    for t in range(2, 8):
        if constrain:
            body.append(f'wx{t} ~ a*wx{t-1} + b*wy{t-1}')
            body.append(f'wy{t} ~ c*wy{t-1} + d*wx{t-1}')
        else:
            body.append(f'wx{t} ~ a{t}*wx{t-1} + b{t}*wy{t-1}')
            body.append(f'wy{t} ~ c{t}*wy{t-1} + d{t}*wx{t-1}')
    body.append('# Within-time covariance')
    for t in range(2, 8):
        body.append(f'wx{t} ~~ wy{t}')
    body.append('# Inter-RI covariance')
    body.append('RIx ~~ RIy')
    return '\n'.join(body)

def fit_save(name, spec, dat):
    print(f'\n=== Fitting RI-CLPM: {name} ===')
    m = semopy.Model(spec)
    try:
        m.fit(dat, obj='FIML')
    except Exception as e:
        print(f'  fit failed: {e}; trying ML')
        m.fit(dat.dropna(), obj='ML')
    est = m.inspect(std_est=True)
    est.to_csv(os.path.join(OUT, f'riclpm_{name}_est.csv'), index=False, encoding='utf-8-sig')
    # Manual fit calc — degrees of freedom from spec
    try:
        info = semopy.calc_stats(m)
        info.to_csv(os.path.join(OUT, f'riclpm_{name}_fit.csv'), index=False, encoding='utf-8-sig')
        print(info.to_string())
    except Exception as e:
        print(f'  calc_stats failed: {e}')
    # Print key paths
    if 'op' in est.columns:
        keys = est[(est['op']=='~') & est['rval'].astype(str).str.startswith(('wx','wy'))]
        print(keys.to_string())
    return m, est

m_c, est_c = fit_save('sm_dep_constrained', make_riclpm_spec(constrain=True), dat)
m_f, est_f = fit_save('sm_dep_free', make_riclpm_spec(constrain=False), dat)

# Compare LL
ll_c = m_c.calc_lf()
ll_f = m_f.calc_lf()
n_par_c = len(m_c.param_vals)
n_par_f = len(m_f.param_vals)
chi2_diff = -2 * (ll_c - ll_f)
df_diff = n_par_f - n_par_c
from scipy import stats as sps
p_diff = 1 - sps.chi2.cdf(chi2_diff, df_diff) if df_diff > 0 else np.nan
print(f'\n=== LRT: free vs constrained ===')
print(f'  LL constrained: {-ll_c:.2f}, n_params: {n_par_c}')
print(f'  LL free: {-ll_f:.2f}, n_params: {n_par_f}')
print(f'  Chi² diff = {chi2_diff:.2f}, df = {df_diff}, p = {p_diff:.4g}')

# Save standardized estimates as a clean summary
summary_rows = []
for label, est in [('constrained', est_c), ('free', est_f)]:
    sub = est[(est['op']=='~') & (est['rval'].astype(str).str.startswith(('wx','wy')))]
    summary_rows.append(sub.assign(model=label))
summary = pd.concat(summary_rows)
summary.to_csv(os.path.join(OUT, 'riclpm_sm_dep_summary.csv'), index=False, encoding='utf-8-sig')

# Plot
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
fig, ax = plt.subplots(figsize=(8, 6))
# Read est_c keys
keys_c = est_c[(est_c['op']=='~') & est_c['rval'].astype(str).str.startswith(('wx','wy'))]
ax.axis('off')
ax.text(0.05, 0.95, 'RI-CLPM (constrained, standardized): Smartphone Addiction <-> Depression',
        fontsize=13, weight='bold', transform=ax.transAxes)
y = 0.85
keys_c2 = keys_c.drop_duplicates(subset=['lval','rval'])
for _, r in keys_c2.iterrows():
    ax.text(0.05, y, f'{r["lval"]} ~ {r["rval"]}: β = {r.get("Est. Std", r["Estimate"]):.3f}, p = {r["p-value"]}',
            fontsize=10, family='monospace', transform=ax.transAxes)
    y -= 0.05
plt.tight_layout()
plt.savefig(os.path.join(OUT, 'riclpm_sm_dep_summary_panel.png'), dpi=130, bbox_inches='tight')
plt.close()
print('Done.')
