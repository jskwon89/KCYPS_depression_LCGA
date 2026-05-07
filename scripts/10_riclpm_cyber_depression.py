# -*- coding: utf-8 -*-
"""
Random-Intercept Cross-Lagged Panel Model (RI-CLPM) - cyber delinquency <-> depression.
Hamaker, Kuiper & Grasman (2015).

Spec (7-wave RI-CLPM):
  X_it = RI_X_i + WX_it     (within-person fluctuation)
  Y_it = RI_Y_i + WY_it
  WX_it = a*WX_(i,t-1) + b*WY_(i,t-1) + e_X_it     -> b is "Y -> X" within-person cross-lag
  WY_it = c*WY_(i,t-1) + d*WX_(i,t-1) + e_Y_it     -> d is "X -> Y" within-person cross-lag

We use semopy if available; otherwise fall back to lavaan-like manual implementation
via OpenMx-equivalent (which we don't have)... we'll implement via semopy.
"""
import pandas as pd
import numpy as np
import os, json
import warnings
warnings.filterwarnings('ignore')

ROOT = r'D:/2026/SCI/청소년패널조사'
DESC = os.path.join(ROOT, 'output', 'desc')
OUT = os.path.join(ROOT, 'output', 'paperB_riclpm')
os.makedirs(OUT, exist_ok=True)

long_df = pd.read_parquet(os.path.join(DESC, 'long_panel.parquet'))

# Pivot to wide for X (cyber_delin) and Y (depression)
def to_wide(var):
    w = long_df.pivot_table(index='ID', columns='wave', values=var)
    w.columns = [f'{var}_w{c}' for c in w.columns]
    return w

X = to_wide('cyber_delin')
Y = to_wide('depression')

# Fix: cyber_delin in long_df is COUNT_ANY which counted NaN as not-1, giving zeros even when missing.
# Re-read raw items and recompute properly
import pyreadstat
ROOT2 = ROOT
for var_root, prefix, n_items in [('cyber', 'YDLQ2A', 15), ('off', 'YDLQ1A', 15)]:
    cols = []
    for w in range(1, 8):
        df, _ = pyreadstat.read_dta(os.path.join(ROOT2, f'KCYPS2018m1Yw{w}.dta'))
        items = [f'{prefix}{i:02d}w{w}' for i in range(1, n_items+1)]
        sub = df.reindex(columns=items).copy()
        # KCYPS coding for delinquency: 1=없다, 2~5=빈도. NaN handling: <0 missing, ≥90 missing
        sub = sub.mask((sub < 1) | (sub > 6), other=np.nan)
        # If majority of items missing -> person did not respond
        n_missing = sub.isna().sum(axis=1)
        any_resp = n_missing < (n_items - 2)
        # count_any: number of items with any occurrence (>=2)
        count = (sub >= 2).sum(axis=1).astype(float)
        count[~any_resp] = np.nan
        cols.append(pd.Series(count.values, index=df['ID'].values, name=f'{var_root}_w{w}'))
    df_var = pd.concat(cols, axis=1)
    df_var.to_parquet(os.path.join(OUT, f'{var_root}_count_wide.parquet'))
    print(f'{var_root} count: shape {df_var.shape}, n nonzero = {(df_var>0).sum().to_dict()}, mean = {df_var.mean().to_dict()}')

# Reload corrected variables
X = pd.read_parquet(os.path.join(OUT, 'cyber_count_wide.parquet'))
X.columns = [f'X{i+1}' for i in range(7)]
# Log1p transform for skewed count data
X = np.log1p(X)
Y_corr = Y.copy()
Y_corr.columns = [f'Y{i+1}' for i in range(7)]

dat = X.join(Y_corr, how='inner').dropna(thresh=8)  # require >=8 of 14 obs
print(f'\nFinal panel: N = {len(dat)}')
print(dat.describe().round(3).to_string())

# Save
dat.to_csv(os.path.join(OUT, 'riclpm_data.csv'), encoding='utf-8-sig')

# Try semopy
try:
    import semopy
    print('\nsemopy version:', semopy.__version__)
    HAVE_SEMOPY = True
except ImportError:
    HAVE_SEMOPY = False
    print('\nsemopy not available - will use manual RI-CLPM via state-space')

if HAVE_SEMOPY:
    # RI-CLPM specification in lavaan-like syntax
    T = 7
    spec = """
    # Random intercepts as latent factors
    RIx =~ 1*X1 + 1*X2 + 1*X3 + 1*X4 + 1*X5 + 1*X6 + 1*X7
    RIy =~ 1*Y1 + 1*Y2 + 1*Y3 + 1*Y4 + 1*Y5 + 1*Y6 + 1*Y7

    # Within-person centered factors
    wx1 =~ 1*X1
    wx2 =~ 1*X2
    wx3 =~ 1*X3
    wx4 =~ 1*X4
    wx5 =~ 1*X5
    wx6 =~ 1*X6
    wx7 =~ 1*X7
    wy1 =~ 1*Y1
    wy2 =~ 1*Y2
    wy3 =~ 1*Y3
    wy4 =~ 1*Y4
    wy5 =~ 1*Y5
    wy6 =~ 1*Y6
    wy7 =~ 1*Y7

    # Constrain observed variances to 0
    X1 ~~ 0*X1
    X2 ~~ 0*X2
    X3 ~~ 0*X3
    X4 ~~ 0*X4
    X5 ~~ 0*X5
    X6 ~~ 0*X6
    X7 ~~ 0*X7
    Y1 ~~ 0*Y1
    Y2 ~~ 0*Y2
    Y3 ~~ 0*Y3
    Y4 ~~ 0*Y4
    Y5 ~~ 0*Y5
    Y6 ~~ 0*Y6
    Y7 ~~ 0*Y7

    # Auto-regressive paths (constrained equal across waves for parsimony)
    wx2 ~ a*wx1 + b*wy1
    wx3 ~ a*wx2 + b*wy2
    wx4 ~ a*wx3 + b*wy3
    wx5 ~ a*wx4 + b*wy4
    wx6 ~ a*wx5 + b*wy5
    wx7 ~ a*wx6 + b*wy6
    wy2 ~ c*wy1 + d*wx1
    wy3 ~ c*wy2 + d*wx2
    wy4 ~ c*wy3 + d*wx3
    wy5 ~ c*wy4 + d*wx4
    wy6 ~ c*wy5 + d*wx5
    wy7 ~ c*wy6 + d*wx6

    # RI uncorrelated with first within-person factor (defaults), but free covariance
    RIx ~~ RIy

    # Within-time covariance of disturbances
    wx2 ~~ wy2
    wx3 ~~ wy3
    wx4 ~~ wy4
    wx5 ~~ wy5
    wx6 ~~ wy6
    wx7 ~~ wy7
    """
    try:
        m = semopy.Model(spec)
        m.fit(dat, obj='FIML')
        est = m.inspect()
        est.to_csv(os.path.join(OUT, 'riclpm_estimates.csv'), index=False, encoding='utf-8-sig')
        stats = semopy.calc_stats(m)
        stats.to_csv(os.path.join(OUT, 'riclpm_fit.csv'), index=False, encoding='utf-8-sig')
        print('\nFit indices:')
        print(stats.to_string())
        print('\nKey paths (a/b/c/d):')
        print(est[est['op'] == '~'].to_string())
    except Exception as e:
        print(f'semopy failed: {e}')
        import traceback; traceback.print_exc()

print('\nRI-CLPM script done')
