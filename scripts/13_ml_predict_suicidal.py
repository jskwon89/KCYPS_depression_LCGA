# -*- coding: utf-8 -*-
"""
ML prediction of late-adolescent / emerging-adult outcomes from W1 (Grade 7) baseline.

Outcomes (W7, age 19):
  (a) Suicidal ideation (YPSY4E04w7 >= 2)  → binary
  (b) High depression (YPSY4E score >= 75th pctile)
  (c) Belonging to "high-risk" trajectory class (C1 or C4 from LCGA K=4)

Models: Logistic-L1, Random Forest, XGBoost, LightGBM
Eval: nested 5-fold CV, AUC, calibration, SHAP global importance.
"""
import pandas as pd
import numpy as np
import os, json
from sklearn.model_selection import StratifiedKFold, cross_val_score
from sklearn.linear_model import LogisticRegressionCV
from sklearn.ensemble import RandomForestClassifier, GradientBoostingClassifier
from sklearn.metrics import roc_auc_score, average_precision_score, brier_score_loss
from sklearn.calibration import calibration_curve
from sklearn.preprocessing import StandardScaler
from sklearn.impute import SimpleImputer
from sklearn.pipeline import Pipeline
import warnings
warnings.filterwarnings('ignore')
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

ROOT = r'D:/2026/SCI/청소년패널조사'
DESC = os.path.join(ROOT, 'output', 'desc')
LCGA_OUT = os.path.join(ROOT, 'output', 'paperA_lcga')
ML_OUT = os.path.join(ROOT, 'output', 'paperC_ml')
os.makedirs(ML_OUT, exist_ok=True)

# Build feature matrix from W1 data
import pyreadstat
df1, meta1 = pyreadstat.read_dta(os.path.join(ROOT, 'KCYPS2018m1Yw1.dta'))
df7, _ = pyreadstat.read_dta(os.path.join(ROOT, 'KCYPS2018m1Yw7.dta'))

long_df = pd.read_parquet(os.path.join(DESC, 'long_panel.parquet'))
class_df = pd.read_csv(os.path.join(LCGA_OUT, 'class_assignments_K4.csv'), encoding='utf-8-sig')

# Get all w1 scale scores from long_df
w1 = long_df[long_df.wave == 1].set_index('ID')

# Add raw item-level features for richer ML (selected validated scales W1)
# We'll use scale-level features primarily for interpretability
feat_cols = ['YGENDER', 'parent_warmth', 'parent_reject', 'parent_coerce', 'parent_chaos',
             'parent_structure', 'parent_autonomy', 'peer_attach', 'teacher_rel',
             'self_esteem', 'sm_addiction', 'depression', 'life_sat', 'happy', 'grit',
             'aggression', 'somatic', 'withdrawal', 'attention', 'YBRT1A', 'YBRT1B']
X = w1[feat_cols].copy()
X['female'] = (X['YGENDER'] == 2).astype(float)
X = X.drop(columns=['YGENDER'])

# Add area, weight
for col in ['ARA1A', 'ARA2A']:
    src_col = f'{col}w1'
    if src_col in df1.columns:
        X[col] = df1.set_index('ID')[src_col].reindex(X.index).values

# Outcomes at W7
df7['suic_w7'] = df7['YPSY4E04w7'].where((df7['YPSY4E04w7']>=1) & (df7['YPSY4E04w7']<=4))
y_suic = (df7.set_index('ID')['suic_w7'] >= 2).astype(float)

# High depression W7
dep_w7 = long_df[long_df.wave==7].set_index('ID')['depression']
high_dep = (dep_w7 >= dep_w7.quantile(0.75)).astype(float)

# High-risk class membership (C1 or C4 from LCGA)
class_dict = class_df.set_index('ID')['class_K4']
hr_class = class_dict.isin([1, 4]).astype(float)

outcomes = {
    'suicidal_ideation_w7': y_suic,
    'high_depression_w7': high_dep,
    'high_risk_trajectory_C1C4': hr_class,
}

# Align indices
common_idx = X.index
results_all = []

for outcome_name, y_series in outcomes.items():
    print(f'\n=== Outcome: {outcome_name} ===')
    y = y_series.reindex(common_idx)
    mask = y.notna() & X.notna().all(axis=1)
    Xm = X[mask].copy()
    ym = y[mask].astype(int).values
    print(f'  N = {len(Xm)}, positive rate = {ym.mean():.3f}')
    if len(Xm) < 100 or ym.sum() < 30:
        print('  too few cases, skip')
        continue
    # CV
    skf = StratifiedKFold(n_splits=5, shuffle=True, random_state=42)
    models = {
        'LogisticL1': Pipeline([
            ('imp', SimpleImputer(strategy='median')),
            ('sc', StandardScaler()),
            ('clf', LogisticRegressionCV(Cs=20, penalty='l1', solver='saga', max_iter=2000, scoring='roc_auc')),
        ]),
        'RandomForest': Pipeline([
            ('imp', SimpleImputer(strategy='median')),
            ('clf', RandomForestClassifier(n_estimators=500, max_depth=10, min_samples_leaf=10, n_jobs=-1, random_state=42)),
        ]),
        'GradientBoosting': Pipeline([
            ('imp', SimpleImputer(strategy='median')),
            ('clf', GradientBoostingClassifier(n_estimators=300, learning_rate=0.05, max_depth=4, random_state=42)),
        ]),
    }
    # Try LightGBM if available
    try:
        from lightgbm import LGBMClassifier
        models['LightGBM'] = Pipeline([
            ('imp', SimpleImputer(strategy='median')),
            ('clf', LGBMClassifier(n_estimators=500, learning_rate=0.05, num_leaves=31, min_child_samples=20, random_state=42, verbose=-1)),
        ])
    except ImportError:
        pass
    try:
        from xgboost import XGBClassifier
        models['XGBoost'] = Pipeline([
            ('imp', SimpleImputer(strategy='median')),
            ('clf', XGBClassifier(n_estimators=400, learning_rate=0.05, max_depth=4, random_state=42, eval_metric='auc', use_label_encoder=False, verbosity=0)),
        ])
    except ImportError:
        pass

    fold_results = []
    oof = {n: np.zeros(len(Xm)) for n in models}
    for name, mod in models.items():
        for fold, (tr, te) in enumerate(skf.split(Xm, ym)):
            mod.fit(Xm.iloc[tr], ym[tr])
            p = mod.predict_proba(Xm.iloc[te])[:, 1]
            oof[name][te] = p
            auc = roc_auc_score(ym[te], p)
            ap = average_precision_score(ym[te], p)
            br = brier_score_loss(ym[te], p)
            fold_results.append({'outcome': outcome_name, 'model': name, 'fold': fold, 'auc': auc, 'ap': ap, 'brier': br})
        print(f'  {name}: AUC = {np.mean([r["auc"] for r in fold_results if r["model"]==name]):.3f}, AP = {np.mean([r["ap"] for r in fold_results if r["model"]==name]):.3f}')
    results_all.extend(fold_results)

    # Save OOF predictions
    pd.DataFrame({**{n: oof[n] for n in oof}, 'y': ym}, index=Xm.index).to_csv(
        os.path.join(ML_OUT, f'oof_{outcome_name}.csv'), encoding='utf-8-sig')

    # Best model -> SHAP
    best_name = pd.DataFrame(fold_results).query('outcome==@outcome_name').groupby('model')['auc'].mean().idxmax()
    print(f'  Best model: {best_name}')
    best_mod = models[best_name]
    best_mod.fit(Xm, ym)
    # Try shap on tree model
    try:
        import shap
        if 'XGBoost' in best_name or 'LightGBM' in best_name or 'GradientBoosting' in best_name or 'RandomForest' in best_name:
            # Use tree explainer
            tree_clf = best_mod.named_steps['clf']
            X_imp = best_mod.named_steps['imp'].transform(Xm)
            X_imp = pd.DataFrame(X_imp, columns=Xm.columns, index=Xm.index)
            explainer = shap.TreeExplainer(tree_clf)
            shap_values = explainer.shap_values(X_imp)
            # For tree classifiers, shap_values for binary class is 2D array per class
            if isinstance(shap_values, list):
                shap_values = shap_values[1]
            # Mean abs SHAP
            mean_shap = np.abs(shap_values).mean(axis=0)
            shap_imp = pd.DataFrame({'feature': Xm.columns, 'mean_abs_shap': mean_shap}).sort_values('mean_abs_shap', ascending=False)
            shap_imp.to_csv(os.path.join(ML_OUT, f'shap_importance_{outcome_name}.csv'), index=False, encoding='utf-8-sig')
            # Plot
            fig, ax = plt.subplots(figsize=(8, 6))
            top = shap_imp.head(15)
            ax.barh(top['feature'][::-1], top['mean_abs_shap'][::-1], color='#1f77b4')
            ax.set_xlabel('Mean |SHAP value|')
            ax.set_title(f'SHAP feature importance ({best_name})\noutcome: {outcome_name}')
            plt.tight_layout()
            plt.savefig(os.path.join(ML_OUT, f'shap_top15_{outcome_name}.png'), dpi=150, bbox_inches='tight')
            plt.close()
            print(f'  Top 5 SHAP features:')
            print(shap_imp.head(5).to_string(index=False))
    except Exception as e:
        print(f'  SHAP failed: {e}')

# Save results summary
res_df = pd.DataFrame(results_all)
res_df.to_csv(os.path.join(ML_OUT, 'cv_results.csv'), index=False, encoding='utf-8-sig')
summary = res_df.groupby(['outcome', 'model']).agg(auc_mean=('auc','mean'), auc_sd=('auc','std'),
                                                     ap_mean=('ap','mean'), brier_mean=('brier','mean')).round(3)
summary.to_csv(os.path.join(ML_OUT, 'summary.csv'), encoding='utf-8-sig')
print('\n=== ML SUMMARY ===')
print(summary.to_string())
