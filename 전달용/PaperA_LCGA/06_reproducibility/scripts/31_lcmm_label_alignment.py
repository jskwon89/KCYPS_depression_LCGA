# -*- coding: utf-8 -*-
"""
Document explicit mapping between R lcmm raw class labels (1-4) and the
manuscript trajectory labels (C1 Stable-Low, C2 Persistent-High, C3 Decreasing,
C4 Late-Increasing). Save aligned summary file.
"""
import pandas as pd
import numpy as np
import os

ROOT = r'D:/2026/SCI/청소년패널조사'
P9   = os.path.join(ROOT, 'output', 'paperA_lcga', 'primary_9item')
RDIR = os.path.join(P9, 'lcmm_R')

# Read R lcmm K=4 predictions
pred = pd.read_csv(os.path.join(RDIR, 'class_pred_K4_9item.csv'))
print('R lcmm K=4 predicted means by raw class label:')
print(pred.to_string(index=False))

# Read posterior probabilities
post = pd.read_csv(os.path.join(RDIR, 'posterior_K4_9item.csv'))
# class column has the raw R class assignment
class_counts = post['class'].value_counts().sort_index()
print('\nR lcmm raw class proportions (modal):')
for k, n in class_counts.items():
    print(f'  R class {k}: n = {n} ({n/len(post)*100:.2f}%)')

# Compute trajectory shape per class -> identify which is which
# Stable-Low: lowest intercept, near-zero slope
# Persistent-High: highest intercept, peak around mid-waves
# Decreasing: high intercept, declining
# Late-Increasing: low-to-mid intercept, rising
shapes = {}
for col in ['Ypred_class1', 'Ypred_class2', 'Ypred_class3', 'Ypred_class4']:
    raw_class = int(col.replace('Ypred_class', ''))
    means = pred[col].values
    intercept_mid = means[3]  # W4 (mid)
    slope_lin = means[6] - means[0]  # W7 - W1
    is_high = intercept_mid > 2.0
    is_rising = slope_lin > 0.3
    is_falling = slope_lin < -0.3
    if intercept_mid < 1.5:
        label = 'C1 Stable-Low'
    elif is_high and not is_falling:
        label = 'C2 Persistent-High'
    elif is_falling:
        label = 'C3 Decreasing'
    elif is_rising:
        label = 'C4 Late-Increasing'
    else:
        label = '?'
    shapes[raw_class] = {'manuscript_label': label,
                         'intercept_mid_w4': intercept_mid,
                         'mean_w1': means[0], 'mean_w7': means[6],
                         'slope_w1_to_w7': slope_lin}

# Build alignment table
align_rows = []
for r_class in [1, 2, 3, 4]:
    info = shapes[r_class]
    n = class_counts.get(r_class, 0)
    align_rows.append({
        'R_lcmm_class': r_class,
        'manuscript_label': info['manuscript_label'],
        'modal_n_R': n,
        'modal_pct_R': round(n / len(post) * 100, 2),
        'mean_w1': round(info['mean_w1'], 3),
        'mean_w4_intercept': round(info['intercept_mid_w4'], 3),
        'mean_w7': round(info['mean_w7'], 3),
        'slope_w1_w7': round(info['slope_w1_to_w7'], 3),
    })
align_df = pd.DataFrame(align_rows)

# Sort by manuscript label (C1, C2, C3, C4)
order = {'C1 Stable-Low': 1, 'C2 Persistent-High': 2,
         'C3 Decreasing': 3, 'C4 Late-Increasing': 4}
align_df['_sort'] = align_df['manuscript_label'].map(order)
align_df = align_df.sort_values('_sort').drop(columns='_sort')

print('\n=== R lcmm <-> manuscript class label alignment ===')
print(align_df.to_string(index=False))
align_df.to_csv(os.path.join(RDIR, 'class_label_alignment_9item.csv'),
                index=False, encoding='utf-8-sig')

# Compare to Python EM modal counts
py_class = pd.read_csv(os.path.join(P9, 'class_assignments_K4_9item.csv'))
py_counts = py_class['class_K4_9item'].value_counts().sort_index()
PY_LABEL = {1: 'C1 Stable-Low', 2: 'C2 Persistent-High',
            3: 'C3 Decreasing', 4: 'C4 Late-Increasing'}
print('\n=== Side-by-side comparison: Python EM vs R lcmm (manuscript label aligned) ===')
comp_rows = []
for ms_label in ['C1 Stable-Low', 'C2 Persistent-High', 'C3 Decreasing', 'C4 Late-Increasing']:
    py_k = next((k for k, v in PY_LABEL.items() if v == ms_label), None)
    r_row = align_df[align_df.manuscript_label == ms_label].iloc[0]
    py_n = py_counts.get(py_k, 0)
    py_pct = py_n / py_counts.sum() * 100
    comp_rows.append({
        'manuscript_label': ms_label,
        'Python_EM_class': py_k, 'Python_modal_n': py_n, 'Python_pct': round(py_pct, 2),
        'R_lcmm_class': int(r_row['R_lcmm_class']),
        'R_modal_n': int(r_row['modal_n_R']), 'R_pct': r_row['modal_pct_R'],
    })
comp = pd.DataFrame(comp_rows)
print(comp.to_string(index=False))
comp.to_csv(os.path.join(RDIR, 'python_vs_R_class_comparison_9item.csv'),
            index=False, encoding='utf-8-sig')

# Generate aligned class predictions (R lcmm with manuscript class labels)
pred_aligned = pred.copy()
rename_map = {f'Ypred_class{r}': f'Ypred_{shapes[r]["manuscript_label"].replace(" ", "_").replace("-", "")}'
              for r in [1, 2, 3, 4]}
pred_aligned = pred_aligned.rename(columns=rename_map)
# Reorder
preferred_order = ['wave',
                   'Ypred_C1_StableLow', 'Ypred_C2_PersistentHigh',
                   'Ypred_C3_Decreasing', 'Ypred_C4_LateIncreasing']
pred_aligned = pred_aligned[[c for c in preferred_order if c in pred_aligned.columns]]
pred_aligned.to_csv(os.path.join(RDIR, 'class_pred_K4_9item_LABELED.csv'),
                    index=False, encoding='utf-8-sig')
print('\nLabeled prediction saved to class_pred_K4_9item_LABELED.csv')

print('\nDone.')
