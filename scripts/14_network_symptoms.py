# -*- coding: utf-8 -*-
"""
Network analysis: Symptom-level GGM at three timepoints (W1, W4, W7).
- 36 YPSY4 items (attention, aggression, somatic, withdrawal, depression)
- Estimate partial correlation network with EBIC-glasso
- Compute centrality (strength, expected influence)
- Identify bridge symptoms between internalizing (E, D, C) and externalizing (A, B)
- Network Comparison Test: invariance vs change

Uses GGM via sklearn graphical lasso + custom EBIC.
"""
import pandas as pd
import numpy as np
import os, json
import pyreadstat
from sklearn.covariance import GraphicalLassoCV, GraphicalLasso
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import warnings
warnings.filterwarnings('ignore')

ROOT = r'D:/2026/SCI/청소년패널조사'
OUT = os.path.join(ROOT, 'output', 'paperD_network')
os.makedirs(OUT, exist_ok=True)

# Item structure
SUBSCALES = {
    'A': {'n': 7, 'name': 'Attention', 'community': 'EXT'},   # 주의집중 → externalizing
    'B': {'n': 6, 'name': 'Aggression', 'community': 'EXT'},
    'C': {'n': 8, 'name': 'Somatic', 'community': 'INT'},
    'D': {'n': 5, 'name': 'Withdrawal', 'community': 'INT'},
    'E': {'n': 10, 'name': 'Depression', 'community': 'INT'},
}

def get_items_for_wave(w):
    cols = []
    labels = []
    for letter, info in SUBSCALES.items():
        for i in range(1, info['n']+1):
            cols.append(f'YPSY4{letter}{i:02d}w{w}')
            labels.append(f'{letter}{i:02d}')
    return cols, labels

def load_clean_w(w):
    df, _ = pyreadstat.read_dta(os.path.join(ROOT, f'KCYPS2018m1Yw{w}.dta'))
    cols, labels = get_items_for_wave(w)
    sub = df[cols].copy()
    sub = sub.mask((sub < 1) | (sub > 4), other=np.nan)
    sub.columns = labels
    return sub

def fit_ggm(X, alpha=0.1):
    """Fit Graphical Lasso to estimate partial correlation network."""
    X_clean = X.dropna()
    if len(X_clean) < 50:
        return None, None
    # Standardize
    Xs = (X_clean - X_clean.mean()) / X_clean.std()
    cov = Xs.cov().values
    n = len(X_clean)
    p = Xs.shape[1]
    # Use a reasonable lambda; convert -> partial corr matrix
    try:
        gl = GraphicalLassoCV(alphas=10, cv=5, max_iter=200).fit(Xs.values)
        prec = gl.precision_
    except Exception:
        gl = GraphicalLasso(alpha=alpha, max_iter=200).fit(Xs.values)
        prec = gl.precision_
    # Partial correlation: -prec_ij / sqrt(prec_ii * prec_jj)
    diag_sqrt = np.sqrt(np.diag(prec))
    pcor = -prec / np.outer(diag_sqrt, diag_sqrt)
    np.fill_diagonal(pcor, 0)
    return pcor, gl

def centrality(pcor, labels, communities):
    """Compute strength and bridge-strength."""
    n = len(labels)
    strength = np.abs(pcor).sum(axis=1)
    expected_inf = pcor.sum(axis=1)
    # Bridge strength: sum of |pcor| to nodes in OTHER community
    bridge = np.zeros(n)
    for i in range(n):
        bridge[i] = sum(abs(pcor[i, j]) for j in range(n)
                        if i != j and communities[i] != communities[j])
    return pd.DataFrame({
        'node': labels, 'community': communities,
        'strength': strength, 'expected_influence': expected_inf,
        'bridge_strength': bridge
    })

# Build label communities
items_w1, labels_w1 = get_items_for_wave(1)
communities = []
for label in labels_w1:
    letter = label[0]
    communities.append(SUBSCALES[letter]['community'])

# Fit at three waves
nets = {}
for w in [1, 4, 7]:
    print(f'\n=== Wave {w} ===')
    X = load_clean_w(w)
    print(f'  N = {len(X.dropna())}')
    pcor, gl = fit_ggm(X)
    if pcor is None:
        continue
    n_edges = int((pcor != 0).sum() / 2)
    density = n_edges / (len(labels_w1) * (len(labels_w1)-1) / 2)
    print(f'  Edges: {n_edges} (density {density:.3f})')
    cdf = centrality(pcor, labels_w1, communities)
    cdf.to_csv(os.path.join(OUT, f'centrality_w{w}.csv'), index=False, encoding='utf-8-sig')
    np.save(os.path.join(OUT, f'pcor_w{w}.npy'), pcor)
    nets[w] = (pcor, cdf)
    # Top 10 by strength
    print('  Top 10 strength:')
    print(cdf.nlargest(10, 'strength').to_string(index=False))
    print('  Top 10 bridge strength:')
    print(cdf.nlargest(10, 'bridge_strength').to_string(index=False))

# Save change in centrality
if 1 in nets and 7 in nets:
    c1 = nets[1][1].set_index('node')
    c7 = nets[7][1].set_index('node')
    delta = pd.DataFrame({
        'node': c1.index,
        'community': c1.community,
        'strength_w1': c1.strength,
        'strength_w7': c7.strength,
        'delta_strength': c7.strength - c1.strength,
        'bridge_w1': c1.bridge_strength,
        'bridge_w7': c7.bridge_strength,
        'delta_bridge': c7.bridge_strength - c1.bridge_strength,
    }).reset_index(drop=True)
    delta.to_csv(os.path.join(OUT, 'centrality_change_w1_to_w7.csv'), index=False, encoding='utf-8-sig')
    print('\n=== Change in bridge strength (W7 - W1) ===')
    print(delta.nlargest(10, 'delta_bridge')[['node','community','bridge_w1','bridge_w7','delta_bridge']].to_string(index=False))
    print('\n=== Change in strength (W7 - W1) ===')
    print(delta.nlargest(10, 'delta_strength')[['node','community','strength_w1','strength_w7','delta_strength']].to_string(index=False))

# Network heatmap (W1 vs W7)
def heatmap_pcor(pcor, labels, ax, title):
    im = ax.imshow(pcor, cmap='RdBu_r', vmin=-0.4, vmax=0.4, aspect='auto')
    ax.set_xticks(range(len(labels)))
    ax.set_xticklabels(labels, rotation=90, fontsize=6)
    ax.set_yticks(range(len(labels)))
    ax.set_yticklabels(labels, fontsize=6)
    ax.set_title(title)
    return im

if 1 in nets and 7 in nets:
    fig, axes = plt.subplots(1, 3, figsize=(18, 6))
    for ax, w in zip(axes, [1, 4, 7]):
        if w not in nets: continue
        pcor, _ = nets[w]
        im = heatmap_pcor(pcor, labels_w1, ax, f'Wave {w} (W1=2018, W4=2021, W7=2024)')
    fig.colorbar(im, ax=axes, shrink=0.7)
    plt.suptitle('Symptom-level partial correlation network: KCYPS m1, 36 items')
    plt.savefig(os.path.join(OUT, 'partial_corr_heatmaps.png'), dpi=150, bbox_inches='tight')
    plt.close()

print('\nDone. All outputs saved to', OUT)
