# -*- coding: utf-8 -*-
"""
Build long-form panel for KCYPS m1 cohort (2018-2024, 7 waves).
Compute key scale scores per wave.
Output: long.parquet, wide.parquet, scale_reliability.csv.
"""
import pandas as pd
import numpy as np
import os, json
import pyreadstat
from sklearn.metrics import pairwise_distances

ROOT = r'D:/2026/SCI/청소년패널조사'
OUT = os.path.join(ROOT, 'output', 'desc')
os.makedirs(OUT, exist_ok=True)

# Define scales (item suffix -> meaning) — use construct dictionaries
SCALES = {
    'life_sat': {       # YPSY1 SWLS
        'items': [f'YPSY1A0{i}' for i in range(1,6)],
        'reverse': [],
        'min_max': (1, 4),  # KCYPS used 1=전혀그렇지않다 4=매우그렇다 (4-pt) — verify
    },
    'happy': {          # YPSY2 Subjective happiness
        'items': [f'YPSY2A0{i}' for i in range(1,5)],
        'reverse': [4],   # 4번 문항은 역코딩 (불행한 사람과 비교)
        'min_max': (1, 7),  # SHS uses 1-7 typically
    },
    'self_esteem': {    # YPSY3 Rosenberg
        'items': [f'YPSY3A{i:02d}' for i in range(1,11)],
        'reverse': [2, 5, 6, 8, 9],   # standard Rosenberg reverse items 2,5,6,8,9
        'min_max': (1, 4),
    },
    'attention': {      # YPSY4A
        'items': [f'YPSY4A0{i}' for i in range(1,8)],
        'reverse': [],
        'min_max': (1, 4),
    },
    'aggression': {     # YPSY4B
        'items': [f'YPSY4B0{i}' for i in range(1,7)],
        'reverse': [],
        'min_max': (1, 4),
    },
    'somatic': {        # YPSY4C
        'items': [f'YPSY4C0{i}' for i in range(1,9)],
        'reverse': [],
        'min_max': (1, 4),
    },
    'withdrawal': {     # YPSY4D 사회적 위축
        'items': [f'YPSY4D0{i}' for i in range(1,6)],
        'reverse': [],
        'min_max': (1, 4),
    },
    'depression': {     # YPSY4E 우울 10문항 (E04 자살사고 포함; sensitivity-only - DO NOT use as primary trajectory due to criterion contamination with SI distal outcome)
        'items': [f'YPSY4E{i:02d}' for i in range(1,11)],
        'reverse': [],
        'min_max': (1, 4),
    },
    'depression9': {    # YPSY4E 우울 9문항 (E04 제외) - PRIMARY trajectory measure to avoid criterion contamination with SI distal outcome
        'items': [f'YPSY4E{i:02d}' for i in [1,2,3,5,6,7,8,9,10]],
        'reverse': [],
        'min_max': (1, 4),
    },
    'grit': {           # YPSY7
        'items': [f'YPSY7A0{i}' for i in range(1,9)],
        'reverse': [1, 3, 5, 6],  # Duckworth Short Grit reverse items (consistency-of-interest items)
        'min_max': (1, 5),
    },
    'parent_warmth': {  # YFAM2A
        'items': [f'YFAM2A0{i}' for i in range(1,5)],
        'reverse': [],
        'min_max': (1, 4),
    },
    'parent_reject': {  # YFAM2B
        'items': [f'YFAM2B0{i}' for i in range(1,5)],
        'reverse': [],
        'min_max': (1, 4),
    },
    'parent_autonomy': { # YFAM2C
        'items': [f'YFAM2C0{i}' for i in range(1,5)],
        'reverse': [],
        'min_max': (1, 4),
    },
    'parent_coerce': {   # YFAM2D 강요
        'items': [f'YFAM2D0{i}' for i in range(1,5)],
        'reverse': [],
        'min_max': (1, 4),
    },
    'parent_structure': {  # YFAM2E
        'items': [f'YFAM2E0{i}' for i in range(1,5)],
        'reverse': [],
        'min_max': (1, 4),
    },
    'parent_chaos': {    # YFAM2F 혼돈/무관심
        'items': [f'YFAM2F0{i}' for i in range(1,5)],
        'reverse': [],
        'min_max': (1, 4),
    },
    'peer_attach': {    # YEDU2 — items 9-13 are negative (반대로 코딩)
        'items': [f'YEDU2A{i:02d}' for i in range(1,14)],
        'reverse': [9, 10, 11, 12, 13],
        'min_max': (1, 4),
    },
    'teacher_rel': {    # YEDU3 — all positive items
        'items': [f'YEDU3A{i:02d}' for i in range(1,15)],
        'reverse': [],
        'min_max': (1, 4),
    },
    'sm_addiction': {   # YMDA1 C series 1-15 — smartphone dependency. Reverse 5,10,15
        'items': [f'YMDA1C{i:02d}' for i in range(1,16)],
        'reverse': [5, 10, 15],
        'min_max': (1, 4),
    },
    # Delinquency: count any-yes (since these are frequencies / 0=never)
    'offline_delin': {  # YDLQ1 - item-level frequency, treat as count of behaviors >=2 (any occurrence)
        'items': [f'YDLQ1A{i:02d}' for i in range(1,16)],
        'reverse': [],
        'min_max': None,
        'mode': 'count_any',  # binarize
    },
    'cyber_delin': {    # YDLQ2 - same
        'items': [f'YDLQ2A{i:02d}' for i in range(1,16)],
        'reverse': [],
        'min_max': None,
        'mode': 'count_any',
    },
}

def load_wave(w):
    path = os.path.join(ROOT, f'KCYPS2018m1Yw{w}.dta')
    df, meta = pyreadstat.read_dta(path)
    return df, meta

def reverse_item(x, lo, hi):
    return (lo + hi) - x

def compute_scale(df_wave, scale_name, items_w, reverse_idx, min_max, mode=None, wave=None):
    """Compute scale score for one wave. items_w are wave-suffixed names (e.g. YPSY1A01w1)."""
    sub = df_wave.reindex(columns=items_w).copy()
    # KCYPS missing codes: negative values, 99 (해당없음)
    sub = sub.mask((sub < 0) | (sub >= 90), other=np.nan)
    if mode == 'count_any':
        # Any occurrence -> 1, never -> 0
        # Items: 1=전혀없다, 2~ = 있다 (frequency)
        # Treat val >=2 as "any" — but if scale is freq days, val>=1 might mean "any" too
        # KCYPS YDLQ uses 1-5 scale: 1=없다, 2=한두번, 3=한달에 1~2번, 4=일주일에 1~2번, 5=거의매일
        any_count = (sub > 1).sum(axis=1)  # number of behaviors with any occurrence
        return any_count, sub.notna().sum(axis=1), pd.Series([np.nan]*len(sub), index=sub.index)  # alpha skipped
    else:
        if reverse_idx and min_max:
            lo, hi = min_max
            for idx in reverse_idx:
                col = items_w[idx-1]
                if col in sub.columns:
                    sub[col] = reverse_item(sub[col], lo, hi)
        # mean of available items (require at least 50% answered)
        n_items = len(items_w)
        n_answered = sub.notna().sum(axis=1)
        ok = n_answered >= max(2, n_items // 2)
        score = sub.mean(axis=1)
        score[~ok] = np.nan
        # Cronbach alpha quick (ignoring missing pairwise)
        try:
            X = sub.dropna()
            if len(X) > 30 and X.shape[1] > 1:
                k = X.shape[1]
                alpha = (k / (k - 1)) * (1 - X.var(axis=0, ddof=1).sum() / X.sum(axis=1).var(ddof=1))
            else:
                alpha = np.nan
        except Exception:
            alpha = np.nan
        return score, n_answered, pd.Series([alpha]*len(sub), index=sub.index)

# Build wide panel: id + per-wave scale scores
print('Building per-wave scale scores ...')
all_long = []
alpha_log = []
for w in range(1, 8):
    df, meta = load_wave(w)
    print(f'  wave {w}: shape={df.shape}')
    ids = df[['ID', 'HID', 'PID']].copy()
    ids['wave'] = w
    out = ids.copy()
    for name, spec in SCALES.items():
        items_w = [f'{it}w{w}' for it in spec['items']]
        # filter to existing
        existing = [c for c in items_w if c in df.columns]
        if not existing:
            out[name] = np.nan
            continue
        # Re-derive mapping in same order
        mapped = [c for c in items_w if c in df.columns]
        score, n_ans, alpha = compute_scale(df, name, mapped, spec.get('reverse', []), spec.get('min_max'), spec.get('mode'), wave=w)
        out[name] = score.values
        out[name + '_n'] = n_ans.values
        # log alpha (one value per wave-scale)
        a = alpha.iloc[0] if len(alpha) else np.nan
        alpha_log.append({'wave': w, 'scale': name, 'alpha': a, 'n_items_used': len(existing)})
    # Add sex, birth year if exists
    for v in ['YGENDER', 'YBRT1A', 'YBRT1B', 'WEIGHTA1', 'WEIGHTA2', 'WEIGHTB1', 'WEIGHTB2']:
        col = next((c for c in df.columns if c.startswith(v)), None)
        if col:
            out[v] = df[col].values
    all_long.append(out)

long_df = pd.concat(all_long, ignore_index=True)
long_df.to_parquet(os.path.join(OUT, 'long_panel.parquet'))
long_df.to_csv(os.path.join(OUT, 'long_panel_preview.csv'), index=False, encoding='utf-8-sig')

alpha_df = pd.DataFrame(alpha_log)
alpha_df.to_csv(os.path.join(OUT, 'scale_alpha_per_wave.csv'), index=False, encoding='utf-8-sig')

# Wide form: pivot per scale by wave
scale_cols = list(SCALES.keys())
wide = long_df.pivot_table(index='ID', columns='wave', values=scale_cols)
wide.columns = [f'{s}_w{w}' for s, w in wide.columns]
wide.to_parquet(os.path.join(OUT, 'wide_panel.parquet'))

# Sample size per wave & attrition
n_per_wave = long_df.groupby('wave').size().to_dict()
print('\nN per wave:', n_per_wave)

# Patterns: number of waves each ID appears
pattern = long_df.groupby('ID')['wave'].apply(lambda s: ''.join('1' if w in set(s) else '0' for w in range(1,8)))
patt_counts = pattern.value_counts().reset_index()
patt_counts.columns = ['pattern', 'n']
patt_counts.to_csv(os.path.join(OUT, 'attrition_pattern.csv'), index=False, encoding='utf-8-sig')
print('Top attrition patterns:')
print(patt_counts.head(15).to_string())

print('\nReliability per wave:')
print(alpha_df.pivot(index='scale', columns='wave', values='alpha').round(3).to_string())

# Descriptive per scale per wave
desc = long_df.groupby('wave')[scale_cols].agg(['mean', 'std', 'count']).round(3)
desc.to_csv(os.path.join(OUT, 'scale_desc_per_wave.csv'), encoding='utf-8-sig')
print('\nMeans (depression, life_sat, self_esteem, parent_warmth):')
print(long_df.groupby('wave')[['depression','life_sat','self_esteem','parent_warmth']].mean().round(3).to_string())
