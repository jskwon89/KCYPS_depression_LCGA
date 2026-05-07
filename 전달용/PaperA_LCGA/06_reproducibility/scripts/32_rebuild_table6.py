# -*- coding: utf-8 -*-
"""
Rebuild Table 6 with cluster-bootstrap per-wave HRs (v4).
Old Table 6 had average HR + Schoenfeld; new version reports the time-varying
class × wave HRs with 95% percentile CIs and bootstrap p-values that the manuscript
§3.6 actually cites.

Layout:
  Class | Schoenfeld p | bootstrap interaction p | HR W1 [95% CI] (boot p) | HR W4 [95% CI] (boot p) | HR W7 [95% CI] (boot p)
"""
import pandas as pd
import os

ROOT = r'D:/2026/SCI/청소년패널조사'
P9   = os.path.join(ROOT, 'output', 'paperA_lcga', 'primary_9item')
LD9  = os.path.join(P9, 'long_distal')
TBLS = os.path.join(P9, 'tables')

boot = pd.read_csv(os.path.join(LD9, 'cox_HR_per_wave_with_CI_9item_BOOT.csv'))
sch = pd.read_csv(os.path.join(LD9, 'cox_PH_schoenfeld_test.csv'))
sch = sch.rename(columns={sch.columns[0]: 'covariate'})
boot_int = pd.read_csv(os.path.join(LD9, 'cox_interaction_bootstrap_dist_9item.csv'))

# Bootstrap interaction p-values from distribution
def boot_p(col):
    vals = boot_int[col].values
    p_left = (vals <= 0).mean()
    p_right = (vals >= 0).mean()
    p_two = 2 * min(p_left, p_right)
    return max(p_two, 1 / (len(vals) + 1))

int_p = {'c2_x_t': boot_p('c2_x_t'),
         'c3_x_t': boot_p('c3_x_t'),
         'c4_x_t': boot_p('c4_x_t')}
print('Bootstrap interaction p-values:', {k: round(v, 4) for k, v in int_p.items()})

def fmt(row):
    p_str = '<.001' if row['boot_p'] < 0.001 else f"{row['boot_p']:.3f}"
    return f"{row['HR']:.2f} [{row['HR_lo']:.2f}, {row['HR_hi']:.2f}] (p={p_str})"

# Build new Table 6
class_info = [
    ('Stable-Low (ref)', 'c1', None, None),
    ('Persistent-High vs Stable-Low', 'c2', 'c2_x_t', 'c2'),
    ('Decreasing vs Stable-Low', 'c3', 'c3_x_t', 'c3'),
    ('Late-Increasing vs Stable-Low', 'c4', 'c4_x_t', 'c4'),
]

LABEL_MAP = {'c2': 'Persistent-High', 'c3': 'Decreasing', 'c4': 'Late-Increasing'}

rows = []
for class_label, ckey, ikey, sch_key in class_info:
    if ckey == 'c1':
        rows.append({
            'Class': class_label,
            'HR W1 [95% CI] (boot p)': '1.00 (ref)',
            'HR W4 [95% CI] (boot p)': '1.00 (ref)',
            'HR W7 [95% CI] (boot p)': '1.00 (ref)',
            'Schoenfeld p (PH)': '—',
            'Bootstrap interaction p (rate-of-change)': '—',
        })
        continue
    sub = boot[boot['class_short'] == ckey]
    w1 = sub[sub.wave == 1].iloc[0]
    w4 = sub[sub.wave == 4].iloc[0]
    w7 = sub[sub.wave == 7].iloc[0]
    s = sch[sch['covariate'] == sch_key].iloc[0]
    sch_p = '<.001' if s['p'] < 0.001 else f"{s['p']:.3f}"
    int_p_val = int_p[ikey]
    int_p_str = '<.001' if int_p_val < 0.001 else f"{int_p_val:.3f}"
    rows.append({
        'Class': class_label,
        'HR W1 [95% CI] (boot p)': fmt(w1),
        'HR W4 [95% CI] (boot p)': fmt(w4),
        'HR W7 [95% CI] (boot p)': fmt(w7),
        'Schoenfeld p (PH)': sch_p,
        'Bootstrap interaction p (rate-of-change)': int_p_str,
    })

# Compute the add-one bootstrap p floor from the actual boot-distribution sample size
B_boot = len(boot_int)
floor = 1.0 / (B_boot + 1)
floor_str = f'{floor:.4f}'.rstrip('0').rstrip('.')

# Build the floor note that will appear at the bottom of Table 6 / Table 6b
NOTE_TEXT = (
    f'Bootstrap p add-one floor = 1/(B+1) = 1/{B_boot + 1} = {floor:.5f} '
    f'(rounds to .003 at 3 decimals); planned B=1000 supplementary bootstrap '
    f'would lower the floor to .000999.'
)

t6 = pd.DataFrame(rows)
# Append a Note row carrying the floor caveat (reproducible, not a manual edit).
note_row6 = {col: '' for col in t6.columns}
note_row6[t6.columns[0]] = 'Note'
note_row6[t6.columns[1]] = NOTE_TEXT
t6 = pd.concat([t6, pd.DataFrame([note_row6])], ignore_index=True)
t6.to_csv(os.path.join(TBLS, 'Table6_cox_HR_9item.csv'), index=False, encoding='utf-8-sig')
print('Table 6 saved (with bootstrap-floor note row).')

# Also save full per-wave table as Table 6b for supplement
fmt_all_rows = []
for cls_short, cls_name in LABEL_MAP.items():
    sub = boot[boot['class_short'] == cls_short]
    row = {'Class (vs Stable-Low)': cls_name}
    for w in range(1, 8):
        s = sub[sub.wave == w].iloc[0]
        p_str = '<.001' if s['boot_p'] < 0.001 else f"{s['boot_p']:.3f}"
        row[f'W{w}'] = f"{s['HR']:.2f} [{s['HR_lo']:.2f}, {s['HR_hi']:.2f}], p={p_str}"
    fmt_all_rows.append(row)
t6b = pd.DataFrame(fmt_all_rows)
note_row6b = {col: '' for col in t6b.columns}
note_row6b[t6b.columns[0]] = 'Note'
note_row6b[t6b.columns[1]] = NOTE_TEXT
t6b = pd.concat([t6b, pd.DataFrame([note_row6b])], ignore_index=True)
t6b.to_csv(os.path.join(TBLS, 'Table6b_cox_HR_per_wave_full_9item.csv'),
           index=False, encoding='utf-8-sig')
print('Table 6b saved (with bootstrap-floor note row).')

# Also write a standalone notes file so the floor caveat is preserved if a
# downstream consumer strips trailing rows from the CSV.
NOTES_FILE = os.path.join(TBLS, 'Table6_notes.txt')
with open(NOTES_FILE, 'w', encoding='utf-8') as f:
    f.write('Notes for Table 6 / Table 6b (cluster-bootstrap Cox per-wave HRs)\n')
    f.write('=' * 70 + '\n')
    f.write(f'Bootstrap B = {B_boot}.\n')
    f.write(f'Add-one bootstrap p floor = 1/(B+1) = 1/{B_boot + 1} = {floor:.5f}\n')
    f.write('At 3 decimal places this rounds to .003.\n')
    f.write('Cells reading p=0.003 are at this floor; a planned supplementary '
            'B=1000 bootstrap would lower the floor to .000999.\n')
    f.write('p-values larger than the floor are exact bootstrap two-sided '
            'percentile p-values.\n')
print('Notes file:', NOTES_FILE)
