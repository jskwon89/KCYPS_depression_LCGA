# -*- coding: utf-8 -*-
"""
KCYPS 2018 m1 cohort (waves 1-7, 2018-2024) -- explore codebook & data.
Output: variable inventory, wave coverage matrix, candidate constructs.
"""
import pandas as pd
import numpy as np
import json
import os

ROOT = r'D:/2026/SCI/청소년패널조사'
CB = os.path.join(ROOT, 'KCYPS 2018_Codebook.xlsx')
OUT = os.path.join(ROOT, 'output', 'codebook')
os.makedirs(OUT, exist_ok=True)

# Read all relevant sheets with header detection
sheets = ['LAYOUT_m', 'LAYOUT_h', 'LAYOUT_p', 'CB_m', 'CB_h', 'CB_p']
xl = pd.ExcelFile(CB)
print('Available sheets:', xl.sheet_names)

for s in sheets:
    df = pd.read_excel(CB, sheet_name=s, header=None)
    df.to_csv(os.path.join(OUT, f'{s}_raw.csv'), index=False, encoding='utf-8-sig')
    print(f'{s}: shape={df.shape}')

# Try to auto-detect header row in LAYOUT/CB sheets
for s in sheets:
    df = pd.read_excel(CB, sheet_name=s, header=None)
    # Find first row with multiple non-NaN string entries that look like headers
    for i in range(min(20, len(df))):
        row = df.iloc[i].dropna().astype(str).tolist()
        if any(any(k in r.lower() for k in ['var', 'variable', 'name', 'label', '변수', '문항', '코드']) for r in row):
            print(f'\n--- {s} row {i} candidate header: {row[:8]}')
            break
