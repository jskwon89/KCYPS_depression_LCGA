# -*- coding: utf-8 -*-
"""
Parse the codebook xlsx properly (Excel preserves unicode internally).
Output: cleaned CB_m, CB_h, CB_p with construct names.
Then build cross-wave variable inventory using stripped variable roots (drop trailing w{N}).
"""
import pandas as pd
import numpy as np
import os, re, json

ROOT = r'D:/2026/SCI/청소년패널조사'
CB = os.path.join(ROOT, 'KCYPS 2018_Codebook.xlsx')
OUT = os.path.join(ROOT, 'output', 'codebook')

# Read CB sheets - try header rows
for name in ['CB_m', 'CB_h', 'CB_p', 'LAYOUT_m', 'LAYOUT_h', 'LAYOUT_p']:
    raw = pd.read_excel(CB, sheet_name=name, header=None)
    # Find header row by looking for likely header content
    hdr_row = None
    for i in range(min(10, len(raw))):
        cells = raw.iloc[i].astype(str).tolist()
        joined = ' | '.join(cells)
        if any(k in joined for k in ['변수', '문항', '코드', '변인', 'Variable', 'NAME', '항목', '구성']):
            hdr_row = i
            break
    if hdr_row is None:
        hdr_row = 4  # guess
    df = pd.read_excel(CB, sheet_name=name, header=hdr_row)
    df.to_csv(os.path.join(OUT, f'{name}.csv'), index=False, encoding='utf-8-sig')
    print(f'{name}: header_row={hdr_row}, shape after header={df.shape}')
    print('  cols:', df.columns.tolist()[:10])
    print('  first row:', df.iloc[0].tolist()[:10] if len(df) else 'empty')
    print()
